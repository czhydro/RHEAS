"""
Code to plot and analyze RHEAS results
------------------------------------------------------------

Plot results and statistics for
 - VIC
     Annual Water balance for the watershed
     Monthly time series plot watershed scale
     Streamflow results for couple of locations w/ stats
     
 - DSSAT
     Yield boxplot by watershed compared with NASS-ish 
     Yield map 
     Max Waterstress per season plots
     Irrigation demand per season per acre
     
@author: Vikalp Mishra
         Sept - 2020
         NASA SERVIR | SPoRT | ESSC-UAH
         
         Updates
         Added VIC streamflow results
         Mar 2021
         
         Added DSSAT results plots
         Oct 2023
         
         Latest code check and updates (libs)
         Aug 2025
         
--------------------------------------------------------------
"""


#import matplotlib.pyplot as plt
#import glob, os, sys
#import seaborn as sns
#import dataretrieval as usgs
#import pandas as pd
#import numpy as np
#import geopandas as gpd
#import rasterio as rio
#from datetime import date, datetime, timedelta
#from matplotlib import colors
#from rasterio.io import MemoryFile
#import requests


def fetch_nass_yields_by_state(state_alpha, year):
    
    import requests
    import pandas as pd
    
    
    url = "https://quickstats.nass.usda.gov/api/api_GET/"
    params = {
        "key": "ED0239DB-2C3E-39CD-A789-D7C19631774B",
        "commodity_desc": "CORN",
        "statisticcat_desc": "YIELD",
        "unit_desc": "BU / ACRE",
        "source_desc": "SURVEY",
        "sector_desc": "CROPS",
        "state_alpha": state_alpha,
        "agg_level_desc": "COUNTY",
        "year": str(year),
        "format": "JSON"
    }

    response = requests.get(url, params=params)
    response.raise_for_status()
    data = response.json()

    # Convert to DataFrame and clean
    df = pd.DataFrame(data["data"])
    df = df[["state_fips_code", "county_code", "county_name", "Value"]].copy()
    df["FIPS"] = df["state_fips_code"].str.zfill(2) + df["county_code"].str.zfill(3)
    df["Yield_kg_ha"] = pd.to_numeric(df["Value"], errors="coerce") / 0.0159
    df = df[["FIPS", "Yield_kg_ha"]]
    return df



def get_fips_from_latlon(lat, lon):
    import requests
    
    
    url = f"https://geo.fcc.gov/api/census/block/find?latitude={lat}&longitude={lon}&format=json"
    response = requests.get(url).json()
    fips = response['County']['FIPS']
    county = response['County']['name']
    state = response['State']['code']
    return fips, county, state


def rast_read(infile):
    import rasterio as rio
    import numpy as np    
    import warnings
        
    warnings.filterwarnings("ignore")
        
    with rio.open(infile) as src:
        d = src.read(1)
        d[d<0] = np.nan
        return d

# function to compute daily stats over all years
def calc_stats(flist, outf):
    import rasterio as rio
    import numpy as np
    import warnings
    
    warnings.filterwarnings("ignore")
    
    var = 'x'
    with rio.open(flist[0]) as src:
        profile = src.profile
        
    array_list = [rast_read(x) for x in flist]

    for k in range(0,len(array_list)):
        array_list[k][array_list[k]<0.0] = np.nan
    if 'SM' in var: 
        array_outm = np.nanmean(array_list, axis=0)
    else:
        array_outm = np.nansum(array_list,axis=0)
    
    with rio.open(outf,'w',**profile) as dst:
        dst.write(array_outm,1)
    return 
    

def res_analysis(self, vic=False,dssat=False, rout = False):
    
    import matplotlib.pyplot as plt
    import glob, os, sys
    import pandas as pd
    import numpy as np
    import geopandas as gpd
    import rasterio as rio
    from datetime import date, datetime, timedelta
    from matplotlib import colors
    from rasterio.io import MemoryFile
    import requests
    import warnings
    
    warnings.filterwarnings("ignore", message=".*CPLE_AppDefined.*")

    
    if vic ==None and dssat == None:
        print('None of the models selected for results analysis.')
        print('Ending result analysis section....')
        sys.exit()
    elif vic == True:
        res_dir = self.dbpath 
        
        # water balance
        bf_files = glob.glob(res_dir+'/baseflow*.tif')
        if len(bf_files)>0:
            
            startyear = self.startyear
            startmonth = self.startmonth
            startday = self.startday 
            startdate = datetime(startyear, startmonth, startday)
            endyear = self.endyear 
            endmonth = self.endmonth 
            endday = self.endday 
            enddate = datetime(endyear, endmonth, endday)
            tmppath = self.model_path
            
            if not os.path.exists(res_dir+'/PLOTS/'):
                os.makedirs(res_dir+'/PLOTS/')
            
            bf = []
            rf = []
            ppt = []
            et = []
            st = []
            sm = []
            dt = []
            yrr = []
            
            dates = pd.date_range(startdate, enddate).date
            for yr in range(startyear,endyear+1):
                sdate = datetime(yr,1,1)
                edate = datetime(yr,12,31)
                str_sdate = datetime.strftime(sdate,'%Y%m%d')
                str_edate = datetime.strftime(edate,'%Y%m%d')
                
                if enddate < datetime(yr,12,31):
                    str_edate = datetime.strftime(enddate,'%Y%m%d')
                
                sm11 = rast_read(res_dir+'/soil_moist_'+str_sdate+'_01_00.tif')
                sm12 = rast_read(res_dir+'/soil_moist_'+str_sdate+'_02_00.tif')
                sm13 = rast_read(res_dir+'/soil_moist_'+str_sdate+'_03_00.tif')
                sm21 = rast_read(res_dir+'/soil_moist_'+str_edate+'_01_00.tif')
                sm22 = rast_read(res_dir+'/soil_moist_'+str_edate+'_02_00.tif')
                sm23 = rast_read(res_dir+'/soil_moist_'+str_edate+'_03_00.tif')
                
                sm_init = sm11+sm12+sm13
                sm_fin = sm21+sm22+sm23
                
                st_tot = np.nanmean(sm_fin - sm_init)
                yrr.append(str_edate)
                st.append(st_tot)
                
            df_ann = pd.DataFrame({'Year':yrr,'dS':st})
            df_ann = df_ann.set_index('Year')
            df_ann.index = pd.to_datetime(df_ann.index)
                
            for dts in dates:
                cdate = datetime.strftime(dts,'%Y%m%d')
                
                bf_tmp = rast_read(res_dir+'/baseflow_'+str(cdate)+'_01_00.tif')
                rf_tmp = rast_read(res_dir+'/runoff_'+str(cdate)+'_01_00.tif')
                ppt_tmp = rast_read(res_dir+'/rainf_'+str(cdate)+'_01_00.tif')
                et_tmp = rast_read(res_dir+'/evap_'+str(cdate)+'_01_00.tif')
                sm1_tmp= rast_read(res_dir+'/soil_moist_'+str(cdate)+'_01_00.tif')
                sm2_tmp= rast_read(res_dir+'/soil_moist_'+str(cdate)+'_02_00.tif')
                sm3_tmp= rast_read(res_dir+'/soil_moist_'+str(cdate)+'_03_00.tif')

                sm_tmp = sm1_tmp+sm2_tmp+sm3_tmp
                
                bf.append(np.nanmean(bf_tmp))
                rf.append(np.nanmean(rf_tmp))
                et.append(np.nanmean(et_tmp))
                ppt.append(np.nanmean(ppt_tmp))
                sm.append(np.nanmean(sm_tmp))
                dt.append(cdate)
                
            df_daily = pd.DataFrame({"Date": dt,"P": ppt,"ET": et,
                                     "Runoff": rf, "Baseflow": bf, "SM":sm})
            
            df_daily = df_daily.sort_values("Date").reset_index(drop=True)
            df_daily['Date'] = pd.to_datetime(df_daily.Date)
            df_daily = df_daily.set_index('Date')
            
            df_mon = df_daily.resample('M').sum()
            df_mon.index = pd.to_datetime(df_mon.index)
                    
            # Derived columns
            df_ann_tmp = df_daily.resample('Y').sum()
            df_ann = pd.merge(df_ann,df_ann_tmp,left_index=True, right_index=True)
            
            df_ann["Q"] = df_ann["Runoff"] + df_ann["Baseflow"]
            df_ann["Residual"] = df_ann["P"] - (df_ann["ET"] + df_ann["Q"] + df_ann["dS"])
            df_ann["Residual_pctP"] = 100 * df_ann["Residual"] / df_ann["P"]
            
            
            ######################
            #  FIGURE SECTION    #
            ######################
            fig, ax = plt.subplots(figsize=(10,7))
    
            x = np.arange(len(df_ann))
            width = 0.35  # bar width
            # positions: precip bar on left, stacked bar on right
            xp = x - width/2
            xs = x + width/2
                  
            # Precipitation bar
            bp = ax.bar(xp, df_ann['P'].values, width, label="Precip", 
                        color='lightgray', hatch='/',edgecolor = 'k')
        
            # Stacked bar: ET (bottom), Q (middle), ΔS (top)
            b_et = ax.bar(xs, df_ann['ET'].values, width, label="ET",
                          color = 'cornflowerblue',edgecolor = 'k')
            b_q  = ax.bar(xs, df_ann["Q"].values, width, bottom=df_ann['ET'].values, 
                          label="Q (Runoff+Baseflow)", color = 'mediumturquoise',edgecolor = 'k')
            b_ds = ax.bar(xs, df_ann['dS'].values, width, 
                          bottom=df_ann['ET'].values + df_ann["Q"].values, 
                          label="ΔS", color = 'khaki',edgecolor = 'k')
        
            # Secondary axis: residual percent of P
            ax2 = ax.twinx()
            ax2.plot(x, df_ann["Residual_pctP"].values, marker="o", linewidth=1.5, label="Residual (% of P)", color = 'k')
            ax2.axhline(0, linestyle="--", linewidth=1, alpha=0.5, c= 'k')
        
            units = "mm/yr"
            title = "Annual Water Balance"
        
            ax.set_title(title, fontsize=14)
            ax.set_ylabel(units, fontsize=14)
            ax2.set_ylabel("Residual (% of P)", fontsize=14)
            ax.set_xticks(x)
            ax.set_xticklabels(df_ann.index.year.astype(str).tolist())
        
            
            ax.legend()
            plt.savefig(res_dir+'/PLOTS/WatBal_Plot.png',dpi=400)
            plt.close()
    

            ######################################################################################
            # PRECIPITATION HEAT MAP + ANOMALY
            # -----------------------------------------------------------------------------------
                
            precip = df_daily['P'].resample('M').sum()
            
            # Create Year/Month table
            mat = precip.to_frame().assign(
                Year=precip.index.year,
                Month=precip.index.month
            ).pivot(index='Year', columns='Month', values='P')
            
            clim = mat.mean(axis=0)
            anom = mat - clim
            
            ######################
            #  FIGURE SECTION    #
            ######################
            fig, axes = plt.subplots(1, 2, figsize=(14, 4))
            
            # Heatmap of values
            im0 = axes[0].imshow(mat, aspect='auto', origin='lower', cmap='Blues')
            axes[0].set_title('Monthly Precipitation (mm)')
            axes[0].set_xlabel('Month')
            axes[0].set_ylabel('Year')
            axes[0].set_xticks(range(12))
            axes[0].set_xticklabels(range(1, 13))
            axes[0].set_yticks(range(len(mat.index)))
            axes[0].set_yticklabels(mat.index)
            plt.colorbar(im0, ax=axes[0], label='mm')
            
            # Heatmap of anomalies
            vmax = np.nanmax(np.abs(anom.values))
            im1 = axes[1].imshow(anom, aspect='auto', origin='lower', cmap='RdBu',
                                 vmin=-vmax, vmax=vmax)
            axes[1].set_title('Precipitation Anomaly (mm)')
            axes[1].set_xlabel('Month')
            axes[1].set_ylabel('Year')
            axes[1].set_xticks(range(12))
            axes[1].set_xticklabels(range(1, 13))
            axes[1].set_yticks(range(len(anom.index)))
            axes[1].set_yticklabels(anom.index)
            plt.colorbar(im1, ax=axes[1], label='mm')  
            plt.savefig(res_dir+'/PLOTS/Precip_HeatMap.png',dpi=400)
            plt.close()
            
            
            
            ######################################################################################
            # P VS ET PLOT - BY MONTH AND YEAR
            # -----------------------------------------------------------------------------------
            p_col = 'P'
            et_col = 'ET' 

            x = df_mon[p_col].to_numpy()
            y = df_mon[et_col].to_numpy()
            m = df_mon.index.month
            yrs = df_mon.index.year
            mask = np.isfinite(x) & np.isfinite(y)
        
            x, y, m, yrs = x[mask], y[mask], m[mask], yrs[mask]
        
            xymax = np.nanmax([np.nanmax(x), np.nanmax(y)])
            lim = xymax * 1.05 if np.isfinite(xymax) else 1.0
            
            ######################
            #  FIGURE SECTION    #
            ######################
            fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)
        
            # --- Left: color by month ---
            sc1 = axes[0].scatter(x, y, c=m, cmap="twilight", s=90, edgecolor="k", alpha=0.85)
            axes[0].plot([0, lim], [0, lim], "--", color="gray", alpha=0.7)
            axes[0].set_title("P vs ET (Month)",fontsize =14)
            axes[0].set_xlabel("P (mm/month)",fontsize =14)
            axes[0].set_ylabel("ET (mm/month)",fontsize =14)
            axes[0].set_xlim(0, lim); axes[0].set_ylim(0, lim)
            cbar1 = plt.colorbar(sc1, ax=axes[0])
            cbar1.set_label("Month", fontsize =14)
            cbar1.set_ticks(np.arange(1, 13))
            cbar1.set_ticklabels(["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"])
        

            uniq_years = np.unique(yrs)
            cmap_years = plt.get_cmap("Set2", len(uniq_years))  # discrete colormap with N colors
            norm_years = colors.BoundaryNorm(
                boundaries=np.arange(len(uniq_years)+1)-0.5,
                ncolors=len(uniq_years))
            
            sc2 = axes[1].scatter(
                x, y,
                c=[np.where(uniq_years == yr)[0][0] for yr in yrs],  # index of year in uniq_years
                cmap=cmap_years,
                norm=norm_years,
                s=90, edgecolor="k", alpha=0.85)
            
            axes[1].plot([0, lim], [0, lim], "--", color="gray", alpha=0.7)
            axes[1].set_title("P vs ET — (Year)",fontsize =14)
            axes[1].set_xlabel("P (mm/month)",fontsize =14)
            axes[1].set_ylabel("ET (mm/month)",fontsize =14)
            axes[1].set_xlim(0, lim)
            axes[1].set_ylim(0, lim)
            
            cbar2 = plt.colorbar(sc2, ax=axes[1], ticks=np.arange(len(uniq_years)))
            cbar2.set_label("Year",fontsize =14)
            cbar2.set_ticklabels(uniq_years.astype(str))
            plt.savefig(res_dir+'/PLOTS/P_ET_plot.png',dpi=400)
            plt.close()
            
                    
        
            ######################################################################################
            # RASTER PLOTS FOR THE MONTHS OF MAY AND SEPT
            # -----------------------------------------------------------------------------------
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature
            import geopandas as gpd
            
            gdf = gpd.read_file(self.basin)
            
            for yr in range(startyear, endyear+1):
                mons = ['05','09']
                for mon in mons:
                    yyyymm = str(yr)+mon
                    
                    #vars = ['baseflow','rainf','runoff','evap']
                    bf_files = glob.glob(res_dir+'/baseflow_'+yyyymm+'*.tif')
                    rf_files = glob.glob(res_dir+'/runoff_'+yyyymm+'*.tif')
                    et_files = glob.glob(res_dir+'/evap_'+yyyymm+'*.tif')
                    ppt_files = glob.glob(res_dir+'/rainf_'+yyyymm+'*.tif')
                    
                    tiff_files = [res_dir+'/temp_bf_'+yyyymm+'.tif',res_dir+'/temp_rf_'+yyyymm+'.tif',
                                  res_dir+'/temp_ppt_'+yyyymm+'.tif', res_dir+'/temp_et_'+yyyymm+'.tif']
                    
                    calc_stats(bf_files, tiff_files[0])
                    calc_stats(rf_files, tiff_files[1])
                    calc_stats(ppt_files, tiff_files[2])
                    calc_stats(et_files, tiff_files[3])
                   
        
                    variable_specs = [{"name": "Baseflow","cmap": "Blues","vmin": 0, "vmax": 80, "units": "mm"},
                        {"name": "Runoff", "cmap": "YlGnBu", "vmin": 0, "vmax": 60, "units": "mm"},
                        {"name": "Precip", "cmap": "plasma_r","vmin": 0, "vmax": 400, "units": "mm"},
                        {"name": "Evap", "cmap": "PuBuGn", "vmin": 0, "vmax": 200, "units": "mm"},]
                    
                    ######################
                    #  FIGURE SECTION    #
                    ######################
                    fig, axes = plt.subplots(2, 2,figsize=(12, 13),
                        subplot_kw={'projection': ccrs.PlateCarree()})
                    axes = axes.flatten()
                    fig.tight_layout(rect=[0.03, 0.03, 0.98, 0.97])
                    
                    for ax, tiff_file, spec in zip(axes, tiff_files, variable_specs):
                        with rio.open(tiff_file) as src:
                            img = src.read(1).astype(float)
                            img[img <= 0] = np.nan  # mask non-positive if that's your rule
                            extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
                            vmin = 0
                            vmax = np.nanpercentile(img, 0.95)
                            
                        # If no fixed vmin/vmax, compute robust limits (2–98th percentiles)
                        if spec["vmin"] is None or spec["vmax"] is None:
                            finite_vals = img[np.isfinite(img)]
                            if finite_vals.size > 0:
                                vmin = np.nanpercentile(finite_vals, 2) if spec["vmin"] is None else spec["vmin"]
                                vmax = np.nanpercentile(finite_vals, 98) if spec["vmax"] is None else spec["vmax"]
                            else:
                                vmin, vmax = 0, 1
                        else:
                            vmin, vmax = spec["vmin"], spec["vmax"]
                            
                        # Raster
                        im = ax.imshow(img, cmap=spec["cmap"], extent=extent, 
                                       transform=ccrs.PlateCarree(),vmin=vmin, 
                                       vmax=vmax, alpha=0.9, zorder=1)
                    
                        # Base map layers (land/backdrop)
                        ax.add_feature(cfeature.LAND, color="white", zorder=0)
                        ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=0.6)
                        ax.coastlines(linewidth=0.7)
                    
                        # Country boundaries from your GeoDataFrame (assumes PlateCarree)
                        gdf.plot(ax=ax, facecolor="none", edgecolor="black", linewidth=1.0, transform=ccrs.PlateCarree())
                                       
                        ax.set_title(f'{spec["name"]}', fontsize=14)
                        ax.set_xticks([])
                        ax.set_yticks([])
                    
                        cbar = fig.colorbar(im, ax=ax, orientation="vertical", pad=0.03, fraction=0.046, extend="max")
                        cbar.set_label(spec["units"], fontsize=12)
                    
                    if mon =='05':
                        cmon = 'May' 
                    else: cmon = 'Sept'
                    
                    fig.suptitle(f'VIC Hydrologic Outputs for {cmon}, {yr}', fontsize = 16)
                    # Tighten layout
                    plt.subplots_adjust(hspace=0.1, wspace=0.14, 
                                        top=0.9,bottom=0.042,left=0.043,right=0.937)
                    plt.savefig(f'{res_dir}/PLOTS/VIC_Map_{cmon}, {yr}.png',dpi=400)
                    plt.close()
                    
                    for f in tiff_files:
                        os.remove(f)
            
            
            
            
            
            ######################################################################################
            # PLOT STREAMFLOW RESULTS 
            # -----------------------------------------------------------------------------------
            import dataretrieval.nwis as nwis
            
            tmp_dir = tmppath+'/rout/'
            glist_file =  tmp_dir+'/GaugeList.txt'
            if os.path.exists(glist_file):
                glist = pd.read_table(glist_file, header = None, sep = '\s+')
                for i in range(len(glist)):
                    point = glist[0][i]
                    gauge = glist[1][i]
                    
                    sflow_file = tmp_dir+point+'.day'
                    if os.path.exists(sflow_file):
                        sf = pd.read_table(sflow_file, header=None,sep= '\s+')
                        dd = []
                        for i in range(0,len(sf)):
                            row = sf.iloc[i]
                            yr = int(row[0])
                            mm = int(row[1])
                            dy = int(row[2])
                            tmp = datetime(yr,mm,dy)
                            tmp = datetime.strftime(tmp, '%Y%m%d')
                            dd.append(tmp)
                            
                        sf['Date'] = dd
                        sf = sf.drop(columns= [0,1,2])
                        sf = sf.set_index('Date')
                        sf.index = pd.to_datetime(sf.index)
                        
                        mindate = sf.index.min()
                        maxdate = sf.index.max()
                        
                        #getting USGS gauge data 
                        site = '0'+str(gauge)
                        df = nwis.get_record(sites=site, service='dv', start=mindate.strftime('%Y-%m-%d'), 
                                         end=maxdate.strftime('%Y-%m-%d'),parameterCd= '00060')
                        if len(df)>0:
                            if '00060_Mean' in df.columns:
                                obs = df['00060_Mean']
                            if '00060_2' in df.columns:
                                obs = df['00060_2_Mean']
                                                    
                            obs.index = pd.to_datetime(obs.index.date)
                            obs.index = obs.index.date
                            obs = obs
                        
                        sf = sf.merge(obs, left_index=True, right_index=True, how='right')
                        sf.columns = ['Sim','Obs']
                        
                        ndf = sf.copy()
                        ndf = ndf.dropna()
                        r = np.round(ndf.corr().values[0,1],2)
                        denom = np.sum((ndf.Obs-ndf.Obs.mean())**2)
                        if denom>0:
                            NSE = 1 - np.sum((ndf.Sim - ndf.Obs)**2)/denom
                        else: NSE = np.nan
                        
                        alpha = (sf.Sim.std() / sf.Obs.std()) if sf.Obs.std() != 0 else np.nan
                        beta  = (sf.Sim.mean() / sf.Obs.mean()) if sf.Obs.mean() != 0 else np.nan
                        if np.any(np.isnan([r, alpha, beta])):
                            KGE = np.nan
                        else:
                            KGE = 1.0 - np.sqrt((r - 1.0)**2 + (alpha - 1.0)**2 + (beta - 1.0)**2)
                            
                        KGE = np.round(KGE,3)
                        NSE = np.round(NSE,3)
                        bias = np.round(sf.Sim.mean() - sf.Obs.mean(),2)
                        RMSE = np.round(np.sqrt(((sf.Sim - sf.Obs)**2).mean()),2)
                        
                        textstr = (
                                f"n     = {len(ndf)}\n"
                                f"r     = {r}\n"
                                f"NSE   = {NSE}\n"
                                f"KGE   = {KGE}\n"
                                f"Bias  = {bias}\n"
                                f"RMSE  = {RMSE}"
                            )
                        
                        ######################
                        #  FIGURE SECTION    #
                        ######################
                        fig, axes = plt.subplots(1, 2, figsize=(18, 9))
                        sf.plot(ax = axes[0], lw = 1.1, color= ['orangered','dodgerblue'])
                        axes[0].legend()
                        axes[0].set_ylabel('Streamflow (cfs)', fontsize = 12)
                        axes[0].set_title('Stream Flow', fontsize = 14)
                        axes[0].text( 0.98, 0.98, textstr, transform=axes[0].transAxes,
                            ha='right', va='top', fontsize=11,
                            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.8, boxstyle='round,pad=0.3'))
                        
                        
                        sim_sort = np.sort(sf.Sim)[::-1]
                        exc_sim = np.arange(1,len(sim_sort)+1)/len(sim_sort)*100
                        obs_sort = np.sort(sf.Obs)[::-1]
                        exc_obs = np.arange(1,len(obs_sort)+1)/len(obs_sort)*100
                        
                        axes[1].scatter(exc_sim,sim_sort, label = 'Sim',s = 10, color = 'orangered')
                        axes[1].scatter(exc_obs,obs_sort, label = 'Obs',s = 10, color = 'dodgerblue')
                        axes[1].set_xlabel('Exceedance Probability (%)', fontsize =12)
                        axes[1].set_ylabel('Flow (cfs)',  fontsize = 12)
                        axes[1].set_yscale('log')
                        axes[1].set_title('Flow Duration Curve', fontsize = 14)
                        axes[1].legend()
                        
                        plt.suptitle('USGS Gauge: '+site, fontsize = 14)
                        plt.subplots_adjust(hspace=0.0, wspace=0.14, 
                                            top=0.92,bottom=0.1,left=0.043,right=0.937)
                        plt.savefig(res_dir+'/PLOTS/Streaflow_'+site+'.png',dpi=400)
                        plt.close()



    elif dssat == True:
        res_dir = self.dbpath
        
        print('Plotting yield resutls...')
        
        
        #from geopy.geocoders import Nominatim
        #import requests
        import pandas as pd
        import numpy as np
        import os, sys,glob
        import geopandas as gpd
        #import time
        
        shp = self.shapefile
        startyear = self.startyear
        endyear = self.endyear
        tmp_dir = self.path
          
        gdf = gpd.read_file(shp)
        gdf = gdf.to_crs(epsg=4326)
        
        # get RHEAS yields
        dir_list = os.listdir(tmp_dir)
        if len(dir_list)<1:
            print('  ')
            print('*** NO DSSAT RUNS FOUND... YIELD PLOTS NOT GENERATED ****')
            print('  ')
        else:
            lats = []
            lons = []
            nsea = []
            latlons = []
            for d in dir_list:
                tmp = d.split('_')
                nsea.append(int(tmp[2]))
                latlons.append(tmp[0]+'_'+tmp[1])
            
            nsea = np.array(nsea)
            nsea = np.unique(nsea) 
            latlons = set(latlons)
            
            # get NASS yields for all the States and years DSSAT model was run for
            st = []
            fips = []
            county = []
            yrs = np.arange(startyear,endyear+1)
            
            for latlon in latlons:           
                clat = np.round(float(latlon.split('_')[1]),2)
                clon = np.round(float(latlon.split('_')[0]),2)
                lats.append(clat)
                lons.append(clon)
                fips_code, county_name, state_abbr = get_fips_from_latlon(clat, clon)
                st.append(state_abbr)
                fips.append(fips_code)
                county.append(county_name)
                
            loc_df = pd.DataFrame({'State':st,'County':county,
                                   'FIPS':fips,'Lat':lats, 'Lon':lons})
            all_nass_yields = []
            uniqe_states = list(loc_df.State.unique())
            nass_yield_db = pd.DataFrame()
            
            for yr in yrs:
                for state in uniqe_states:
                    try:
                        nass_df = fetch_nass_yields_by_state(state, yr)
                        nass_df['Year'] = yr
                        nass_df['State'] = state
                        all_nass_yields.append(nass_df)
                    except:
                        continue
                nass_yield_db = pd.concat(all_nass_yields)
            
            nass_yield_db.columns = ['FIPS','NASS_Yield','Year','State']
            
            rh_yield_db = pd.DataFrame()
            cyear = startyear 
            for season in nsea:
                rh_fips = []
                season_yield = []
                for latlon in latlons:
                    clat = np.round(float(latlon.split('_')[1]),2)
                    clon = np.round(float(latlon.split('_')[0]),2)
                    
                    fips_code, county_name, state_abbr = get_fips_from_latlon(clat, clon)
                    rh_fips.append(fips_code)
                    
                    flist = glob.glob(tmp_dir+'/'+latlon+'_'+str(season)+'/PLANTGRO*.OUT')
                    ens_yield = []
                    if len(flist)>0:
                        for f in flist:
                            df = pd.read_table(f, sep = '\s+').GWAD
                            g = df.max()
                            if g < 6000:
                                g = g*1.25
                                
                            ens_yield.append(g)
                        season_yield.append(np.nanmedian(ens_yield))
                    else:
                        season_yield = np.nan
           
                tmp_df = pd.DataFrame({'Lat':lats,'Lon':lons,
                                           'Sim_Yield':season_yield,
                                           'FIPS': rh_fips,
                                           'Year': cyear})
                if season ==0:
                    rh_yield_db = tmp_df
                else:
                    rh_yield_db = pd.concat([rh_yield_db,tmp_df])
                cyear = cyear+1
            
            merged = pd.merge(rh_yield_db, nass_yield_db, on=['FIPS', 'Year'], how='left')
            yield_data = merged.copy()
            
            # fill in nan values with other county or state mean values 
            st_fips = []
            for ips in fips:
                st_fips.append(ips[:2])
            st_fips = set(st_fips)
            
            for yr in yrs:
                for stf in st_fips:
                    #print(yr, stf)
                    t = nass_yield_db.loc[
                        (nass_yield_db['FIPS'] == stf+'998') & (nass_yield_db['Year'] == yr)]
                    other_yield = t.NASS_Yield.median()
                    for idx,row in yield_data.iterrows():
                        if row.Year == yr and np.isnan(row.NASS_Yield):
                            yield_data.loc[idx,'NASS_Yield'] = other_yield
                                   
            yield_gdf = gpd.GeoDataFrame(yield_data, 
                                         geometry=gpd.points_from_xy(yield_data.Lon, yield_data.Lat),
                                         crs = 'EPSG:4326')
            joined = gpd.sjoin(yield_gdf,gdf,  how='right', predicate='within')  
                        
                           
            ######################
            #  FIGURE SECTION    #
            ######################
            import matplotlib as mpl
            for yr in yrs:
                ndf = joined.loc[joined.Year == yr]
                if len(ndf)>0:
                    #-----------
                    fig, axes = plt.subplots(1, 2, figsize=(16, 8), sharex=True, sharey=False)
                    
                    # Normalize color scale across both panels
                    vmin = 3000
                    vmax = 11000
                    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
                    cmap = plt.cm.YlGn
                    
                    # Plot background HUC shapefile in both panels
                    #for ax in axes:
                    #    gdf.plot(ax=ax, facecolor='none', edgecolor='gray', linewidth=0.75)
                    
                    ndf.plot(
                        ax=axes[0],
                        column='Sim_Yield',
                        cmap=cmap,
                        markersize=50,
                        legend=False,
                        norm=norm)
                    
                    axes[0].set_title("Simulated Yields (kg/ha)", fontsize=14)
                    axes[0].set_xlabel("Longitude")
                    axes[0].set_ylabel("Latitude")
                    
                    ndf.plot(
                        ax=axes[1],
                        column='NASS_Yield',
                        cmap=cmap,
                        markersize=50,
                        legend=False,
                        norm=norm)
                    
                    axes[1].set_title("NASS Yield (kg/ha)", fontsize = 14)
                    axes[1].set_xlabel("Longitude")
                    axes[1].set_ylabel("Latitude")
                    
                    for ax in axes:
                        gdf.plot(ax=ax, facecolor='none', edgecolor='gray', linewidth=0.75)
                        
                    plt.suptitle(str(yr), fontsize=16)
                    
                    # Create shared colorbar
                    cbar_ax = fig.add_axes([0.92, 0.25, 0.02, 0.5])  # adjust as needed
                    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
                    sm._A = []  # dummy mappable
                    cbar = fig.colorbar(sm, cax=cbar_ax,extend='max')
                    cbar.set_label("Yield (Kg/ha)")
                    
                    plt.tight_layout()  # leave space for colorbar
                    plt.savefig(res_dir+'/PLOTS/Yield_Map_'+str(yr)+'.png',dpi=400)
                    plt.close()


            
        