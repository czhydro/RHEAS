""" RHEAS module for generating drought products.

.. module:: drought
   :synopsis: Module that contains functionality for generating drought products

.. moduleauthor:: Kostas Andreadis <kandread@umass.edu>

"""

import logging
import threading
from datetime import date, datetime

import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rx
import scipy.stats as stats
from dateutil.relativedelta import relativedelta
from numba import float64, guvectorize

from . import dbio


def calcVCI(model, table="ndvi.modis"):
    """Calculate Vegetation Condition Index."""
    log = logging.getLogger(__name__)
    sdate = date(model.startyear, model.startmonth, model.startday).strftime("%Y-%m-%d")
    edate = date(model.endyear, model.endmonth, model.endday).strftime("%Y-%m-%d")
    if dbio.tableExists(model.dbname, table.split(".")[0], table.split(".")[1]):
        db = dbio.connect(model.dbname)
        cur = db.cursor()
        cur.execute("drop table if exists ndvi_max, ndvi_min, ndvi_max_min, f1")
        db.commit()
        sql = "create table ndvi_max as (select st_union(rast, 'MAX') as rast from {0})".format(table)
        cur.execute(sql)
        sql = "create table ndvi_min as (select st_union(rast, 'MIN') as rast from {0})".format(table)
        cur.execute(sql)
        sql = "create table ndvi_max_min as (select st_mapalgebra(max.rast, 1, min.rast, 1, '[rast1]-[rast2]') as rast from ndvi_max as max, ndvi_min as min)"
        cur.execute(sql)
        db.commit()
        sql = "create table f1 as (select fdate, st_mapalgebra(f.rast, 1, min.rast, 1, '[rast1]-[rast2]') as rast from {0} as f, ndvi_min as min where fdate>=date'{1}' and fdate<=date'{2}' group by fdate,f.rast,min.rast)".format(table, sdate, edate)
        cur.execute(sql)
        db.commit()
        if dbio.tableExists(model.dbname, model.name, "vci"):
            sql = "insert into {0}.vci (fdate, rast) select fdate, st_mapalgebra(f1.rast, 1, mm.rast, 1, '[rast1]/([rast2]+0.0001)') as rast from f1, ndvi_max_min as mm group by fdate,f1.rast,mm.rast".format(model.name)
        else:
            sql = "create table {0}.vci as (select fdate, st_mapalgebra(f1.rast, 1, mm.rast, 1, '[rast1]/([rast2]+0.0001)') as rast from f1, ndvi_max_min as mm group by fdate,f1.rast,mm.rast)".format(model.name)
        cur.execute(sql)
        db.commit()
        cur.execute("drop table ndvi_max")
        cur.execute("drop table ndvi_min")
        cur.execute("drop table ndvi_max_min")
        cur.execute("drop table f1")
        db.commit()
        cur.close()
        db.close()
    else:
        log.warning("No NDVI data were found in database. Cannot calculate VCI!")
    return None


@guvectorize("(float64[:], float64[:])", "(n) -> (n)")
def spri(s, out):
    bound = 3.09
    c0 = 2.515517
    c1 = 0.802583
    c2 = 0.010328
    d1 = 1.432788
    d2 = 0.189269
    d3 = 0.001308
    sf = np.where(s > 0, s, np.nan)
    sf = sf[~np.isnan(sf)]
    p = np.zeros(len(sf))
    p[np.argsort(sf)] = (np.arange(len(sf)) + 1 - 0.44) / (len(sf) + 0.12)
    t = np.where(p <= 0.5, np.sqrt(np.log(1 / p**2)), np.sqrt(np.log(1 / (1 - p)**2)))
    val = t - (c0 + c1 * t + c2 * t**2) / (1 + d1 + d2 * t**2 + d3 * t**3)
    si = np.zeros(len(sf))
    si[p <= 0.5] = -val[p <= 0.5]
    si[p > 0.5] = val[p > 0.5]
    out[:] = np.zeros(len(s))
    out[((~np.isnan(s)) & (s > 0))] = si
    out[np.isnan(s)] = np.nan
    out[s == 0] = -bound
    out[:] = np.clip(out, -bound, bound)


def xr_spri(data, dim):
    fspri = xr.apply_ufunc(
        spri,
        data,
        input_core_dims=[[dim]],
        output_core_dims=[[dim]],
        dask="parallelized",
        output_dtypes=[data.dtype],
    )
    return fspri


@guvectorize("(float64[:], float64[:])", "(n) -> (n)")
def severity(s, out):
    r = np.argsort(s)
    out[:] = np.empty_like(s)
    out[r] = (np.arange(len(s)) + 1) / len(s) * 100


def xr_sev(data, dim):
    fsev = xr.apply_ufunc(
        severity,
        data,
        input_core_dims=[[dim]],
        output_core_dims=[[dim]],
        dask="parallelized",
        output_dtypes=[data.dtype],
    )
    return fsev


def calcSeverity(model, ensemble, varname="soil_moist"):
    """Calculate drought severity from *climatology* table stored in database."""
    startdate = datetime(model.startyear + model.skipyear, model.startmonth, model.startday)
    enddate = datetime(model.endyear, model.endmonth, model.endday)
    if bool(ensemble):
        equery = "and ((ensemble={0} and fdate>=date'{1}') or (ensemble=0 and fdate<date'{1}'))".format(ensemble, startdate.strftime("%Y-%m-%d"))
    else:
        equery = ""
    db = dbio.connect(model.dbname)
    cur = db.cursor()
    if varname == "soil_moist":
        files = {}
        for lyr in range(2):
            sql = "select fdate,filename from {0}.soil_moist where layer={3} and fdate<='{1}' {2} order by fdate".format(model.name, enddate.strftime("%Y-%m-%d"), equery, lyr+1)
            cur.execute(sql)
            results = cur.fetchall()
            for r in results:
                name = "{0}/{1}".format(model.dbpath, r[1])
                files[r[0]] = files[r[0]] + [name] if r[0] in files else [name]
        x = []
        for t in files:
            x1 = rx.open_rasterio(files[t][0], chunks=True, lock=threading.Lock())
            x2 = rx.open_rasterio(files[t][1], chunks=True, lock=threading.Lock())
            x.append(x1 + x2)
        dates = list(files.keys())
    else:
        sql = "select fdate,filename from {0}.runoff where fdate<='{1}' {2} group by fdate order by fdate".format(model.name, enddate.strftime("%Y-%m-%d"), equery)
        cur.execute(sql)
        results = cur.fetchall()
        files = ["{0}/{1}".format(model.dbpath, r[1]) for r in results]
        dates = [r[0] for r in results]
        x = []
        for f in files:
            x.append(rx.open_rasterio(f, chunks=True, lock=threading.Lock()))
    x = xr.concat(x, dim='time')
    x = x.assign_coords({'time': dates})
    x = x.squeeze().reset_coords('band', drop=True)
    x = x.chunk(dict(time=len(files), x='auto', y='auto'))
    y = x.rolling(time=10).mean() # decad mean
    s = xr_sev(y, 'time')
    sev = s.sel(time=slice(startdate.date(), enddate.date())).transpose('time', 'y', 'x').data
    return sev


def calcSI(duration, model, ensemble, varname):
    """Calculate Standardized Drought Index for specified month
    *duration*."""
    log = logging.getLogger(__name__)
    startdate = datetime(model.startyear + model.skipyear, model.startmonth, model.startday)
    enddate = datetime(model.endyear, model.endmonth, model.endday)
    if bool(ensemble):
        equery = "and ((ensemble={0} and fdate>=date'{1}') or (ensemble=0 and fdate<date'{1}'))".format(ensemble, startdate.strftime("%Y-%m-%d"))
    else:
        equery = ""
    dstartdate = startdate - relativedelta(months=duration)
    if dstartdate > startdate:
        dstartdate = startdate
    db = dbio.connect(model.dbname)
    cur = db.cursor()
    sql = "select count(fdate) from {0}.{4} where fdate>=date'{1}' and fdate<=date'{2}' {3}".format(model.name, dstartdate.strftime("%Y-%m-%d"), enddate.strftime("%Y-%m-%d"), equery, varname)
    cur.execute(sql)
    nt = cur.fetchone()[0]
    ndays = (enddate - dstartdate).days + 1
    if duration < 1 or (ndays > nt and nt < duration * 30):
        log.warning("Cannot calculate {1} with {0} months duration.".format(duration, "SPI" if varname == "rainf" else "SRI"))
        spri = None
    else:
        sql = "select fdate,filename from {0}.{3} where fdate<=date'{1}' {2} group by fdate,filename order by fdate".format(model.name, enddate.strftime("%Y-%m-%d"), equery, varname)
        cur.execute(sql)
        results = cur.fetchall()
        files = ["{0}/{1}".format(model.dbpath, r[1]) for r in results]
        dates = [r[0] for r in results]
        x = []
        for f in files:
            x.append(rx.open_rasterio(f, chunks=True, lock=threading.Lock()))
        x = xr.concat(x, dim='time')
        x = x.assign_coords({'time': dates})
        x = x.squeeze().reset_coords('band', drop=True)
        x = x.chunk(dict(time=len(files), x='auto', y='auto'))
        y = x.rolling(time=duration*30).sum() # assume each month is 30 days
        s = xr_spri(y, 'time')
        spri = s.sel(time=slice(startdate.date(), enddate.date())).transpose('time', 'y', 'x').data
    return spri


def calcSPI(duration, model, ensemble):
    """Calculate Standardized Precipitation Index for specified month
    *duration*."""
    spi = calcSI(duration, model, ensemble, "rainf")
    return spi


def calcSRI(duration, model, ensemble):
    """Calculate Standardized Runoff Index for specified month
    *duration*."""
    sri = calcSI(duration, model, ensemble, "runoff")
    return sri


@guvectorize("(float64[:], float64[:])", "(n) -> (n)")
def dryspells(s, out):
    # criterion for defining a dry day
    drought_thresh = 0.0
    # number of consecutive days above drought threshold to break dry spell
    recovduration = 2
    # amount of daily rainfall to break a dry spell
    break_thresh = 10.0
    break_thresh_s = 5.0
    rdt = recovduration - 1
    out[:] = np.zeros(len(s))
    drydays = (s <= drought_thresh).astype('int')
    for t in range(1, len(s)):
        if t >= rdt and ((drydays[t-rdt:t+1].sum() == 0 and s[t-rdt:t+1].sum() >= break_thresh_s) or s[t] >= break_thresh):
            out[t] = 0
        else:
            out[t] = (out[t-1] + 1)


def xr_dsp(data, dim):
    fdsp = xr.apply_ufunc(
        dryspells,
        data,
        input_core_dims=[[dim]],
        output_core_dims=[[dim]],
        dask="parallelized",
        output_dtypes=[data.dtype],
    )
    return fdsp


def calcDrySpells(model, ensemble):
    """Calculate maps of number of dry spells during simulation period."""
    # FIXME: Currently only uses precipitation to identify dry spells. Need to change it to also use soil moisture and runoff
    startdate = datetime(model.startyear + model.skipyear, model.startmonth, model.startday)
    enddate = datetime(model.endyear, model.endmonth, model.endday)
    if bool(ensemble):
        equery = "and ((ensemble={0} and fdate>=date'{1}') or (ensemble=0 and fdate<date'{1}'))".format(ensemble, startdate.strftime("%Y-%m-%d"))
    else:
        equery = ""
    db = dbio.connect(model.dbname)
    cur = db.cursor()
    sql = "select fdate,filename from {0}.rainf where fdate>=date'{1}-{2}-{3}' and fdate<=date'{4}-{5}-{6}' {7} order by fdate".format(model.name, model.startyear, model.startmonth, model.startday, model.endyear, model.endmonth, model.endday, equery)
    cur.execute(sql)
    results = cur.fetchall()
    files = ["{0}/{1}".format(model.dbpath, r[1]) for r in results]
    dates = [r[0] for r in results]
    x = []
    for f in files:
        x.append(rx.open_rasterio(f, chunks=True, lock=threading.Lock()))
    x = xr.concat(x, dim='time')
    x = x.assign_coords({'time': dates})
    x = x.squeeze().reset_coords('band', drop=True)
    x = x.chunk(dict(time=len(files), x='auto', y='auto'))
    s = xr_dsp(x, 'time')
    dsp = s.transpose('time', 'y', 'x').data
    return dsp


@guvectorize("(float64[:], int32[:], float64[:])", "(n),(n) -> (n)")
def smdi(s, w, out):
    out[:] = np.empty_like(s)
    msw = np.zeros(len(s))
    maxsw = np.zeros(len(s))
    minsw = np.zeros(len(s))
    for i in range(1, max(w)+1):
        sw = s[w == i]
        msw[w == i] = np.median(sw)
        maxsw[w == i] = np.max(sw)
        minsw[w == i] = np.min(sw)
    sd = (s - msw) / ((s <= msw) * (msw - minsw) + (s > msw) * (maxsw - msw) + 1e-6) * 100
    out[0] = sd[0] / 50
    for i in range(1, len(s)):
        out[i] = 0.5 * out[i-1] + sd[i] / 50


def xr_smdi(data, dim):
    w = np.array(pd.to_datetime(data.time).isocalendar().week, 'int32')
    xw = xr.zeros_like(data, 'int32')
    xw[:, :, :] = w[:, None, None]
    fsmdi = xr.apply_ufunc(
        smdi,
        data,
        xw,
        input_core_dims=[[dim], [dim]],
        output_core_dims=[[dim]],
        dask="parallelized",
        output_dtypes=[data.dtype],
    )
    return fsmdi


def calcSMDI(model, ensemble):
    """Calculate Soil Moisture Deficit Index (Narasimhan & Srinivasan, 2005)."""
    startdate = datetime(model.startyear + model.skipyear, model.startmonth, model.startday)
    enddate = datetime(model.endyear, model.endmonth, model.endday)
    if bool(ensemble):
        equery = "and ((ensemble={0} and fdate>=date'{1}') or (ensemble=0 and fdate<date'{1}'))".format(ensemble, startdate.strftime("%Y-%m-%d"))
    else:
        equery = ""
    db = dbio.connect(model.dbname)
    cur = db.cursor()
    files = {}
    for lyr in range(2):
        sql = "select fdate,filename from {0}.soil_moist where layer={3} and fdate<='{1}' {2} order by fdate".format(model.name, enddate.strftime("%Y-%m-%d"), equery, lyr+1)
        cur.execute(sql)
        results = cur.fetchall()
        for r in results:
            name = "{0}/{1}".format(model.dbpath, r[1])
            files[r[0]] = files[r[0]] + [name] if r[0] in files else [name]
    x = []
    for t in files:
        x1 = rx.open_rasterio(files[t][0], chunks=True, lock=threading.Lock())
        x2 = rx.open_rasterio(files[t][1], chunks=True, lock=threading.Lock())
        x.append(x1 + x2)
    dates = list(files.keys())
    x = xr.concat(x, dim='time')
    x = x.assign_coords({'time': dates})
    x = x.squeeze().reset_coords('band', drop=True)
    x = x.chunk(dict(time=len(files), x='auto', y='auto'))
    y = x.rolling(time=7).mean().fillna(x)
    s = xr_smdi(y, 'time')
    sd = s.sel(time=slice(startdate.date(), enddate.date())).transpose('time', 'y', 'x').data
    return sd


@guvectorize("(float64[:], float64, float64, float64, float64, float64[:])", "(n),(),(),(),() -> (n)")
def suctionhead(data, z, n, b, psi, out):
    pf = np.log(psi * ((data / z) / n)**(-b))
    out[:] = (pf - pf.mean()) / pf.std()


def xr_suctionhead(data, za, na, ba, psia, dim):
    fsh = xr.apply_ufunc(
        suctionhead,
        data,
        za,
        na,
        ba,
        psia,
        input_core_dims=[[dim], [], [], [], []],
        output_core_dims=[[dim]],
        dask="parallelized",
        output_dtypes=[data.dtype],
    )
    return fsh


def calcSuctionHead(model, ensemble, nlayers=3):
    """Calculate soil suction from soil moisture using the Clapp
    and Hornberger (1978) model and parameters."""
    startdate = datetime(model.startyear + model.skipyear, model.startmonth, model.startday)
    enddate = datetime(model.endyear, model.endmonth, model.endday)
    if bool(ensemble):
        equery = "and ((ensemble={0} and fdate>=date'{1}') or (ensemble=0 and fdate<date'{1}'))".format(ensemble, startdate.strftime("%Y-%m-%d"))
    else:
        equery = ""
    Ksat = np.array([63.36, 56.16, 12.49, 2.59, 2.5, 2.27, 0.612, 0.882, 0.781, 0.371, 0.461])
    Ksat *= (10 * 24.)  # convert from cm/hr to mm/day
    n = np.array([.395, .41, .435, .485, .451, .42, .477, .476, .426, .492, .482])
    psi_a = np.array([121., 90., 218., 786., 478., 299., 356., 630., 153., 490., 405.])
    b = np.array([4.05, 4.38, 4.9, 5.3, 5.39, 7.12, 7.75, 8.52, 10.4, 10.4, 11.4])
    db = dbio.connect(model.dbname)
    cur = db.cursor()
    files = {}
    for lyr in range(2):
        sql = "select fdate,filename from {0}.soil_moist where layer={3} and fdate<='{1}' {2} order by fdate".format(model.name, enddate.strftime("%Y-%m-%d"), equery, lyr+1)
        cur.execute(sql)
        results = cur.fetchall()
        for r in results:
            name = "{0}/{1}".format(model.dbpath, r[1])
            files[r[0]] = files[r[0]] + [name] if r[0] in files else [name]
    x = []
    for t in files:
        x1 = rx.open_rasterio(files[t][0], chunks=True, lock=threading.Lock())
        x2 = rx.open_rasterio(files[t][1], chunks=True, lock=threading.Lock())
        x.append(x1 + x2)
    dates = list(files.keys())
    x = xr.concat(x, dim='time')
    x = x.assign_coords({'time': dates})
    x = x.squeeze().reset_coords('band', drop=True)
    x = x.chunk(dict(time=len(files), x='auto', y='auto'))
    za = xr.zeros_like(x.isel(time=0))
    psia = xr.zeros_like(x.isel(time=0))
    na = xr.zeros_like(x.isel(time=0))
    ba = xr.zeros_like(x.isel(time=0))
    soil = pd.read_csv("{0}/soil.txt".format(model.model_path), header=None, delim_whitespace=True, usecols=[2, 3]+list(range(4*nlayers+10, 5*nlayers+10))+list(range(9+nlayers,nlayers+12))).values
    xul = min(x.x.data)
    yul = max(x.y.data)
    lati = ((yul - soil[:, 0]) / model.res).astype('int')
    lonj = ((soil[:, 1] - xul) / model.res).astype('int')
    c = np.ravel_multi_index((lati, lonj), za.shape)
    za_ = np.zeros(za.shape[0] * za.shape[1])
    za_[c] = soil[:, 2:2+nlayers].sum(axis=1)
    za.data[:] = za_.reshape(za.shape)
    del za_
    k = soil[:, 5:].mean(axis=1)
    ki = np.argmin(np.abs(k[:, None] - Ksat), axis=1)
    na_ = np.zeros(na.shape[0] * na.shape[1])
    na_[c] = n[ki]
    na.data[:] = na_.reshape(na.shape)
    del na_
    psia_ = np.zeros(psia.shape[0] * psia.shape[1])
    psia_[c] = psi_a[ki]
    psia.data[:] = psia_.reshape(psia.shape)
    del psia_
    ba_ = np.zeros(ba.shape[0] * ba.shape[1])
    ba_[c] = b[ki]
    ba.data[:] = ba_.reshape(ba.shape)
    del ba_
    y = x.rolling(time=10).mean().fillna(x)
    pf = xr_suctionhead(y, za, na, ba, psia, 'time')
    return pf.sel(time=slice(startdate.date(), enddate.date())).transpose('time', 'y', 'x').data


def calcFpar(model, ensemble):
    """Retrieve the Photosynthetically Active Radiation from the model simulation."""
    if bool(ensemble):
        equery = "where (ensemble={0} or ensemble=0)".format(ensemble)
    else:
        equery = ""
    startdate = datetime(model.startyear + model.skipyear, model.startmonth, model.startday)
    enddate = datetime(model.endyear, model.endmonth, model.endday)
    db = dbio.connect(model.dbname)
    cur = db.cursor()
    sql = "select fdate,filename from {0}.par {1} order by fdate".format(model.name, equery)
    cur.execute(sql)
    results = cur.fetchall()
    files = ["{0}/{1}".format(model.dbpath, r[1]) for r in results]
    dates = [r[0] for r in results]
    x = []
    for f in files:
        x.append(rx.open_rasterio(f, chunks=True, lock=threading.Lock()))
    x = xr.concat(x, dim='time')
    x = x.assign_coords({'time': dates})
    x = x.squeeze().reset_coords('band', drop=True)
    x = x.chunk(dict(time=len(files), x='auto', y='auto'))
    y = x.rolling(time=10).mean().fillna(x)
    fpar = (y - y.mean()) / y.std()
    return fpar.sel(time=slice(startdate.date(), enddate.date())).transpose('time', 'y', 'x').data


def calcCDI(model, ensemble):
    """Calculate Combined Drought Index as a monthly time series. The index is
    categorical with the values corresponding to:
    0 = No drought
    1 = Watch (Precipitation deficit)
    2 = Warning (Soil moisture deficit)
    3 = Alert 1 (Vegetation stress following precipitation deficit)
    4 = Alert 2 (Vegetation stress following precipitation/soil moisture deficit)."""
    log = logging.getLogger(__name__)
    spi = calcSPI(3, model, ensemble)
    sma = calcSuctionHead(model, ensemble)
    fapar = calcFpar(model, ensemble)
    if all(v is not None for v in [spi, sma, fapar]):
        cdi = np.zeros(spi.shape, dtype='int')
        cdi[spi < -1] = 1
        cdi[(fapar > 1) & (spi < -1)] = 2
        cdi[(fapar < -1) & (spi < -1)] = 3
        cdi[(fapar < -1) & (sma > 1) & (spi < -1)] = 4
    else:
        log.warning("Error in calculating SPI-3, SMA or PAR. Cannot calculate CDI!")
        cdi = None
    return cdi


def calc(varname, model, ensemble):
    """Calculate drought-related variable."""
    if varname.startswith("spi"):
        duration = int(varname[3:])
        output = calcSPI(duration, model, ensemble)
    elif varname.startswith("sri"):
        duration = int(varname[3:])
        output = calcSRI(duration, model, ensemble)
    elif varname == "severity":
        output = calcSeverity(model, ensemble)
    elif varname == "cdi":
        output = calcCDI(model, ensemble)
    elif varname == "smdi":
        output = calcSMDI(model, ensemble)
    elif varname == "dryspells":
        output = calcDrySpells(model, ensemble)
    elif varname == "vci":
        output = calcVCI(model)
    return output
