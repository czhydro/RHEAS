""" Class definition for the VIC model interface

.. module:: vic
   :synopsis: Definition of the VIC model class

.. moduleauthor:: Kostas Andreadis <kandread@umass.edu>

"""

import decimal
import logging
import multiprocessing as mp
import os
import random
import shutil
import string
import subprocess
import sys
import tempfile
from collections import OrderedDict
from datetime import date, datetime, timedelta

import numpy as np
import pandas as pd

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

from .. import dbio, drought
from ..raster import TileReader
from . import output as vicoutput


class VIC:

    def __init__(self, path, dbname, resolution, startyear, startmonth, startday,
                 endyear, endmonth, endday, name="", savestate="", nlayer=3):
        log = logging.getLogger(__name__)
        self.model_path = path
        self.data_path = os.path.dirname(os.path.abspath(__file__)) + "/../../data"
        if not os.path.isdir(self.data_path):
            try:
                self.data_path = os.environ['RHEASDATA']
            except KeyError:
                log.error("Data directory is not accessible. Please set RHEASDATA environment variable. Exiting...")
                sys.exit()
        # FIXME: pass this as an option from the configuration file
        self.nodata = -9999.
        if bool(name):
            self.name = name
        else:
            #self.name = None    # OLD CODE
            self.name = "unnamed"  # use a default name for the simulation
        self.dbpath = "./output/{0}".format(self.name)
        if not os.path.exists(self.dbpath):
            os.makedirs(self.dbpath)
        self.startyear = startyear
        self.startmonth = startmonth
        self.startday = startday
        self.startdate = datetime(startyear, startmonth, startday)
        self.endyear = endyear
        self.endmonth = endmonth
        self.endday = endday
        self.enddate = datetime(endyear, endmonth, endday)
        self.nlayers = nlayer
        self.dbname = dbname
        db = dbio.connect(dbname)
        cur = db.cursor()
        cur.execute(
            "select schema_name from information_schema.schemata where schema_name='{0}'".format(self.name))
        if not bool(cur.rowcount):
            cur.execute("create schema {0}".format(self.name))
            db.commit()
        cur.execute(
            "select resolution from vic.input order by abs(resolution - {0})".format(resolution))
        if not bool(cur.rowcount):
            log.error("No appropriate VIC input files found in the database. Exiting!")
            sys.exit()
        self.res = cur.fetchone()[0]
        cur.close()
        self.grid_decimal = -(decimal.Decimal(str(self.res)).as_tuple().exponent - 1)
        self.lat = []
        self.lon = []
        self.gid = OrderedDict()
        self.lgid = OrderedDict()
        self.depths = OrderedDict()
        self.skipyear = 0
        self.elev = OrderedDict()
        self.statefile = ""

    def paramFromDB(self):
        """Retrieve file parameters from database."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        # cur = self.db.cursor()
        cur.execute(
            'select veglib,vegparam,snowbandfile from vic.input where resolution=%f;' % self.res)
        veglib, vegparam, snowbands = cur.fetchone()
        cur.close()
        db.close()
        return veglib, vegparam, snowbands

    def _getSnowbands(self, snowbands):
        """Find number of snow bands from file."""
        filename = "{0}/{1}".format(self.data_path, snowbands)
        with open(filename) as f:
            line = f.readline()
        return int((len(line.split()) - 1) / 3)

    def writeSoilFile(self, shapefile, toast=False):
        """Write soil parameter file for current simulation based on basin shapefile.
        Optionally handle TOAST database tables (see https://postgis.net/docs/performance_tips.html#small_tables_large_objects).
        """
        log = logging.getLogger(__name__)
        temptable = ''.join(random.SystemRandom().choice(string.ascii_letters) for _ in range(8))
        with tempfile.NamedTemporaryFile('w') as fout:
            subprocess.call(["shp2pgsql", "-s", "4326", "-d", "-I", shapefile, temptable], stdout=fout)
            cmd = ["psql", "-d", self.dbname, "-f", fout.name]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            for line in iter(proc.stdout.readline, b''):
                log.debug(line.strip())
            proc.wait()
        # overwrite basin table if it exists
        dbio.deleteTable(self.dbname, self.name, "basin")
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        db.commit()
        if toast:
            cur.execute("alter table {0} add column bbox geometry".format(temptable))
            cur.execute("update {0} set bbox=st_envelope(st_force2d(geom))".format(temptable))
            db.commit()
            sql = "create table {0}.basin as (select distinct v.id gid, v.elev, v.depths, v.geom, v.line from vic.soils v,{1} t where t.bbox && v.geom and v.resolution={2})".format(self.name, temptable, self.res)
        else:
            sql = "create table {0}.basin as (select distinct v.id gid, v.elev, v.depths, v.geom, v.line from vic.soils v, {1} t where st_intersects(v.geom, t.geom) and v.resolution={2})".format(self.name, temptable, self.res)
        cur.execute(sql)
        db.commit()
        cur.execute("drop table {0}".format(temptable))
        cur.execute("create index basin_i on {0}.basin(gid)".format(self.name))
        cur.execute("create index basin_s on {0}.basin using gist(geom)".format(self.name))
        db.commit()
        sql = "select line,gid,st_y(geom),st_x(geom),elev,depths from {0}.basin order by gid".format(
            self.name)
        cur.execute(sql)
        lines = cur.fetchall()
        with open(self.model_path + '/soil.txt', 'w') as fout:
            for line in lines:
                gid, lat, lon, elev, depths = line[1:]
                fout.write("{0}\n".format(line[0]))
                self.lat.append(lat)
                self.lon.append(lon)
                self.gid[gid] = (lat, lon)
                self.lgid[(lat, lon)] = gid
                self.depths[gid] = depths
                self.elev[gid] = elev
        cur.execute("alter table {0}.basin drop column line".format(self.name))
        cur.close()
        db.close()

    def stateFile(self):
        """Retrieve state file path from database."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        sql = "select filename from {0}.state where fdate = date '{1}-{2}-{3}'".format(
            self.name, self.startyear, self.startmonth, self.startday)
        try:
            cur.execute(sql)
            result = cur.fetchone()
        except:
            result = False
        if bool(result):
            filename = result[0]
        else:
            filename = None
        cur.close()
        db.close()
        return filename

    def _stateToDb(self, statefilepath):
        """Add path to state file into database."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        cur.execute(
            "select schema_name from information_schema.schemata where schema_name='{0}'".format(self.name))
        if not bool(cur.rowcount):
            cur.execute("create schema {0}".format(self.name))
            db.commit()
        cur.execute(
            "select table_name from information_schema.tables where table_schema='{0}' and table_name='state'".format(self.name))
        if not bool(cur.rowcount):
            sql = "create table {0}.state (filename text, fdate date)".format(
                self.name)
            cur.execute(sql)
            db.commit()
        statefile = "{0}/vic.state_{1:04d}{2:02d}{3:02d}".format(
            statefilepath, self.endyear, self.endmonth, self.endday)
        statedate = "{0}-{1}-{2}".format(self.endyear,
                                         self.endmonth, self.endday)
        cur.execute(
            "select * from {0}.state where fdate=date '{1}'".format(self.name, statedate))
        if bool(cur.rowcount):
            sql = "update {0}.state set filename='{1}' where fdate=date '{2}'".format(
                self.name, statefile, statedate)
        else:
            sql = "insert into {0}.state values ('{1}', date '{2}')".format(
                self.name, statefile, statedate)
        cur.execute(sql)
        db.commit()
        cur.close()
        db.close()

    def writeParamFile(self, nodes=3, time_step=24, save_state="", init_state=False, state_file="", save_state_to_db=False):
        """Write VIC global parameter file for current simulation."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        cur.execute(
            'select rootzones from vic.input where resolution=%f;' % self.res)
        root_zones = cur.fetchone()[0]
        cur.close()
        db.close()
        # self.nlayers = nlayer
        fout = open(self.model_path + '/global.txt', 'w')
        fout.write("NLAYER\t{0:d}\n".format(self.nlayers))
        fout.write("NODES\t{0:d}\n".format(nodes))
        if time_step < 24:
            fout.write("TIME_STEP\t{0:d}\nSNOW_STEP\t{1:d}\n".format(
                time_step, time_step))
            fout.write("FULL_ENERGY\tTRUE\nFROZEN_SOIL\tTRUE\n")
        else:
            fout.write("TIME_STEP\t24\nSNOW_STEP\t3\n")
            fout.write("FULL_ENERGY\tFALSE\nFROZEN_SOIL\tFALSE\n")
        fout.write("STARTYEAR\t{0:04d}\n".format(self.startyear))
        fout.write("STARTMONTH\t{0:02d}\n".format(self.startmonth))
        fout.write("STARTDAY\t{0:02d}\n".format(self.startday))
        fout.write("ENDYEAR\t{0:04d}\n".format(self.endyear))
        fout.write("ENDMONTH\t{0:02d}\n".format(self.endmonth))
        fout.write("ENDDAY\t{0:02d}\n".format(self.endday))
        fout.write("IMPLICIT\tFALSE\nTFALLBACK\tTRUE\n")
        fout.write("SNOW_ALBEDO\tUSACE\nSNOW_DENSITY\tDENS_SNTHRM\n")
        fout.write("BLOWING\tFALSE\nCOMPUTE_TREELINE\tFALSE\n")
        fout.write("DIST_PRCP\tFALSE\nPREC_EXPT\t0.6\nCORRPREC\tFALSE\n")
        fout.write("MAX_SNOW_TEMP\t0.5\nMIN_RAIN_TEMP\t-0.5\n")
        fout.write("MIN_WIND_SPEED\t0.1\nAERO_RESIST_CANSNOW\tAR_406_FULL\n")
        if bool(state_file):
            statefile = state_file
        elif init_state:
            statefile = self.stateFile()
        else:
            statefile = None
        if statefile:
            fout.write("INIT_STATE\t{0:s}\n".format(statefile))
        if bool(save_state):
            if isinstance(save_state, str):
                if not os.path.isdir(save_state):
                    os.mkdir(save_state)
                fout.write("STATENAME\t{0}/vic.state\n".format(save_state))
            else:
                fout.write(
                    "STATENAME\t{0}/vic.state\n".format(self.model_path))
            fout.write("STATEYEAR\t{0:04d}\n".format(self.endyear))
            fout.write("STATEMONTH\t{0:02d}\n".format(self.endmonth))
            fout.write("STATEDAY\t{0:02d}\n".format(self.endday))
            self.statefile = "vic.state_{0:04d}{1:02d}{2:02d}".format(
                self.endyear, self.endmonth, self.endday)
            if save_state_to_db:
                self._stateToDb(save_state)
        fout.write("BINARY_STATE_FILE\tFALSE\n")
        fout.write(
            "FORCING1\t{0:s}/data_\n".format(self.model_path + "/forcings"))
        fout.write("FORCE_FORMAT\tASCII\nFORCE_ENDIAN\tLITTLE\nN_TYPES\t4\n")
        fout.write("FORCE_TYPE\tPREC\n")
        fout.write("FORCE_TYPE\tTMAX\n")
        fout.write("FORCE_TYPE\tTMIN\n")
        fout.write("FORCE_TYPE\tWIND\n")
        fout.write("FORCE_DT\t24\n")
        fout.write("FORCEYEAR\t{0:04d}\n".format(self.startyear))
        fout.write("FORCEMONTH\t{0:02d}\n".format(self.startmonth))
        fout.write("FORCEDAY\t{0:02d}\n".format(self.startday))
        fout.write("FORCEHOUR\t0\n")
        fout.write("GRID_DECIMAL\t{0:d}\n".format(self.grid_decimal))
        fout.write("WIND_H\t10.0\nMEASURE_H\t2.0\nALMA_INPUT\tFALSE\n")
        fout.write("SOIL\t{0:s}\n".format(self.model_path + "/soil.txt"))
        veglib, vegparam, snowbands = self.paramFromDB()
        fout.write("VEGLIB\t{0}/{1}\n".format(self.data_path, veglib))
        fout.write("VEGPARAM\t{0}/{1}\n".format(self.data_path, vegparam))
        fout.write("VEGPARAM_LAI\tTRUE\n")
        fout.write("ROOT_ZONES\t{0:d}\n".format(root_zones))
        fout.write("LAI_SRC\tLAI_FROM_VEGPARAM\n")
        nbands = self._getSnowbands(snowbands)
        fout.write(
            "SNOW_BAND\t{0:d}\t{1}/{2}\n".format(nbands, self.data_path, snowbands))
        fout.write("RESULT_DIR\t{0}/output\n".format(self.model_path))
        fout.write("OUT_STEP\t24\n")
        fout.write("BINARY_OUTPUT\tFALSE\n")
        fout.write("MOISTFRACT\tFALSE\n")
        fout.write(
            "COMPRESS\tFALSE\nALMA_OUTPUT\tFALSE\nPTR_HEADER\tFALSE\nPRT_SNOW_BAND\tFALSE\n")
        fout.write(vicoutput.template(["eb", "wb", "sub", "sur", "csp", "eva"]))
        fout.close()

    def createIndexTable(self, dataset):
        """Creates index table from raster row, column, and tile for each grid cell."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        sname, tname = dataset.split(".")
        rtable = "{0}.{1}_{2}".format(sname, tname, int(1.0 / self.res))
        dbio.deleteTable(self.dbname, "public", "{0}_xy".format(sname))
        sql = "create table {0}_xy as (select gid,st_worldtorastercoordx(rast,geom) as x,st_worldtorastercoordy(rast,geom) as y,rid as tile from {4},{5}.basin where fdate=date'{1}-{2}-{3}' and st_intersects(rast,geom))".format(
            sname, self.startyear, self.startmonth, self.startday, rtable, self.name)
        cur.execute(sql)
        db.commit()
        cur.execute("create index {0}_xy_r on {0}_xy(tile)".format(sname))
        cur.close()
        db.close()
        return rtable

    def _getTiles(self, itable):
        """Get raster tile IDs for the domain."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        cur.execute("select distinct(tile) from {0}".format(itable))
        tiles = [int(r[0]) for r in cur.fetchall()]
        cur.close()
        db.close()
        return tiles

    def _dropIndexTable(self, sname):
        """Deletes index table."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        cur.execute("drop table {0}_xy".format(sname))
        db.commit()
        cur.close()
        db.close()

    def _getTileData(self, rtable, t):
        """Retrieve data from *rtable* for specific tile *t*."""
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        var = rtable.split(".")[0]
        sql = "select gid,fdate,st_value(rast,x,y) from {0},{1}_xy where rid=tile and tile={8} and fdate>=date'{2}-{3}-{4}' and fdate<=date'{5}-{6}-{7}' order by gid,fdate".format(
            rtable, var, self.startyear, self.startmonth, self.startday, self.endyear, self.endmonth, self.endday, t)
        cur.execute(sql)
        data = cur.fetchall()
        return data

    def getForcings(self, options):
        """Get meteorological forcings from database."""
        log = logging.getLogger(__name__)
        if not ('precip' in options and 'temperature' in options and 'wind' in options):
            log.error("No data source provided for VIC forcings")
            sys.exit()
        datasets = ["precip." + options["precip"], "tmax." + options["temperature"],
                    "tmin." + options["temperature"], "wind." + options["wind"]]
        if 'lai' in options:
            datasets.append("lai" + options["lai"])
            self.lai = options['lai']
        options['tmax'] = options['temperature']
        options['tmin'] = options['temperature']
        rtables = {}
        for v in ['precip', 'tmax', 'tmin', 'wind']:
            rtables[v] = self.createIndexTable("{0}.{1}".format(v, options[v]))
        tiles = {v: self._getTiles("{0}_xy".format(v))
                 for v in ['precip', 'tmax', 'tmin', 'wind']}
        data = {}
        nprocs = mp.cpu_count()
        p = mp.Pool(nprocs)
        for v in tiles:
            reader = TileReader(self.dbname, rtables[
                                v], self.startyear, self.startmonth, self.startday, self.endyear, self.endmonth, self.endday)
            data[v] = p.map_async(reader, tiles[v])
        data = {v: [i for s in data[v].get() for i in s if i[2] is not None] for v in data}
        p.close()
        p.join()
        for s in data:
            self._dropIndexTable(s)
        self.precip = options['precip']
        self.temp = options['temperature']
        self.wind = options['wind']
        return data['precip'], data['tmax'], data['tmin'], data['wind']

    def _handleMissingData(self, data, table):
        """Handle missing data when running the VIC model. Most likley the missing
        data have to do with the windowing in the spatial queries. Therefore, we will
        query the entire raster for the pixels with missing data."""
        log = logging.getLogger(__name__)
        t0 = date(self.startyear, self.startmonth, self.startday)
        t1 = date(self.endyear, self.endmonth, self.endday)
        try:
            ndays = (t1 - t0).days + 1
            assert len(data) == len(self.lat) * ndays
        except AssertionError:
            log.warning("Missing meteorological data for {0} in database. Filling with nearest values!".format(table))
            # Check for missing data due to raster tiling
            pixels = list(set(self.gid.keys()) - set([r[0] for r in data]))
            db = dbio.connect(self.dbname)
            cur = db.cursor()
            for p in pixels:
                sql = "select gid, fdate, st_nearestvalue(rast,geom) from {0},{1}.basin where gid={2} and fdate>=date'{3}' and fdate<=date'{4}' order by fdate".format(table, self.name, p, t0.strftime("%Y-%m-%d"), t1.strftime("%Y-%m-%d"))
                cur.execute(sql)
                res = cur.fetchall()
                data += res
            cur.close()
            db.close()
        # Now check for missing values in time series of each pixel
        try:
            assert len(data) == len(self.lat) * ndays
        except AssertionError:
            pdata = OrderedDict({d[0]: [] for d in data})
            for d in data:
                pdata[d[0]].append(d)
            data = {}
            for gid in pdata:
                dates = pd.to_datetime([d[1] for d in pdata[gid]], infer_datetime_format=True)
                p = pd.Series([d[2] for d in pdata[gid]], dates)
                if len(p) < ndays:
                    p = p.reindex(pd.date_range(self.startdate, self.enddate)).interpolate(method='pad')
                for t, v in p.iteritems():
                    data.append([gid, t.date(), v])
        return data

    def writeForcings(self, prec, tmax, tmin, wind, lai=None):
        """Write VIC meteorological forcing data files."""
        log = logging.getLogger(__name__)
        if not os.path.exists(self.model_path + '/forcings'):
            os.mkdir(self.model_path + '/forcings')
        if all([v in self.__dict__ for v in ['precip', 'tmax', 'tmin', 'wind']]):
            prec = self._handleMissingData(prec, "precip.{0}".format(self.precip))
            tmax = self._handleMissingData(tmax, "tmax.{0}".format(self.temp))
            tmin = self._handleMissingData(tmin, "tmin.{0}".format(self.temp))
            wind = self._handleMissingData(wind, "wind.{0}".format(self.wind))
        def met2df(data):
            pdata = {}
            dt = {}
            for i in range(len(data)):
                k = data[i][0]
                if k in pdata:
                    pdata[k].append(data[i][2])
                    dt[k].append(data[i][1])
                else:
                    pdata[k] = [data[i][2]]
                    dt[k] = [data[i][1]]
            for k in pdata:
                pdata[k] = pd.Series(pdata[k], dt[k])
            return pd.DataFrame(pdata)
        prec = met2df(prec)
        tmax = met2df(tmax)
        tmin = met2df(tmin)
        wind = met2df(wind)
        for gid in prec.columns:
            filename = "data_{0:.{2}f}_{1:.{2}f}".format(
                    self.gid[gid][0], self.gid[gid][1], self.grid_decimal)
            out = pd.concat([prec[gid], tmax[gid], tmin[gid], wind[gid]], axis=1).values
            with open("{0}/forcings/{1}".format(self.model_path, filename), 'w') as fout:
                for t in range(out.shape[0]):
                    fout.write("{0:f} {1:.2f} {2:.2f} {3:.1f}\n".format(out[t, 0], out[t, 1], out[t, 2], out[t, 3]))

    def run(self, vicexec):
        """Run VIC model."""
        log = logging.getLogger(__name__)
        log.info("Running VIC...")
        if not os.path.exists(self.model_path + '/output'):
            os.mkdir(self.model_path + '/output')
        print(vicexec, self.model_path)
        proc = subprocess.Popen([vicexec, "-g", "{0}/global.txt".format(self.model_path)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        with proc.stdout:
            for line in iter(proc.stdout.readline, b''):
                log.debug(line.strip())

    def getOutputStruct(self, globalfile):
        """Creates a dictionary with output variable-file pairs."""
        fin = open(globalfile)
        prefix = None
        skipyear = 0
        c = 3  # Assumes daily output
        out = {}
        for line in fin:
            if line.find("OUTFILE") == 0:
                prefix = line.split()[1]
                c = 3
            elif line.find("SKIPYEAR") == 0:
                skipyear = int(line.split()[1])
            else:
                if len(line) > 1 and line[0] != "#" and prefix:
                    varname = line.split()[1].replace("OUT_", "").lower()
                    out[varname] = ("output/" + prefix, c)
                    if varname in ["soil_moist", "soil_temp", "smliqfrac", "smfrozfrac"]:
                        c += self.nlayers
                    else:
                        c += 1
        fin.close()
        out['tmax'] = ("forcings/data", 1)
        out['tmin'] = ("forcings/data", 2)
        out['rainf'] = ("forcings/data", 0)
        self.skipyear = skipyear
        return out

    def saveToDB(self, args, initialize=True, skipsave=0, ensemble=0):
        """Reads VIC output for selected variables."""
        log = logging.getLogger(__name__)
        droughtvars = vicoutput.droughtVariables()
        layervars = ["soil_moist", "soil_temp", "smliqfrac", "smfrozfrac"]
        outvars = self.getOutputStruct(self.model_path + "/global.txt")
        outdata = {}
        if self.lat and self.lon:
            nrows = int(np.round((max(self.lat) - min(self.lat)) / self.res) + 1)
            ncols = int(np.round((max(self.lon) - min(self.lon)) / self.res) + 1)
            mask = np.zeros((nrows, ncols), dtype='bool')
            nt = (date(self.endyear, self.endmonth, self.endday) -
                  date(self.startyear + self.skipyear, self.startmonth, self.startday)).days + 1
            args = vicoutput.variableGroup(args)
            if args:
                for var in args:
                    if var in outvars or var in droughtvars:
                        if var in layervars:
                            outdata[var] = np.zeros((nt, self.nlayers, nrows, ncols)) + self.nodata
                        else:
                            outdata[var] = np.zeros((nt, 1, nrows, ncols)) + self.nodata
                    else:
                        log.warning("Variable {0} not found in output files. Skipping import.".format(var))
                prefix = set([outvars[v][0] for v in outdata.keys() if v not in droughtvars])
                startdate = "{0}-{1}-{2}".format(self.startyear, self.startmonth, self.startday)
                enddate = "{0}-{1}-{2}".format(self.endyear, self.endmonth, self.endday)
                dates = pd.date_range(startdate, enddate).values
                for c in range(len(self.lat)):
                    pdata = {}
                    for p in prefix:
                        filename = "{0}/{1}_{2:.{4}f}_{3:.{4}f}".format(self.model_path, p, self.lat[c], self.lon[c], self.grid_decimal)
                        print(filename)
                        pdata[p] = pd.read_csv(filename, delim_whitespace=True, header=None).values
                    i = int((max(self.lat) + self.res / 2.0 - self.lat[c]) / self.res)
                    j = int((self.lon[c] - min(self.lon) + self.res / 2.0) / self.res)
                    mask[i, j] = True
                    for v in [v for v in outdata if v not in droughtvars]:
                        if v in layervars:
                            for lyr in range(self.nlayers):
                                outdata[v][:, lyr, i, j] = pdata[outvars[v][0]][:, outvars[v][1] + lyr]
                        else:
                            outdata[v][:, 0, i, j] = pdata[outvars[v][0]][:, outvars[v][1]]
                    log.info("Read output for {0}|{1}".format(self.lat[c], self.lon[c]))
                for var in args:
                    if var not in droughtvars:
                        if outdata[var] is not None:
                            self.writeToDB(outdata[var], dates, "{0}".format(var), initialize, skipsave=skipsave)
                del outdata
                for var in args:
                    if var in droughtvars:
                        outdata = np.zeros((nt, 1, nrows, ncols)) + self.nodata
                        dout = drought.calc(var, self, ensemble)
                        if dout is not None:
                            outdata[:, 0, :, :] = dout
                            self.writeToDB(outdata, dates, "{0}".format(var), initialize, skipsave=skipsave)
        else:
            log.info("No pixels simulated, not saving any output!")

    def _writeRaster(self, data, filename):
        """Writes GeoTIFF raster temporarily so that it can be imported into the database."""
        nrows, ncols = data.shape
        driver = gdal.GetDriverByName("GTiff")
        ods = driver.Create(filename, ncols, nrows, 1, gdal.GDT_Float32, ["COMPRESS=LZW"])
        ods.SetGeoTransform([min(self.lon) - self.res / 2.0, self.res,
                             0, max(self.lat) + self.res / 2.0, 0, -self.res])
        srs = osr.SpatialReference()
        srs.SetWellKnownGeogCS("WGS84")
        ods.SetProjection(srs.ExportToWkt())
        ods.GetRasterBand(1).WriteArray(data)
        ods.GetRasterBand(1).SetNoDataValue(self.nodata)
        ods = None

    def writeToDB(self, data, dates, tablename, initialize, ensemble=False, skipsave=0):
        """Writes output data into database."""
        log = logging.getLogger(__name__)
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        if dbio.tableExists(self.dbname, self.name, tablename) and not dbio.columnExists(self.dbname, self.name, tablename, "ensemble"):
            cur.execute("alter table {0}.{1} add column ensemble int".format(self.name, tablename))
            cur.execute("update {0}.{1} set ensemble=0".format(self.name, tablename))
            db.commit()
        if dbio.tableExists(self.dbname, self.name, tablename):
            if initialize:
                for dt in [self.startdate + timedelta(t) for t in range((self.enddate - self.startdate).days+1)]:
                    dbio.deleteRasters(self.dbname, "{0}.{1}".format(self.name, tablename), dt)
        else:
            sql = "create table {0}.{1} (rid serial not null primary key, fdate date not null, rast raster, filename text)".format(self.name, tablename)
            #sql = "create table {0}.{1} (id serial not null primary key, rid int not null, fdate date not null, rast raster, filename text)".format(
            #    self.name, tablename)
            cur.execute(sql)
            if data.shape[1] > 1:
                cur.execute("alter table {0}.{1} add column layer int".format(self.name, tablename))
            # always add ensemble column (0 for deterministic run)
            cur.execute("alter table {0}.{1} add column ensemble int".format(self.name, tablename))
            db.commit()
        startyear, startmonth, startday = self.startyear, self.startmonth, self.startday
        if skipsave > 0:
            ts = date(self.startyear, self.startmonth,
                      self.startday) + timedelta(skipsave)
            data = data[skipsave:]
            startyear, startmonth, startday = ts.year, ts.month, ts.day
        tiffiles = []
        for t in range(data.shape[0]):
            dt = date(startyear, startmonth, startday) + timedelta(t)
            for lyr in range(data.shape[1]):
                filename = "{0}/{1}_{2}{3:02d}{4:02d}_{5:02d}_{6:02d}.tif".format(
                    os.path.abspath(self.dbpath), tablename, dt.year, dt.month, dt.day, lyr + 1, int(ensemble))
                self._writeRaster(data[t, lyr, :, :], filename)
                tiffiles.append(filename)
        # Write rasters to database in chunks of 365 days to avoid "argument list too long" error
        with tempfile.NamedTemporaryFile('a') as fout:
            subprocess.call(["raster2pgsql", "-R", "-s", "4326", "-F", "-d"] + tiffiles[:min(365, len(tiffiles))] + ["temp"], stdout=fout)
            for i in range(1, int(len(tiffiles) / 365)+1):
                fnames = tiffiles[i*365:min((i+1)*365, len(tiffiles))]
                fout.flush()
                subprocess.call(["raster2pgsql", "-R", "-s", "4326", "-F", "-a"] + fnames + ["temp"], stdout=fout)
            cmd = ["psql", "-d", self.dbname, "-f", fout.name]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            for line in iter(proc.stdout.readline, b''):
                log.debug(line.strip())
            proc.wait()
        cur.execute("alter table temp add column fdate date")
        cur.execute("update temp set fdate = date (concat_ws('-',substring(filename from {0} for 4),substring(filename from {1} for 2),substring(filename from {2} for 2)))".format(
            len(tablename) + 2, len(tablename) + 6, len(tablename) + 8))
        if data.shape[1] > 1:
            cur.execute("alter table temp add column layer int")
            cur.execute("update temp set layer=(substring(filename from {0} for 2))::int".format(
                len(tablename) + 11))
        if data.shape[1] > 1:
            cur.execute("insert into {0}.{1} (fdate,layer,rast,filename) select fdate,layer,rast,filename from temp".format(self.name, tablename))
        else:
            cur.execute("insert into {0}.{1} (fdate,rast,filename) select fdate,rast,filename from temp".format(self.name, tablename))
        sql = "update {0}.{1} set ensemble = {2} where ensemble is null"
        if bool(ensemble):
            cur.execute(sql.format(self.name, tablename, int(ensemble)))
        else:
            cur.execute(sql.format(self.name, tablename, 0))
        cur.execute("drop index if exists {0}.{1}_dtidx".format(
            self.name, tablename))
        cur.execute("create index {1}_dtidx on {0}.{1}(fdate)".format(
            self.name, tablename))
        cur.execute("drop index if exists {0}.{1}_spidx".format(
            self.name, tablename))
        cur.execute("create index {1}_spidx on {0}.{1} using gist(st_convexhull(rast))".format(
            self.name, tablename))
        db.commit()
        cur.close()
        db.close()

    def save(self, saveto, args, initialize=True, skipsave=0, ensemble=0):
        """Reads and saves selected output data variables into the database or a user-defined directory."""
        if saveto == "db":
            self.saveToDB(args, initialize=initialize, skipsave=skipsave, ensemble=ensemble)
        else:
            #if initialize:
                #if os.path.isdir(saveto):
                #    shutil.rmtree(saveto)
            #     elif os.path.isfile(saveto):
            #         os.remove(saveto)
            # os.makedirs(saveto)
            # shutil.move(self.model_path+"/output", saveto)
            # shutil.move(self.model_path+"/forcings", saveto)
            shutil.copytree(self.model_path, saveto,
                            ignore=shutil.ignore_patterns("*.txt"))
                            
    def _findExecutable(exefile, path=None):
        """Find if *exefile* is available. """
        if path is None:
            path = os.environ['PATH']
        paths = path.split(os.pathsep)
        extlist = ['']
        for ext in extlist:
            execname = exefile + ext
            if os.path.isfile(execname):
                return execname
            else:
                for p in paths:
                    f = os.path.join(p, execname)
                    if os.path.isfile(f):
                       return f
        else:
            return None

    def rout(self,routfilepath):
        """
        Performs Lohman Routing on the resultant VIC runoff and baseflow
        by defulat VIC and routing simulation duration is same
           
        ADDED BY VIKALP MISHRA 
        AUG 2021
        SERVIR SCO
        """
        print('VIC Routing')
        routexe = 'external/vic/src/rout'
        rout_conf_path = routfilepath
        mod_path = self.model_path.split('.')[-1]
        mod_path = '/home/rheas2/RHEAS/'+mod_path
        rout_out_path = mod_path+'/rout/'
        rout_flx_path = mod_path+'/output/fluxes_'
        
        rcnf = open(rout_conf_path+'/rout.conf','w')
        rcnf.write('#INPUT FILE FOR THE {0} BASIN\n'.format(self.name))
        rcnf.write('#PATH OF THE FLOW DIRECTION FILE\n')
        rcnf.write(rout_conf_path + '/fdr.txt\n')
        rcnf.write('#NAME OF FLOW VELOCITY FILE\n')
        rcnf.write('.false.\n')
        rcnf.write('.35\n')
        rcnf.write('#NAME OF DIFFUSION FILE\n')
        rcnf.write('.false.\n')
        rcnf.write('1200\n')
        rcnf.write('#NAME OF XMASK FILE\n')
        rcnf.write('.false.\n')
        rcnf.write('7000\n')
        rcnf.write('#NAME OF FRACTION FILE\n')
        rcnf.write('.true.\n')
        rcnf.write(rout_conf_path + '/fract.asc\n')
        rcnf.write('#NAME OF STATION FILE\n')
        rcnf.write(rout_conf_path+'/StationLoc.txt\n')
        rcnf.write('#PATH TO FORCING FILES & its precision digits\n')
        rcnf.write('{0}\n'.format(rout_flx_path))
        rcnf.write('4\n')
        rcnf.write('#PATH TO OUTPUT DIRECTORY\n')
        rcnf.write('{0}\n'.format(rout_out_path))
        rcnf.write('#START AND END DATES (VIC & ROUTING) - START(YYYY M) END(YYYY MM)\n')
        rcnf.write('{0} {1} {2} {3}\n'.format(self.startyear, self.startmonth,self.endyear,self.endmonth))
        rcnf.write('{0} {1} {2} {3}\n'.format(self.startyear, self.startmonth,self.endyear,self.endmonth)) 
        rcnf.write('#NAME OF UNIT HYDROGRAPH FILE\n')
        rcnf.write(rout_conf_path+'/uh.txt\n')
        rcnf.close()


        if not os.path.exists(rout_out_path): 
            os.makedirs(rout_out_path)   
        
        os.system(routexe+' '+rout_conf_path+'rout.conf')

    def run_parallel(self, vicexec, ncores):
        """Run VIC model in parallel."""
        nprocs = mp.cpu_count()
        try:
            ncores = min(int(ncores), nprocs)
        except ValueError:
            ncores = nprocs
        log = logging.getLogger(__name__)
        log.info("Running VIC in Parallel...")
        if not os.path.exists(self.model_path + '/output'):
            os.mkdir(self.model_path + '/output')
        # Create different soil files
        with open(os.path.join(self.model_path, 'soil.txt'), 'r') as f:
            lines_soil = f.readlines()
        chunk_size = len(lines_soil) // ncores
        break_idx = [chunk_size*i for i in range(ncores)]
        break_idx.append(len(lines_soil))
        for core, (init, end) in enumerate(zip(break_idx[:-1], break_idx[1:])):
            with open(os.path.join(self.model_path, f'soil_{core}.txt'), 'w') as f:
                f.writelines(lines_soil[init:end])
        # Create different global files. Soil file and state must be different
        with open(os.path.join(self.model_path, 'global.txt'), 'r') as f:
            lines_global = f.readlines()
        for core in range(ncores):
            with open(os.path.join(self.model_path, f'global_{core}.txt'), 'w') as f:
                for line in lines_global:
                    if 'SOIL\t' in line[:5]:
                        line = f'SOIL\t{self.model_path}/soil_{core}.txt\n'
                    if 'STATENAME\t' in line[:10]:
                        line = f'STATENAME\tstate/vic.state_{core}\n'
                    if 'INIT_STATE\t' in line[:12]:
                        state_file = line.split('\t')[-1]
                        state_file = f'_{core}_'.join(state_file.split('_'))
                        state_file = state_file.replace('\n', '')
                        line = 'INIT_STATE\t' + state_file + '\n'
                        assert os.path.exists(state_file), \
                            f'{state_file} does not exist, there must be one state file per core'
                    f.write(line)
        proc_list = []
        for core in range(ncores):
            proc_args = [vicexec, "-g", "{0}/global_{1}.txt".format(self.model_path, core)]
            proc_list.append(
                subprocess.Popen(proc_args)
            )
        exit_codes = [p.wait() for p in proc_list]
        assert all([i == 0 for i in exit_codes]), "VIC model failed"

