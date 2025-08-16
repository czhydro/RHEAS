""" RHEAS module for retrieving meteorological forecasts/hindcasts
from the NMME model suite.

.. module:: nmme
   :synopsis: Retrieve NMME meteorological forecast data

.. moduleauthor:: Kostas Andreadis <kandread@umass.edu>

"""

import logging
import os
import random
import shutil
import string
import subprocess
import sys
import tempfile
import zipfile
import climateserv.api
from datetime import datetime, timedelta

import numpy as np

from .. import dbio
from . import datasets

table = "precip.nmme"


def dates(dbname):
    dts = datasets.dates(dbname, table)
    return dts


def _setEnsemble(dbname, sname, ens):
    """Set ensemble column in NMME data table."""
    db = dbio.connect(dbname)
    cur = db.cursor()
    cur.execute("select distinct(resolution) from vic.soils")
    if bool(cur.rowcount):
        resolutions = [r[0] for r in cur.fetchall()]
        for res in resolutions:
            ext = int(1.0 / res)
            if not dbio.columnExists(dbname, sname, "nmme_{0}".format(ext), "ensemble"):
                cur.execute("alter table {0}.nmme_{1} add column ensemble int".format(sname, ext))
                db.commit()
            sql = "update {0}.nmme_{1} set ensemble = {2} where ensemble is null".format(sname, ext, ens)
            cur.execute(sql)
    db.commit()
    cur.close()
    db.close()


def ingest(dbname, varname, filename, dt, ens):
    """Imports Geotif *filename* into database *dbname*."""
    schema = {'Precipitation': 'precip', 'Temperature': 'tmax'}
    db = dbio.connect(dbname)
    cur = db.cursor()
    cur.execute(
        "select * from information_schema.tables where table_schema='{0}' and table_name='nmme'".format(schema[varname]))
    if not bool(cur.rowcount):
        cur.execute("create table {0}.nmme (rid serial not null primary key, fdate date, ensemble int, rast raster)".format(
            schema[varname]))
        db.commit()
    cur.execute("select * from {0}.nmme where fdate='{1}' and ensemble={2}".format(schema[varname], dt.strftime("%Y-%m-%d"), ens))
    if bool(cur.rowcount):
        cur.execute("delete from {0}.nmme where fdate='{1}' and ensemble={2}".format(schema[varname], dt.strftime("%Y-%m-%d"), ens))
        db.commit()
    dbio.ingest(dbname, filename, dt, "{0}.nmme".format(schema[varname]), False, False)
    sql = "update {0}.nmme set ensemble = {1} where ensemble is null".format(schema[varname], ens)
    cur.execute(sql)
    db.commit()
    cur.execute("select distinct(resolution) from vic.soils")
    if bool(cur.rowcount):
        resolutions = [r[0] for r in cur.fetchall()]
        for res in resolutions:
            if dbio.tableExists(dbname, schema[varname], "nmme_{}".format(int(1.0 / res))):
                cur.execute("select * from {0}.nmme_{1} where fdate='{2}' and ensemble={3}".format(schema[varname], int(1.0 / res), dt.strftime("%Y-%m-%d"), ens))
                if bool(cur.rowcount):
                    cur.execute("delete from {0}.nmme_{1} where fdate='{2}' and ensemble={3}".format(schema[varname], int(1.0 / res), dt.strftime("%Y-%m-%d"), ens))
                    db.commit()
    tilesize = (10, 10)
    dbio.createResampledTables(dbname, schema[varname], "nmme", dt, tilesize, False, "and ensemble={0}".format(ens))
    _setEnsemble(dbname, schema[varname], ens)
    cur.close()
    db.close()


def download(dbname, dts, bbox=None):
    """Downloads NMME ensemble forecast data from the SERVIR ClimateSERV
    data server, and imports them into the database *dbname*. Optionally uses
    a bounding box to limit the region with [minlon, minlat, maxlon, maxlat]."""
    log = logging.getLogger(__name__)
    nens = 34
    modelname = lambda e: "CCSM4" if e < 10 else "CFSV2"
    ens_num = lambda e: e+1 if e < 10 else e - 9
    varnames = ["Precipitation", "Temperature"]
    outpath = tempfile.mkdtemp()
    if bbox:
        coords = [[bbox[0], bbox[1]], [bbox[2], bbox[1]], [bbox[2], bbox[3]], [bbox[0], bbox[3]], [bbox[0], bbox[1]]]
    else:
        coords = [[-180, -90], [180, -90], [180, 90], [-180, 90], [-180, -90]]
    dt1 = dts[0].strftime("%m/%d/%Y")
    dt2 = dts[1].strftime("%m/%d/%Y")
    for varname in varnames:
        for e in range(nens):
            zfilename = "{0}/{1}_{2}.zip".format(outpath, varname, e+1)
            climateserv.api.request_data(modelname(e), "Download", dt1, dt2, coords, "ens{0:02d}".format(ens_num(e)), varname, zfilename)
            f = zipfile.ZipFile(zfilename)
            filenames = [name for name in f.namelist() if name.endswith(".tif")]
            f.extractall(outpath, filenames)
            for filename in filenames:
                dt = datetime.strptime(filename.split("_")[-1][0:-4], "%Y%m%d")
                if varname == "Temperature":
                    # convert from Kelvin to Celsius
                    proc = subprocess.Popen(["gdal_calc.py", "-A", "{0}/{1}".format(outpath, filename), "--calc=A-273.15", "--outfile={0}/C{1}".format(outpath, filename)])
                    out, err = proc.communicate()
                    log.debug(out)
                    filename = "C" + filename
                ingest(dbname, varname, "{0}/{1}".format(outpath, filename), dt, e+1)
    shutil.rmtree(outpath)


def _queryDataset(dbname, tablename, name, startyear, startmonth, startday, endyear, endmonth, endday, ens=None):
    """Retrieve meteorological forcing dataset from database."""
    temptable = ''.join(random.SystemRandom().choice(string.ascii_letters) for _ in range(8))
    if ens is None:
        sql = "create table {0}_xy as (select gid,st_worldtorastercoordx(rast,geom) as x,st_worldtorastercoordy(rast,geom) as y,rid as tile from {4},{5}.basin where fdate=date'{1}-{2}-{3}' and st_intersects(rast,geom))".format(temptable, startyear, startmonth, startday, tablename, name)
    else:
        sql = "create table {0}_xy as (select gid,st_worldtorastercoordx(rast,geom) as x,st_worldtorastercoordy(rast,geom) as y,rid as tile from {4},{5}.basin where fdate=date'{1}-{2}-{3}' and st_intersects(rast,geom) and ensemble={6})".format(temptable, startyear, startmonth, startday, tablename, name, ens)
    db = dbio.connect(dbname)
    cur = db.cursor()
    cur.execute(sql)
    cur.execute("create index {0}_xy_r on {0}_xy(tile)".format(temptable))
    db.commit()
    if ens is None:
        sql = "select gid,fdate,st_nearestvalue(rast,x,y) from {0},{1}_xy where rid=tile and fdate>=date'{2}-{3}-{4}' and fdate<=date'{5}-{6}-{7}' order by gid,fdate".format(tablename, temptable, startyear, startmonth, startday, endyear, endmonth, endday)
    else:
        sql = "select gid,fdate,st_nearestvalue(rast,x,y) from {0},{1}_xy where rid=tile and fdate>=date'{2}-{3}-{4}' and fdate<=date'{5}-{6}-{7}' and ensemble={8} order by gid,fdate".format(tablename, temptable, startyear, startmonth, startday, endyear, endmonth, endday, ens)
    cur.execute(sql)
    data = [r for r in cur.fetchall()]
    cur.execute("drop table {0}_xy".format(temptable))
    db.commit()
    cur.close()
    db.close()
    return data


def _getForcings(options, models, res, e):
    """Retrieve meteorological forcings for ensemble."""
    nens = len(models)
    db = dbio.connect(models.dbname)
    cur = db.cursor()
    rtables = dbio.getResampledTables(models.dbname, options, res)
    rsmp = rtables['precip'].split("_")[1]
    prec = _queryDataset(models.dbname, "precip.nmme_{0}".format(rsmp), models.name, models.startyear, models.startmonth, models.startday, models.endyear, models.endmonth, models.endday, e+1)
    temp = _queryDataset(models.dbname, "tmax.nmme_{0}".format(rsmp), models.name, models.startyear, models.startmonth, models.startday, models.endyear, models.endmonth, models.endday, e+1)
    sql = "select distinct(date_part('year',fdate)) from tmax.{0}".format(rtables['tmax'])
    cur.execute(sql)
    years = [r[0] for r in cur.fetchall()]
    if len(years) > 2:
        years.remove(min(years))
        years.remove(max(years))
    if len(years) > 0:
        ndays = (datetime(models.endyear, models.endmonth, models.endday) - datetime(models.startyear, models.startmonth, models.startday)).days
        yr = int(np.random.choice(years))
        t0 = datetime(yr, models.startmonth, models.startday)
        t1 = t0 + timedelta(ndays)
        vtmax = _queryDataset(models.dbname, "tmax.{0}".format(rtables['tmax']), models.name, t0.year, t0.month, t0.day, t1.year, t1.month, t1.day)
        vtmin = _queryDataset(models.dbname, "tmin.{0}".format(rtables['tmin']), models.name, t0.year, t0.month, t0.day, t1.year, t1.month, t1.day)
        wind = _queryDataset(models.dbname, "wind.{0}".format(rtables['wind']), models.name, t0.year, t0.month, t0.day, t1.year, t1.month, t1.day)
        tmax = [(vtmax[i][0], vtmax[i][1], temp[i][2] + 0.5 * (vtmax[i][2] - vtmin[i][2])) for i in range(len(vtmax))]
        tmin = [(vtmin[i][0], vtmin[i][1], temp[i][2] - 0.5 * (vtmax[i][2] - vtmin[i][2])) for i in range(len(vtmin))]
    else:
        prec = tmax = tmin = wind = None
    return prec, tmax, tmin, wind


def generate(options, models):
    """Generate meteorological forecast forcings from downscaled NMME data."""
    log = logging.getLogger(__name__)
    options['vic']['tmax'] = options['vic']['temperature']
    options['vic']['tmin'] = options['vic']['temperature']
    db = dbio.connect(models.dbname)
    cur = db.cursor()
    dt0 = datetime(models.startyear, models.startmonth, models.startday)
    dt1 = datetime(models.endyear, models.endmonth, models.endday)
    # check if forecast period exists in NMME data
    sql = "select count(distinct(fdate)) from precip.nmme where fdate>=date'{0}' and fdate<=date'{1}'".format(dt0.strftime("%Y-%m-%d"), dt1.strftime("%Y-%m-%d"))
    cur.execute(sql)
    ndata = cur.fetchone()[0]
    if ndata == (dt1 - dt0).days + 1:
        for e in range(len(models)):
            prec, tmax, tmin, wind = _getForcings(options, models, models.res, e)
            if tmax is None or tmin is None or wind is None:
                log.error("No data found to generate VIC forcings for NMME forecast. Exiting...")
                sys.exit()
            else:
                models[e].precip = options['vic']['precip']
                models[e].temp = options['vic']['temperature']
                models[e].wind = options['vic']['wind']
                models[e].writeForcings(prec, tmax, tmin, wind)
    else:
        log.error("Not enough data found for requested forecast period! Exiting...")
        sys.exit()
    cur.close()
    db.close()
