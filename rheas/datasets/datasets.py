""" Definition for RHEAS Datasets package.

.. module:: datasets
   :synopsis: Definition of the Datasets package

.. moduleauthor:: Kostas Andreadis <kandread@umass.edu>

"""

import gzip
import logging
import os
import sys
import zipfile
from configparser import ConfigParser
from datetime import datetime, timedelta
from functools import partial

import numpy as np

from .. import dbio
from .decorators import geotiff, path


def uncompress(filename, outpath):
    """Uncompress archived files."""
    if filename.endswith("gz"):
        f = gzip.open("{0}/{1}".format(outpath, filename), 'rb')
        contents = f.read()
        f.close()
        lfilename = filename.replace(".gz", "")
        with open("{0}/{1}".format(outpath, lfilename), 'wb') as f:
            f.write(contents)
    elif filename.endswith("zip"):
        f = zipfile.ZipFile("{0}/{1}".format(outpath, filename))
        lfilename = filter(lambda s: s.endswith("tif"), f.namelist())[0]
        f.extract(lfilename, outpath)
    else:
        lfilename = filename
    return lfilename


def readDatasetList(filename):
    """Read list of datasets to be fetched and imported into
    the RHEAS database."""
    log = logging.getLogger(__name__)
    conf = ConfigParser()
    try:
        conf.read(filename)
    except FileNotFoundError:
        log.error("File not found: {}".format(filename))
        sys.exit()
    return conf


def dates(dbname, tablename):
    """Check what dates need to be imported for a specific dataset."""
    dts = None
    db = dbio.connect(dbname)
    cur = db.cursor()
    sname, tname = tablename.split(".")
    cur.execute(
        "select * from information_schema.tables where table_name='{0}' and table_schema='{1}'".format(tname, sname))
    if bool(cur.rowcount):
        sql = "select max(fdate) from {0}".format(tablename)
        cur.execute(sql)
        te = cur.fetchone()[0]
        try:
            te = datetime(te.year, te.month, te.day)
            if te < datetime.today():
                dts = (te + timedelta(1), datetime.today())
        except (AssertionError, AttributeError):
            dts = None
    else:
        dts = None
    return dts


def spatialSubset(lat, lon, res, bbox):
    """Subsets arrays of latitude/longitude based on bounding box *bbox*."""
    log = logging.getLogger(__name__)
    if bbox is None:
        i1 = 0
        i2 = len(lat)-1
        j1 = 0
        j2 = len(lon)-1
    else:
        i1 = np.where(bbox[3] <= lat+res/2)[0]
        i2 = np.where(bbox[1] >= lat-res/2)[0]
        j1 = np.where(bbox[0] >= lon-res/2)[0]
        j2 = np.where(bbox[2] <= lon+res/2)[0]
        if any([len(idx) <= 0 for idx in [i1, i2, j1, j2]]):
            log.warning("Requested domain larger than dataset extent. Ingesting entire raster.")
        i1 = i1[-1] if len(i1) > 0 else 0
        i2 = i2[0] if len(i2) > 0 else len(lat)-1
        j1 = j1[-1] if len(j1) > 0 else 0
        j2 = j2[0] if len(j2) > 0 else len(lon)-1
    return i1, i2+1, j1, j2+1


def download(dbname, dts, bbox, conf, name):
    """Download a generic dataset based on user-provided information."""
    log = logging.getLogger(__name__)
    nens = conf.getint(name, 'ensemble size', fallback=1)
    try:
        url = conf.get(name, 'path')
        res = conf.getfloat(name, 'res')
        table = conf.get(name, 'table')
    except:
        url = res = table = None

    @geotiff
    @path
    def fetch(dbname, dt, bbox):
        return url, bbox, dt
    base_url = url
    sname, tname = table.split(".")
    for e in range(nens):
        url = partial(base_url.format, ens=e+1)
        if url is not None and res is not None and table is not None:
            for dt in [dts[0] + timedelta(tt) for tt in range((dts[-1] - dts[0]).days + 1)]:
                data, lat, lon, t = fetch(dbname, dt, bbox)
                if dbio.tableExists(dbname, sname, tname):
                    if dbio.columnExists(dbname, sname, tname, "ensemble"):
                        dbio.deleteRasters(dbname, table, t, "and ensemble={0}".format(e+1))
                    else:
                        dbio.deleteRasters(dbname, table, t)
                ingest(dbname, table, data, lat, lon, res, t, overwrite=False)
        else:
            log.warning("Missing options for local dataset {0}. Nothing ingested!".format(name))
        dbio.setAttribute(dbname, table, "ensemble", e+1)


def ingest(dbname, table, data, lat, lon, res, t, resample=True, overwrite=True):
    """Import data into RHEAS database."""
    log = logging.getLogger(__name__)
    if data is not None:
        if len(data.shape) > 2:
            data = data[0, :, :]
        filename = dbio.writeGeotif(lat, lon, res, data)
        dbio.ingest(dbname, filename, t, table, resample, overwrite)
        os.remove(filename)
    else:
        log.warning("No data were available to import into {0} for {1}.".format(table, t.strftime("%Y-%m-%d")))


def validate(dbname, table, dt):
    """Validates that data were correctly downloaded."""
    #edits by Vikalp Oct 2022
    if type(table) == dict:
        table = list(table.values())
    else:
        table = [table]
    for i in range(0,len(table)):
        tb = table[i]
        invalid = []
        db = dbio.connect(dbname)
        cur = db.cursor()
        ndays = (dt[1] - dt[0]).days + 1
        sql = "select count(*) from {0} where fdate>=date'{1}' and fdate<=date'{2}'".format(tb, dt[0].strftime("%Y-%m-%d"), dt[1].strftime("%Y-%m-%d"))
        cur.execute(sql)
        results = cur.fetchone()
        if results[0] < ndays:
            for d in range(ndays):
                t = dt[0] + timedelta(d)
                sql = "select fdate from {0} where fdate=date'{1}'".format(tb, t.strftime("%Y-%m-%d"))
                cur.execute(sql)
                results = cur.fetchone()
                if results is None:
                    invalid.append(t)
        return invalid
