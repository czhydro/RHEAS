""" RHEAS module for retrieving maximum and minimum
temperature and wind speed from the NCEP Reanalysis.

.. module:: ncep
   :synopsis: Retrieve NCEP meteorological data

.. moduleauthor:: Kostas Andreadis <kandread@umass.edu>

"""

from datetime import timedelta, datetime

import numpy as np

from . import datasets
from .decorators import netcdf

table = ["tmax.ncep", "tmin.ncep", "wind.ncep"]


def dates(dbname):
    dts = datasets.dates(dbname, "wind.ncep")
    return dts


@netcdf
def fetch_tmax(dbname, dt, bbox):
    """Downloads maximum temperature from NCEP Reanalysis."""
    url = "http://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis/Dailies/surface_gauss/tmax.2m.gauss.{0}.nc".format(dt[0].year)
    varname = "tmax"
    return url, varname, bbox, dt


@netcdf
def fetch_tmin(dbname, dt, bbox):
    """Downloads minimum temperature from NCEP Reanalysis."""
    url = "https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis/Dailies/surface_gauss/tmin.2m.gauss.{0}.nc".format(dt[0].year)
    varname = "tmin"
    return url, varname, bbox, dt


@netcdf
def fetch_uwnd(dbname, dt, bbox):
    """Downloads U-component wind speed from NCEP Reanalysis."""
    url = "https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis/Dailies/surface_gauss/uwnd.10m.gauss.{0}.nc".format(dt[0].year)
    varname = "uwnd"
    return url, varname, bbox, dt


@netcdf
def fetch_vwnd(dbname, dt, bbox):
    """Downloads U-component wind speed from NCEP Reanalysis."""
    url = "https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis/Dailies/surface_gauss/vwnd.10m.gauss.{0}.nc".format(dt[0].year)
    varname = "vwnd"
    return url, varname, bbox, dt


def download(dbname, dts, bbox=None):
    """Downloads NCEP Reanalysis data."""
    res = 1.875
    for yr in range(dts[0].year, dts[-1].year + 1):
        t0 = max(datetime(yr, 1, 1), dts[0])
        t1 = min(datetime(yr, 12, 31), dts[-1])
        dt = [t0, t1]
        tmax, lat, lon, _ = fetch_tmax(dbname, dt, bbox)
        tmin, _, _, _ = fetch_tmin(dbname, dt, bbox)
        uwnd, _, _, _ = fetch_uwnd(dbname, dt, bbox)
        vwnd, _, _, _ = fetch_vwnd(dbname, dt, bbox)
        wnd = np.sqrt(uwnd**2 + vwnd**2)
        tmax -= 273.15
        tmin -= 273.15
        for ti, t in enumerate([t0 + timedelta(tt) for tt in range((t1 - t0).days + 1)]):
            datasets.ingest(dbname, "tmax.ncep", tmax[ti, :, :], lat, lon, res, t)
            datasets.ingest(dbname, "tmin.ncep", tmin[ti, :, :], lat, lon, res, t)
            datasets.ingest(dbname, "wind.ncep", wnd[ti, :, :], lat, lon, res, t)
