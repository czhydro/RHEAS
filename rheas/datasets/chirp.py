""" RHEAS module for retrieving rainfall data from the Climate Hazard Group
    InfraRed Precipitation (CHIRP) data archive.

.. module:: chirp
   :synopsis: Retrieve CHIRP rainfall data

.. moduleauthor:: Kostas Andreadis <kandread@umass.edu>

"""

from datetime import timedelta

from . import datasets
from .decorators import geotiff, http

table = "precip.chirp"


@geotiff
@http
def fetch(dbname, dt, bbox):
    """Downloads CHIRPS rainfall data from the data server."""
    url = "https://data.chc.ucsb.edu/products/CHIRP/daily/{0:04d}/chirp.{0:04d}.{1:02d}.{2:02d}.tif.gz"
    return url, bbox, dt


def download(dbname, dts, bbox=None):
    res = 0.05
    for dt in [dts[0] + timedelta(tt) for tt in range((dts[-1] - dts[0]).days + 1)]:
        data, lat, lon, t = fetch(dbname, dt, bbox)
        datasets.ingest(dbname, table, data, lat, lon, res, t)


def dates(dbname):
    dts = datasets.dates(dbname, table)
    return dts
