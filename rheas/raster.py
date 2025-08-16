""" RHEAS module for raster manipulation within the PostGIS database

.. module:: raster
   :synopsis: Manipulate PostGIS rasters

.. moduleauthor:: Kostas Andreadis <kandread@umass.edu>

"""

import logging

from . import dbio


class TileReader:
    """Helper class to retrieve raster tile from database."""

    def __init__(self, dbname, rtable, startyear, startmonth, startday, endyear, endmonth, endday):
        self.dbname = dbname
        self.rtable = rtable
        self.startyear = startyear
        self.startmonth = startmonth
        self.startday = startday
        self.endyear = endyear
        self.endmonth = endmonth
        self.endday = endday

    def __call__(self, t):
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        var = self.rtable.split(".")[0]
        sql = "select gid,fdate,st_nearestvalue(rast,x,y) from {0},{1}_xy where rid=tile and tile={8} and fdate>=date'{2}-{3}-{4}' and fdate<=date'{5}-{6}-{7}' order by gid,fdate".format(
            self.rtable, var, self.startyear, self.startmonth, self.startday, self.endyear, self.endmonth, self.endday, t)
        cur.execute(sql)
        data = cur.fetchall()
        return data
