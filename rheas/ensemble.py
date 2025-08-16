""" Class definition for the ensemble interface

.. module:: ensemble
   :synopsis: Definition of the ensemble class

.. moduleauthor:: Kostas Andreadis <kandread@umass.edu>

"""

import logging
import os
import random
import shutil
import sys
import tempfile
import xarray as xr
import rioxarray as rx
import pandas as pd
import threading
from datetime import date, timedelta, datetime
from multiprocessing import Process
from importlib import import_module


from . import dbio, vic


class Ensemble:

    def __init__(self, nens, dbname, resolution, startyear, startmonth, startday,
                 endyear, endmonth, endday, name=""):
        """Create an ensemble of models with size *nens*."""
        self.nens = nens
        self.models = []
        self.name = name
        self.statefiles = []
        self.res = resolution
        self.startyear, self.startmonth, self.startday = startyear, startmonth, startday
        self.endyear, self.endmonth, self.endday = endyear, endmonth, endday
        self.dbname = dbname
        for _ in range(nens):
            modelpath = tempfile.mkdtemp(dir=".")
            model = vic.VIC(modelpath, dbname, resolution, startyear, startmonth, startday,
                            endyear, endmonth, endday, name=name)
            self.models.append(model)

    def __getitem__(self, m):
        """Return a model instance."""
        return self.models[m]

    def __len__(self):
        """Return ensemble size."""
        return len(self.models)

    def __iter__(self):
        """Return an iterator to model ensemble members."""
        return iter(self.models)

    def writeParamFiles(self, savestate=False, statefile=""):
        """Write model parameter file for each ensemble member."""
        for e, model in enumerate(self.models):
            if len(self.statefiles) > 0:
                model.writeParamFile(state_file=self.statefiles[e], save_state=savestate)
            else:
                model.writeParamFile(state_file=statefile, save_state=savestate)

    def writeSoilFiles(self, shapefile):
        """Write soil parameter files based on domain shapefile."""
        self.models[0].writeSoilFile(shapefile)
        for model in self.models[1:]:
            shutil.copy(
                "{0}/soil.txt".format(self.models[0].model_path), "{0}/".format(model.model_path))
            model.lat = self.models[0].lat
            model.lon = self.models[0].lon
            model.gid = self.models[0].gid
            model.lgid = self.models[0].lgid
            model.depths = self.models[0].depths
            model.elev = self.models[0].elev

    def writeForcings(self, method, options):
        """Write forcings for the ensemble based on method (ESP, BCSD)."""
        log = logging.getLogger(__name__)
        if method.lower() == "esp":
            self._ESP(options)
        elif method.lower() == "nmme":
            self.__fromDataset("nmme", options)
        else:
            log.error("No appropriate method for generating meteorological forecast ensemble, exiting!")
            sys.exit()

    def __fromDataset(self, dataset, options):
        """Generate and write forcings by using a dataset-specific function."""
        dsmod = import_module("rheas.datasets.{0}".format(dataset))
        dsmod.generate(options, self)

    def _ESP(self, options):
        """Generate meteorological forcings using the Ensemble Streamflow Prediction method."""
        log = logging.getLogger(__name__)
        ndays = (date(self.endyear, self.endmonth, self.endday) -
                 date(self.startyear, self.startmonth, self.startday)).days
        db = dbio.connect(self.models[0].dbname)
        cur = db.cursor()
        if self.startmonth < self.endmonth:
            sql = "select distinct (date_part('year', fdate)) as year from precip.{0} where date_part('month', fdate) >= {1} and date_part('month', fdate) <= {2}".format(options['vic']['precip'], self.startmonth, self.endmonth)
        else:
            sql = "select distinct (date_part('year', fdate)) as year from precip.{0} where date_part('month', fdate) >= {1} or date_part('month', fdate) <= {2}".format(options['vic']['precip'], self.startmonth, self.endmonth)
        cur.execute(sql)
        years = map(lambda y: int(y[0]), cur.fetchall())
        ## Check if resample years have adequate record lengths
        for year in years:
            t = date(year, self.startmonth, self.startday) + timedelta(ndays)
            sql = "select * from precip.{0} where fdate=date'{1}'".format(options['vic']['precip'], t.strftime("%Y-%m-%d"))
            cur.execute(sql)
            if not bool(cur.rowcount):
                years.remove(year)
        if len(years) < 1:
            log.error("Not enough years in climatology to produce ESP forecast for selected dates, exiting!")
            sys.exit()
        random.shuffle(years)
        while len(years) < self.nens:
            years += years
        for e in range(self.nens):
            model = self.models[e]
            model.startyear = years[e]
            t = date(model.startyear, model.startmonth,
                     model.startday) + timedelta(ndays)
            model.endyear, model.endmonth, model.endday = t.year, t.month, t.day
            prec, tmax, tmin, wind = model.getForcings(options['vic'])
            model.writeForcings(prec, tmax, tmin, wind)
        cur.close()
        db.close()

    def setDates(self, startyear, startmonth, startday, endyear, endmonth, endday):
        """Set simulation dates for entire ensemble."""
        self.startyear, self.startmonth, self.startday = startyear, startmonth, startday
        self.endyear, self.endmonth, self.endday = endyear, endmonth, endday
        for m in self.models:
            m.startyear = startyear
            m.startmonth = startmonth
            m.startday = startday
            m.endyear = endyear
            m.endmonth = endmonth
            m.endday = endday

    def _ensembleTable(self, write, e):
        def write_wrapper(data, dates, tablename, initialize, skipsave):
            return write(data, dates, tablename, initialize, e, skipsave=skipsave)
        return write_wrapper

    def save(self, saveto, args, initialize=True):
        """Reads and saves selected output data variables from the ensemble into the database
        or a user-defined directory."""
        def ensembleTable(write, e):
            def write_wrapper(data, dates, tablename, initialize, skipsave):
                return write(data, dates, tablename, initialize, e)
            return write_wrapper
        for e in range(self.nens):
            model = self.models[e]
            if getattr(model.writeToDB, '__name__') != 'write_wrapper':
                # decorate function to add ensemble information
                model.writeToDB = ensembleTable(model.writeToDB, e + 1)
            if saveto == "db":
                if e > 0:
                    initialize = False
                model.save(saveto, args, initialize, ensemble=e)
            else:
                if e < 1:
                    if os.path.isdir(saveto):
                        shutil.rmtree(saveto)
                    elif os.path.isfile(saveto):
                        os.remove(saveto)
                    os.makedirs(saveto)
                model.save(saveto + "/{0}".format(e + 1), args, False)

    def run(self, vicexe):
        """Run ensemble of VIC models using multi-threading."""
        procs = [Process(target=self.models[e].run, args=(vicexe,))
                 for e in range(self.nens)]
        for p in procs:
            p.start()
        for p in procs:
            p.join()

    def statistics(self, name, stats=['mean', 'std']):
        """Calculate ensemble statistics."""
        log = logging.getLogger(__name__)
        db = dbio.connect(self.dbname)
        cur = db.cursor()
        if dbio.columnExists(self.dbname, self.name, name, "layer"):
            cur.execute("select max(layer) from {0}.{1}".format(self.name, name))
            nlayers = cur.fetchone()[0]
            lquery = "and layer={0}"
        else:
            nlayers = 1
            lquery = ""
        startdate = datetime(self.startyear, self.startmonth, self.startday)
        enddate = datetime(self.endyear, self.endmonth, self.endday)
        if dbio.columnExists(self.dbname, self.name, name, "ensemble"):
            cur.execute("select max(ensemble) from {0}.{1}".format(self.name, name))
            nens = cur.fetchone()[0]
            for lyr in range(1, nlayers+1):
                for e in range(nens):
                    xe = []
                    sql = "select fdate,filename from {0}.{1} where fdate>=date'{2}' and fdate<=date'{3}' and ensemble={4} {5} order by fdate".format(self.name, name, startdate.strftime("%Y-%m-%d"), enddate.strftime("%Y-%m-%d"), e+1, lquery.format(lyr))
                    cur.execute(sql)
                    results = cur.fetchall()
                    files = ["{0}/{1}".format(self[0].dbpath, r[1]) for r in results]
                    dates = [r[0] for r in results]
                    x = []
                    for f in files:
                        x.append(rx.open_rasterio(f, chunks=True, lock=threading.Lock()))
                    x = xr.concat(x, dim='time')
                    x = x.assign_coords({'time': pd.to_datetime(dates)})
                    x = x.squeeze().reset_coords('band', drop=True)
                    x = x.chunk(dict(time=len(files), x='auto', y='auto'))
                    xe.append(x)
                    xe = xr.concat(xe, dim='ensemble')
                    for stat in stats:
                        xs = eval("xe.{0}('ensemble').load()".format(stat))
                        for t in pd.date_range(startdate, enddate):
                            filename = "{0}/{1}_{2}_{3:02d}_{4}.tif".format(self[0].dbpath, name, t.strftime("%Y%m%d"), lyr, stat)
                            self[e]._writeRaster(xs.sel(time=t).data, filename)
        cur.close()
        db.close()

