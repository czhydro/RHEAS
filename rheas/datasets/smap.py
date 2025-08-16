"""Class definition for the SMAP Soil Mositure data type.

.. module:: smap
   :synopsis: Definition of the SMAP class

.. moduleauthor:: Kostas Andreadis <kandread@umass.edu>

"""

import logging
import subprocess
import tempfile
from datetime import timedelta

import numpy as np

from .soilmoist import Soilmoist

from .. import dbio
from . import datasets, earthdata

table = "soilmoist.smap"


def dates(dbname):
    dts = datasets.dates(dbname, table)
    return dts


def download(dbname, dts, bbox=None, enhanced=False):
    """Downloads SMAP soil mositure data for a set of dates *dt*
    and imports them into the PostGIS database *dbname*. Optionally
    uses a bounding box to limit the region with [minlon, minlat, maxlon, maxlat]."""
    log = logging.getLogger(__name__)
    if enhanced:
        res = 0.09
        url = "https://n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP_E.005"
    else:
        res = 0.36
        url = "https://n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP.008"
    for dt in [dts[0] + timedelta(tt) for tt in range((dts[-1] - dts[0]).days + 1)]:
        try:
            outpath, fname = earthdata.download("{0}/{1}".format(url, dt.strftime("%Y.%m.%d")), "SMAP_L3_SM_P_\S*.h5")
            h5file = "{0}/{1}".format(outpath, fname)
            varname = None
            for vgroup in ["Soil_Moisture_Retrieval_Data_AM", "Soil_Moisture_Retrieval_Data_PM"]:
                ext = "" if vgroup.endswith("AM") else "_pm"
                varname = "soil_moisture{0}".format(ext)
                ftmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
                filename = ftmp.name
                ftmp.close()
                cmd = "gdal_translate -a_srs epsg:6933 -a_ullr -17367530.45 7314540.11 17367530.45 -7314540.11 HDF5:'{0}'://{1}/{2} {3}".format(h5file, vgroup, varname, filename)
                subprocess.call(cmd, shell=True)
                proc = subprocess.Popen(["gdalwarp", "-t_srs", "epsg:4326", filename, filename.replace(".tif", "_ll.tif")], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                out, err = proc.communicate()
                log.debug(out)
                if bbox is None:
                    pstr = []
                else:
                    pstr = ["-projwin", str(bbox[0]), str(bbox[3]), str(bbox[2]), str(bbox[1])]
                proc = subprocess.Popen(["gdal_translate"] + pstr + ["-a_srs", "epsg:4326", filename.replace(".tif", "_ll.tif"), filename.replace(".tif", "_s.tif")], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                out, err = proc.communicate()
                log.debug(out)
                sname, tname = table.split(".")
                overpass = 2 if ext.endswith("pm") else 1
                db = dbio.connect(dbname)
                cur = db.cursor()
                # check if same date and overpass exists and delete it
                if dbio.tableExists(dbname, sname, tname):
                    cur.execute("delete from {0} where fdate=date'{1}' and overpass={2}".format(table, dt.strftime("%Y-%m-%d"), overpass))
                    db.commit()
                dbio.ingest(dbname, filename.replace(".tif", "_s.tif"), dt, table, False, False)
                if not dbio.columnExists(dbname, sname, tname, "overpass"):
                    cur.execute("alter table {0} add overpass integer".format(table))
                    db.commit()
                    # 1 is descending (AM), 2 is ascending (PM) overpass
                cur.execute("update {0} set overpass={1} where overpass is null".format(table, overpass))
                db.commit()
                cur.close()
                db.close()
        except:
            log.warning("No SMAP data available for {0}.".format(dt.strftime("%Y-%m-%d")))


class Smap(Soilmoist):

    def __init__(self, uncert=None):
        """Initialize SMAP soil moisture object."""
        super(Smap, self).__init__(uncert)
        self.res = 0.36
        self.stddev = 0.03
        self.tablename = "soilmoist.smap"
