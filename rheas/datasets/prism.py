""" RHEAS module for retrieving PRISM meteorological data.

.. module:: prism
   :synopsis: Retrieve PRISM meteorological data

.. moduleauthor:: Kostas Andreadis <kandread@umass.edu>

"""


import logging, os
import subprocess
import tempfile, glob
import zipfile
import numpy as np
from datetime import datetime
from ftplib import FTP


from .. import dbio
from . import datasets

table = {"ppt": "precip.prism", "tmax": "tmax.prism", "tmin": "tmin.prism"}


def dates(dbname):
    dts = datasets.dates(dbname, table['ppt'])#---------needs tmin, tmax?
    return dts


def _downloadVariable(varname, dbname, dts, bbox):
    """Downloads the PRISM data products for a specific variable and a set of
    dates *dt*. *varname* can be ppt, tmax or tmin."""
    log = logging.getLogger(__name__)
    url = "prism.oregonstate.edu"
    ftp = FTP(url)
    ftp.login()
    ftp.cwd(f"/time_series/us/an/4km/{varname}/daily/")
    print(ftp.dir())

    #ftp.cwd("daily/{0}".format(varname)) ---------original
    outpath = tempfile.mkdtemp()
    print(outpath)

    years = np.sort(list(set([t.year for t in dts])))
    print(years)
      
    for yr in range(years[0],years[1]+1):
        ftp.cwd("{0}".format(yr))
        full_flist = ftp.nlst()
        filenames = []
        for f in full_flist:
            tmp = f.split('_')[-1].split('.')[0]
            if datetime.strptime(tmp,"%Y%m%d") >= dts[0] and datetime.strptime(tmp,'%Y%m%d')<=dts[-1]:
                filenames.append(f)

        #filenames = [f for f in ftp.nlst() if datetime.strptime(f.split("_")[-1], "%Y%m%d") >= dts[0] and datetime.strptime(f.split("_")[-1], "%Y%m%d") <= dts[-1]]
        filename = np.sort(filenames)

        for fname in filename:
            print(fname)
            dt = datetime.strptime(fname.split("_")[-1].split('.')[0], "%Y%m%d")
            print(dt)

            with open("{0}/{1}".format(outpath, fname), 'wb') as f:
                ftp.retrbinary("RETR {0}".format(fname), f.write)
            if fname.endswith("zip"):
                fz = zipfile.ZipFile("{0}/{1}".format(outpath, fname))
                for ls in fz.namelist():
                    if ls[-4:] == '.tif':
                        lfilename = ls
                #lfilename = filter(lambda s: s.endswith("bil"), fz.namelist())[0]
                fz.extractall(outpath)
            #else:
            #    lfilename = fname[:-4]
            #tfilename = lfilename.replace(".bil", ".tif")
            tfilename = fname.replace('.zip','_bb.tif')
            print(tfilename)

            if bbox is not None:
                proc = subprocess.Popen(["gdal_translate", "-projwin", "{0}".format(bbox[0]), "{0}".format(bbox[3]), "{0}".format(bbox[2]), "{0}".format(bbox[1]), "{0}/{1}".format(outpath, lfilename), "{0}/{1}".format(outpath, tfilename)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                out, err = proc.communicate()
                log.debug(out)
                #dbio.writeGeotif(lat, lon, res, data, tfilename)
                dbio.ingest(dbname, "{0}/{1}".format(outpath, tfilename), dt, table[varname], True)
                files2rm = glob.glob(outpath+'/*')
                if len(files2rm)>0:
                  for file2rm in files2rm:
                    os.remove(file2rm)
            else:
                dbio.ingest(dbname, "{0}/{1}".format(outpath, lfilename), dt, table[varname], True)
                os.remove(outpath+'/'+lfilename)
        ftp.cwd("..")


def download(dbname, dts, bbox):
    """Downloads the PRISM data products for a set of
    dates *dt* and imports them into the PostGIS database *dbname*."""
    for varname in ["ppt", "tmax", "tmin"]:
        _downloadVariable(varname, dbname, dts, bbox)
