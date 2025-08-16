""" Definition for RHEAS Earthdata module.

.. module:: earthdata
   :synopsis: Definition of the Earthdata module

.. moduleauthor:: Kostas Andreadis <kandread@umass.edu>

"""

import re
import os
import tempfile
import logging

import requests
from lxml import html


def download(url, filepattern):
    """Download data files from Earthdata search."""
    log = logging.getLogger(__name__)
    try:
        username = os.environ['EARTHDATA_USER']
        password = os.environ['EARTHDATA_PASS']
    except KeyError:
        username = password = ""
        log.error("Earthdata username and/or password need to be set as environment variables!")
    with requests.Session() as session:
        r1 = session.request('get', url)
        r = session.get(r1.url, auth=(username, password))
        if r.ok:
            links = html.fromstring(r.content).xpath('//a/@href')
            matches = [re.match(filepattern, link) for link in links]
            filenames = [m.group(0) for m in matches if m]
            # make sure list has unique filenames
            filenames = list(set(filenames))
            tmppath = tempfile.mkdtemp()
            for filename in filenames:
                rfile = session.get("{0}/{1}".format(url, filename), auth=(username, password))
                with open("{0}/{1}".format(tmppath, filename), 'wb') as fout:
                    for chunk in rfile:
                        fout.write(chunk)
                        filename = filenames[0]
        else:
            tmppath, filename = None
    return tmppath, filename
