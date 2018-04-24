#!/usr/bin/env python
from __future__ import print_function


import pycurl
import tarfile
import zipfile

def download_to_file(url,file):
    with open(file, 'wb') as f:
        c = pycurl.Curl()
        # Redirects to https://www.python.org/.
        c.setopt(c.URL, url)
        # Follow redirect.
        c.setopt(c.FOLLOWLOCATION, True)
        c.setopt(c.WRITEDATA, f)
        c.perform()
        c.close()

def extract(fname,outdir):
    if (fname.endswith("tar.gz")):
        tar = tarfile.open(fname, "r:gz")
        tar.extractall()
        tar.close()
    elif (fname.endswith("tar")):
        tar = tarfile.open(fname, "r:")
        tar.extractall(outdir)
        tar.close()
    elif(fname.endswith("zip")):
        zip_ref = zipfile.ZipFile(fname, 'r')
        zip_ref.extractall(outdir)
        zip_ref.close()


download_to_file('https://ndownloader.figshare.com/articles/5841882?private_link=a4c92dd84f17b2cefea6','/Users/jrobertson/Desktop/out.txt')
