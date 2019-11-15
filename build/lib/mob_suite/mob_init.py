#!/usr/bin/env python3
from mob_suite.version import __version__
import os, pycurl, tarfile, zipfile, gzip, multiprocessing, sys
from mob_suite.blast import BlastRunner
from mob_suite.wrappers import mash
from os import listdir
from os.path import isfile, join
import shutil
import datetime
import logging

LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'

def init_console_logger(lvl):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    report_lvl = logging_levels[lvl]

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)
    return logging


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
    elif(fname.endswith("gz")):
        outfile = os.path.join(outdir,fname.replace('.gz',''))
        with gzip.open(fname, 'rb') as f_in:
            with open(outfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            f_in.close()
            f_out.close()
            os.remove(fname)

def main():
    logging = init_console_logger(2)
    logging.info('Initilizating databases...this will take some time')

    #Find available threads and use the maximum number available for mash sketch but cap it at 32
    num_threads = multiprocessing.cpu_count()
    if num_threads > 32:
        num_threads = 32

    database_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)),'databases/')
    zip_file = os.path.join(database_directory,'data.zip')
    plasmid_database_fasta_file = os.path.join(database_directory,'ncbi_plasmid_full_seqs.fas')
    repetitive_fasta_file = os.path.join(database_directory,'repetitive.dna.fas')
    mash_db_file =  os.path.join(database_directory,'ncbi_plasmid_full_seqs.fas.msh')
    logging.info('Downloading databases...this will take some time')
    download_to_file('https://ndownloader.figshare.com/articles/5841882?private_link=a4c92dd84f17b2cefea6',zip_file)
    if (not os.path.isfile(zip_file)):
        logging.error('Downloading databases failed, please check your internet connection and retry')
        sys.exit(-1)
    else:
        logging.info('Downloading databases successful, now building databases')
    extract(zip_file,database_directory)
    os.remove(zip_file)
    files = [f for f in listdir(database_directory) if isfile(join(database_directory, f))]
    for file in files:

        if file.endswith('gz'):
            extract(os.path.join(database_directory,file), database_directory)

    #Initilize blast and mash daatabases
    logging.info('Building repetive mask database')
    blast_runner = BlastRunner(repetitive_fasta_file, database_directory)
    blast_runner.makeblastdb(repetitive_fasta_file, 'nucl')
    logging.info('Building complete plasmid database')
    blast_runner = BlastRunner(plasmid_database_fasta_file, database_directory)
    blast_runner.makeblastdb(plasmid_database_fasta_file, 'nucl')
    logging.info('Sketching complete plasmid database')
    mObj = mash()
    mObj.mashsketch(plasmid_database_fasta_file,mash_db_file,num_threads=num_threads)
    status_file = os.path.join(database_directory,'status.txt')
    with open(status_file, 'w') as f:
        f.write("Download date: {}".format(datetime.datetime.today().strftime('%Y-%m-%d')))
    f.close()

# call main function
if __name__ == '__main__':
    main()