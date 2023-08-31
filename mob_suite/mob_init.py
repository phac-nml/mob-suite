#!/usr/bin/env python3
from mob_suite.version import __version__
import os, pycurl, tarfile, zipfile, gzip, multiprocessing, sys
import argparse
import hashlib
import json
import functools
from mob_suite.blast import BlastRunner
from mob_suite.wrappers import mash
import shutil
import datetime
import time #waiting for other processes
from mob_suite.utils import init_console_logger
from mob_suite.constants import  default_database_dir
from ete3 import NCBITaxa
from pathlib import Path

logger = init_console_logger(3)
config_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config.json')

with open(config_path, 'r') as configfile:
    config = json.load(configfile)

def arguments():


    parser = argparse.ArgumentParser(
        description="MOB-INIT: initialize databases version: {}".format(
            __version__))
    parser.add_argument('-d', '--database_directory',
                        default=default_database_dir,
                        help='Directory to download databases to. Defaults to {}'.format(
                            default_database_dir))

    parser.add_argument('-v', '--verbose',
                        default=0,
                        action='count',
                        help='Set the verbosity level. Can by used multiple times')

    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__ + "")
    args, unknown_args = parser.parse_known_args()
    return args


def check_hash_or_size(filepath, hashsums):
    """
    Calculate the SHA256 of the given file and compare it to the known good
    value. Returns True if the hashes match, and False if the file is missing or
    the hashes do not match.
    Otherwise uses minimum file size of >= 215 MB
    FigShare mirror returns non hash identical files for some reason

    :param filepath: Path to the file
    :param hashsums:  Expected SHA256 hashes (given different database versions)
    :return: True if the hashes match, False otherwise
    """

    try:
        with open(filepath, 'rb') as f:
            sha = hashlib.sha256()
            contents = f.read()
            sha.update(contents)

    except FileNotFoundError:
        return False

    h_obs = sha.hexdigest()
    logger.info("Download data.zip sha256 checksum is {}".format(h_obs))
    logger.info("Download data.zip size in bytes is {}".format(os.path.getsize(filepath)))
    hash_check_bool = any([refhash == h_obs for refhash in hashsums])
    size_check_bool = (os.path.getsize(filepath) >= 225377917)
    return any([hash_check_bool, size_check_bool])


def download_to_file(url, file):
    """
    Downloads a file from a URL using curl.

    :param url:  Source URL from which the file is downloaded
    :param file: The destination file
    :return: None
    """

    with open(file, 'wb') as f:
        c = pycurl.Curl()
        c.setopt(c.URL, url)
        # Follow redirect.
        c.setopt(c.FOLLOWLOCATION, True)
        c.setopt(c.WRITEDATA, f)
        c.setopt(c.NOPROGRESS, False)
        c.perform()
        c.close()


def extract(fname, outdir):
    """
    Decompress a zip or gzip archive. Decompression method is selected based
    on file extension. Following extraction, the original archive is deleted.

    :param fname:  Path to the archive to be extracted
    :param outdir: Directory into which the results are placed
    :return: None
    """

    logger.info(f'Decompressing {fname}')
    shutil.unpack_archive(fname, outdir)
    dir_name = os.path.join(outdir,os.path.basename(fname))
    for ext in ['.tar.gz','.zip','.gz']:
        dir_name = dir_name.replace(ext,'')

    src_files = os.listdir(dir_name)
    for file_name in src_files:
        full_file_name = os.path.join(dir_name, file_name)
        if os.path.isfile(full_file_name):
            shutil.copyfile(full_file_name, os.path.join(outdir,file_name))
    shutil.rmtree(dir_name)
    os.remove(fname)

def main():
    
    args = arguments()
    

    database_directory = os.path.abspath(args.database_directory)


    if os.path.exists(database_directory) == False:
        os.makedirs(database_directory)
        logger.info("Database directory folder created at {}".format(database_directory))
    else:
        logger.info("Database directory folder already exists at {}".format(database_directory))

    
    # Helper function to simplify adding database_directory to everything
    prepend_db_dir = functools.partial(os.path.join, database_directory)

    lockfilepath=os.path.join(database_directory,".lock")
    status_file = prepend_db_dir('status.txt')

    if os.path.exists(lockfilepath) == False:
        try:
            open(file=lockfilepath, mode="w").close()
            logger.info("Placed lock file at {}".format(lockfilepath))
        except Exception as e:
            logger.error("Failed to place a lock file at {}. Database diretory can not be accessed. Wrong path?".format(lockfilepath))
            logger.error("{}".format(e))
            pass
    else:
        while os.path.exists(lockfilepath):
            elapsed_time = time.time() - os.path.getmtime(lockfilepath)
            logger.info("Lock file found at {}. Waiting for other processes to finish database init ...".format(lockfilepath))
            logger.info("Elapsed time {} min. Will continue processing after 16 min mark.".format(int(elapsed_time/60)))
            if elapsed_time >= 1000:
                logger.info("Elapsed time {} min. Assuming previous process completed all init steps. Continue ...".format(int(elapsed_time/60)))
                try: #if previous process failed, no processes are running and > 16 min passed since the lock was created
                    os.remove(lockfilepath)
                except: #continue if file was removed by other process
                    pass
                break
            time.sleep(60) #recheck every 1 min if lock file was removed
        logger.info("Lock file no longer exists. Assuming init process completed successfully")
        return 0



    logger.info('Initializing databases...this will take some time')
    # Find available threads and use the maximum number available for mash sketch but cap it at 32
    try:
        num_threads = len(os.sched_getaffinity(0))
    except AttributeError:
        num_threads = multiprocessing.cpu_count()

    if num_threads > 32:
        num_threads = 32
    if num_threads < 1:
        num_threads = 1

    if not os.path.exists(database_directory):
        os.makedirs(database_directory)

    zip_file = prepend_db_dir('data.tar.gz')
    plasmid_database_fasta_file = prepend_db_dir('ncbi_plasmid_full_seqs.fas')
    repetitive_fasta_file = prepend_db_dir('repetitive.dna.fas')
    mash_db_file =  prepend_db_dir('ncbi_plasmid_full_seqs.fas.msh')

    logger.info('Downloading databases...this will take some time')

    for db_mirror in config['db_mirrors']:
        try:
            logger.info('Trying mirror {}'.format(db_mirror))
            download_to_file(db_mirror, zip_file)
            break
        except Exception as e:
            logger.error("Download failed with error {}. Removing lock file".format(str(e)))
            os.remove(lockfilepath)
            sys.exit(-1)


    logger.info("Downloading databases successful, now building databases at {}".format(database_directory))
    extract(zip_file, database_directory)

    files = [prepend_db_dir(f)
             for f in os.listdir(database_directory)
             if f.endswith('.gz')]

    for file in files:

        extract(file, database_directory)

    #Initialize blast and mash databases
    try:
        logger.info('Building repetitive mask database')
        blast_runner = BlastRunner(repetitive_fasta_file, database_directory)
        blast_runner.makeblastdb(repetitive_fasta_file, 'nucl',logger)

        logger.info('Building complete plasmid database')
        blast_runner = BlastRunner(plasmid_database_fasta_file, database_directory)
        blast_runner.makeblastdb(plasmid_database_fasta_file, 'nucl',logger,True)

        logger.info('Sketching complete plasmid database')
        mObj = mash()
        mObj.mashsketch(plasmid_database_fasta_file,
                        mash_db_file,
                        num_threads=num_threads)
    except Exception as e:
        logger.error('Downloading databases failed, please check your internet connection and retry')
        logger.error("Process failed with error {}. Removing lock file".format(e))
        os.remove(lockfilepath)
        sys.exit(-1)

    try:
        logger.info("Init ete3 library ...")
        ete3taxadbpath = os.path.abspath(os.path.join(database_directory,"taxa.sqlite"))
        NCBITaxa(dbfile=ete3taxadbpath) #the creatuib if NCBITaxa class triggers update_taxonomy_database()
    except Exception as e:
        logger.error("Init of ete3 library failed with error {}. Removing lock file".format(e))
        os.remove(lockfilepath)
        sys.exit(-1)

    try:
        os.remove(os.path.join(os.getcwd(), "taxdump.tar.gz"))
        logger.info("Removed residual taxdump.tar.gz as ete3 is not doing proper cleaning job.")
    except:
        pass

    with open(status_file, 'w') as f:
        download_date = datetime.datetime.today().strftime('%Y-%m-%d')
        f.write("Download date: {}. Removing lock file.".format(download_date))
        try:
            os.remove(lockfilepath)
        except:
            logger.warning("Lock file is already removed by some other process.")
            pass


    logger.info("MOB init completed successfully")
    return 0

# call main function
if __name__ == '__main__':
    main()
