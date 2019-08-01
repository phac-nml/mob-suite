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

from mob_suite.utils import default_database_dir, init_console_logger

config_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config.json')

with open(config_path, 'r') as configfile:
    config = json.load(configfile)

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--database_directory',
                        default=default_database_dir,
                        help='Directory to download databases to. Defaults to {}'.format(
                            default_database_dir))

    parser.add_argument('-v', '--verbose',
                        default=0,
                        action='count',
                        help='Set the verbosity level. Can by used multiple times')

    args = parser.parse_args()

    return args


def check_hash(filepath, hashsum):
    """
    Calculate the SHA256 of the given file and compare it to the known good
    value. Returns True if the hashes match, and False if the file is missing or
    the hashes do not match.

    :param filepath: Path to the file
    :param hashsum:  Expected SHA256 hash
    :return: True if the hashes match, False otherwise
    """

    try:
        with open(filepath, 'rb') as f:

            sha = hashlib.sha256()
            contents = f.read()
            sha.update(contents)

    except FileNotFoundError:

        return False

    h = sha.hexdigest()
    return h == hashsum


def download_to_file(url, file):
    """
    Downloads a file from a URL using curl.

    :param url:  Source URL from which the file is downloaded
    :param file: The destination file
    :return: None
    """

    with open(file, 'wb') as f:
        c = pycurl.Curl()
        # Redirects to https://www.python.org/.
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

    if fname.endswith(".zip"):

        with zipfile.ZipFile(fname, 'r') as zip_ref:
            zip_ref.extractall(outdir)

    elif fname.endswith(".gz"):

        outfile = os.path.join(outdir, fname.replace('.gz',''))

        with gzip.open(fname, 'rb') as f_in, open(outfile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    os.remove(fname)


def main():


    args = arguments()

    global logger
    logger = init_console_logger(args.verbose)

    database_directory = os.path.abspath(args.database_directory)


    # Helper function to simplify adding database_directory to everything
    prepend_db_dir = functools.partial(os.path.join, database_directory)

    logger.info('Initializing databases...this will take some time')

    # Find available threads and use the maximum number available for mash sketch but cap it at 32
    num_threads = min(multiprocessing.cpu_count(), 32)


    if not os.path.exists(database_directory):
        os.makedirs(database_directory)

    zip_file = prepend_db_dir('data.zip')
    plasmid_database_fasta_file = prepend_db_dir('ncbi_plasmid_full_seqs.fas')
    repetitive_fasta_file = prepend_db_dir('repetitive.dna.fas')
    mash_db_file =  prepend_db_dir('ncbi_plasmid_full_seqs.fas.msh')

    logger.info('Downloading databases...this will take some time')

    for db_mirror in config['db_mirrors']:

        logger.info('Trying mirror {}'.format(db_mirror))
        download_to_file(db_mirror, zip_file)


        if check_hash(zip_file, config['db_hash']):
            break   #do not try other mirror

    else:  # no break
        logger.error('Downloading databases failed, please check your internet connection and retry')
        sys.exit(-1)


    logger.info('Downloading databases successful, now building databases')

    extract(zip_file, database_directory)

    files = [prepend_db_dir(f)
             for f in os.listdir(database_directory)
             if f.endswith('.gz')]

    for file in files:

        extract(file, database_directory)

    #Initialize blast and mash databases
    logger.info('Building repetitive mask database')
    blast_runner = BlastRunner(repetitive_fasta_file, database_directory)
    blast_runner.makeblastdb(repetitive_fasta_file, 'nucl')

    logger.info('Building complete plasmid database')
    blast_runner = BlastRunner(plasmid_database_fasta_file, database_directory)
    blast_runner.makeblastdb(plasmid_database_fasta_file, 'nucl')

    logger.info('Sketching complete plasmid database')
    mObj = mash()
    mObj.mashsketch(plasmid_database_fasta_file,
                    mash_db_file,
                    num_threads=num_threads)

    status_file = prepend_db_dir('status.txt')

    with open(status_file, 'w') as f:
        download_date = datetime.datetime.today().strftime('%Y-%m-%d')

        f.write("Download date: {}".format(download_date))


# call main function
if __name__ == '__main__':
    main()
