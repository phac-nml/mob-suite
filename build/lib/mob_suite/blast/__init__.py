from datetime import datetime
import logging
import shutil

from subprocess import Popen, PIPE
import os

import pandas as pd
from pandas.io.common import EmptyDataError
import re


BLAST_TABLE_COLS = '''
qseqid
sseqid
qlen
slen
qstart
qend
sstart
send
length
mismatch
pident
qcovhsp
qcovs
sstrand
evalue
bitscore
'''.strip().split('\n')


class BlastRunner:

    def __init__(self, fasta_path, tmp_work_dir):
        self.fasta_path = fasta_path

    def makeblastdb(self,fasta_path,dbtype):
        p = Popen(['makeblastdb',
                  '-in', fasta_path,
                  '-dbtype',dbtype],
                  stdout=PIPE,
                  stderr=PIPE)
        p.wait()
        stdout = p.stdout.read()
        stderr = p.stderr.read()



    def run_tblastn(self, query_fasta_path, blast_task, db_path, db_type, min_cov, min_ident, evalue,blast_outfile,num_threads=1):

        p = Popen(['tblastn',
                   '-query', query_fasta_path,
                   '-num_threads','{}'.format(num_threads),
                   '-db', '{}'.format(db_path),
                   '-evalue', '{}'.format(evalue),
                   '-out', blast_outfile,
                   '-outfmt', '6 {}'.format(' '.join(BLAST_TABLE_COLS))],
                  stdout=PIPE,
                  stderr=PIPE)

        p.wait()
        stdout = p.stdout.read()
        stderr = p.stderr.read()
        if stdout is not None and stdout != '':
            logging.debug('blastn on db {} and query {} STDOUT: {}'.format(query_fasta_path, db_path, stdout))

        if stderr is not None and stderr != '':
            logging.debug('blastn on db {} and query {} STDERR: {}'.format(query_fasta_path, db_path, stderr))
            if os.path.exists(blast_outfile):
                return blast_outfile
            else:
                ex_msg = 'tblastn on db {} and query {} did not produce expected output file at {}'.format(
                    query_fasta_path,
                    db_path,
                    blast_outfile)
                logging.error(ex_msg)
                raise Exception(ex_msg)

    def run_blast(self, query_fasta_path, blast_task, db_path, db_type, min_cov, min_ident, evalue,blast_outfile,num_threads=1,word_size=11):

        p = Popen(['blastn',
                   '-task', blast_task,
                   '-query', query_fasta_path,
                   '-db', '{}'.format(db_path),
                   '-num_threads','{}'.format(num_threads),
                   '-evalue', '{}'.format(evalue),
                   '-dust', 'yes',
                   '-perc_identity', '{}'.format(min_ident),
                   '-out', blast_outfile,
                   '-outfmt', '6 {}'.format(' '.join(BLAST_TABLE_COLS))],
                  stdout=PIPE,
                  stderr=PIPE)

        p.wait()
        stdout = p.stdout.read()
        stderr = p.stderr.read()
        if stdout is not None and stdout != '':
            logging.debug('blastn on db {} and query {} STDOUT: {}'.format(query_fasta_path, db_path, stdout))

        if stderr is not None and stderr != '':
            logging.debug('blastn on db {} and query {} STDERR: {}'.format(query_fasta_path, db_path, stderr))
            if os.path.exists(blast_outfile):
                return blast_outfile
            else:
                ex_msg = 'blastn on db {} and query {} did not produce expected output file at {}'.format(
                    query_fasta_path,
                    db_path,
                    blast_outfile)
                logging.error(ex_msg)
                raise Exception(ex_msg)


class BlastReader:
    df = None


    def __init__(self, blast_outfile):
        """Read BLASTN output file into a pandas DataFrame
        Sort the DataFrame by BLAST bitscore.
        If there are no BLASTN results, then no results can be returned.
        Args:
            blast_outfile (str): `blastn` output file path
        Raises:
            EmptyDataError: No data could be parsed from the `blastn` output file
        """
        self.blast_outfile = blast_outfile
        try:

            self.df = pd.read_table(self.blast_outfile, header=None)
            self.df.columns = BLAST_TABLE_COLS

            logging.debug(self.df.head())
            self.is_missing = False

        except EmptyDataError as exc:
            logging.warning('No BLASTN results to parse from file %s', blast_outfile)
            self.is_missing = True
            self.df = pd.DataFrame(index=['A'], columns='A')


    def df_dict(self):
        if not self.is_missing:
            return self.df.to_dict()

