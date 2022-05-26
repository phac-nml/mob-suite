from datetime import datetime
import logging, sys
from pathlib import Path

from subprocess import Popen, PIPE
import os

import pandas as pd
from pandas.errors import EmptyDataError



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

#inherit logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.getLogger().getEffectiveLevel())



class BlastRunner:

    def __init__(self, fasta_path, tmp_work_dir):
        self.fasta_path = fasta_path

    def makeblastdb(self,fasta_path,dbtype,logging,parse_seqids=False):
        if parse_seqids:
            p = Popen(['makeblastdb',
                      '-in', Path(fasta_path),
                       '-parse_seqids',
                      '-dbtype',dbtype],
                      stdout=PIPE,
                      stderr=PIPE)
        else:
            p = Popen(['makeblastdb',
                      '-in', Path(fasta_path),
                      '-dbtype',dbtype],
                      stdout=PIPE,
                      stderr=PIPE)

        p.wait()
        stdout = str(p.stdout.read())
        stderr = str(p.stderr.read())
        if stderr is not None and stderr != '' and stderr != "b''" :
            logging.error('makeblastdb on {} had the following messages STDERR: {}'.format(fasta_path, stderr))
            return False
        return True




    def run_tblastn(self, query_fasta_path, blast_task, db_path, db_type, min_cov, min_ident, evalue,blast_outfile,logging,num_threads=1,max_target_seqs=100000000,seq_id_file=None):
        if seq_id_file:
            p = Popen(['tblastn',
                       '-query', Path(query_fasta_path),
                       '-seqidlist', '{}'.format(seq_id_file),
                       '-num_threads','{}'.format(num_threads),
                       '-db', '{}'.format(db_path),
                       '-evalue', '{}'.format(evalue),
                       '-out', blast_outfile,
                       '-max_target_seqs','{}'.format(max_target_seqs),
                       '-outfmt', '6 {}'.format(' '.join(BLAST_TABLE_COLS))],
                      stdout=PIPE,
                      stderr=PIPE)
        else:
            p = Popen(['tblastn',
                       '-query', Path(query_fasta_path),
                       '-num_threads','{}'.format(num_threads),
                       '-db', '{}'.format(Path(db_path)),
                       '-evalue', '{}'.format(evalue),
                       '-out', blast_outfile,
                       '-max_target_seqs','{}'.format(max_target_seqs),
                       '-outfmt', '6 {}'.format(' '.join(BLAST_TABLE_COLS))],
                      stdout=PIPE,
                      stderr=PIPE)
        p.wait()
        stdout = str(p.stdout.read())
        stderr = str(p.stderr.read())

        if stdout is not None and stdout != '':
            logging.debug('blastn on db {} and query {} STDOUT: {}'.format(query_fasta_path, db_path, stdout))

        if stderr is None or str(stderr) != '' or str(stderr) != "b''" :
            if os.path.exists(blast_outfile):
                return blast_outfile

        if stderr is not None and (str(stderr) != '' and str(stderr) != "b''" ) :
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

    def run_blast(self, query_fasta_path, blast_task, db_path, db_type, min_cov, min_ident, evalue,blast_outfile,logging,num_threads=1,word_size=11,max_target_seqs=100000000,seq_id_file=None):

        if seq_id_file :
            p = Popen(['blastn',
                       '-task', blast_task,
                       '-query', Path(query_fasta_path),
                       '-db', '{}'.format(Path(db_path)),
                       '-seqidlist', '{}'.format(seq_id_file),
                       '-num_threads','{}'.format(num_threads),
                       '-evalue', '{}'.format(evalue),
                       '-dust', 'yes',
                       '-perc_identity', '{}'.format(min_ident),
                       '-max_target_seqs','{}'.format(max_target_seqs),
                       '-out', Path(blast_outfile),
                       '-outfmt', '6 {}'.format(' '.join(BLAST_TABLE_COLS))],
                      stdout=PIPE,
                      stderr=PIPE)
        else:
            p = Popen(['blastn',
                       '-task', blast_task,
                       '-query', '{}'.format(Path(query_fasta_path)),
                       '-db', '{}'.format(Path(db_path)),
                       '-num_threads','{}'.format(num_threads),
                       '-evalue', '{}'.format(evalue),
                       '-dust', 'yes',
                       '-perc_identity', '{}'.format(min_ident),
                       '-max_target_seqs','{}'.format(max_target_seqs),
                       '-out', '{}'.format(Path(blast_outfile)),
                       '-outfmt', '6 {}'.format(' '.join(BLAST_TABLE_COLS))],
                      stdout=PIPE,
                      stderr=PIPE)
        p.wait()
        stdout = str(p.stdout.read())
        stderr = str(p.stderr.read())


        if stdout is not None and stdout != '' and str(stdout) != "b''":
            logging.info('blastn on db {} and query {} STDOUT: {}'.format(query_fasta_path, db_path, stdout))

        if stderr is None and stderr != '' and str(stderr) != "b''" :
            if os.path.exists(blast_outfile):
                return blast_outfile

        if stderr is not None and stderr != '' and stderr != "b''":
            logging.error('blastn on db {} and query {} STDERR: {}'.format(query_fasta_path, db_path, stderr))



class BlastReader:
    df = None


    def __init__(self, blast_outfile,logging):
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
            #self.df = pd.read_table(self.blast_outfile, header=None)
            if not os.path.isfile(blast_outfile):
                logging.warning('No BLASTN results to parse from file %s', blast_outfile)
                self.df = pd.DataFrame()
                return

            if os.path.getsize(blast_outfile) == 0:
                logging.warning('No BLASTN results to parse from file %s', blast_outfile)
                self.df = pd.DataFrame()
                return

            self.df = pd.read_csv(self.blast_outfile,sep="\t",header=None)

            self.df.columns = BLAST_TABLE_COLS

            logger.debug(self.df.head())
            self.is_missing = False


        except EmptyDataError as exc:
            logging.warning('No BLASTN results to parse from file %s', blast_outfile)
            self.is_missing = True
            self.df = pd.DataFrame(index=['A'], columns='A')


    def df_dict(self):
        if not self.is_missing:
            return self.df.to_dict()

