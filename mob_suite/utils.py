from Bio import SeqIO
from Bio.SeqUtils import GC
from mob_suite.blast import BlastRunner
from mob_suite.blast import BlastReader
import os, re
from subprocess import Popen, PIPE, STDOUT
import shutil,sys
import logging

default_database_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'databases')

def check_dependencies(logger):
    external_programs = ['blastn', 'makeblastdb', 'tblastn', 'circlator']
    missing = 0
    for program in external_programs:
        path = shutil.which(program)
        if path is None:
            missing += 1
            logger.error("ERROR: Missing program: {}".format(program,))
        else:
            logger.info("SUCCESS: Found program {} at {}".format(program,path))
    if missing > 0 :
        logger.error("Error, you are missing needed programs for mob-suite, please install them and retry")
        sys.exit(-1)





def fixStart(blast_df):
    for index, row in blast_df.iterrows():
        sstart = blast_df.at[index, 'sstart']
        send = blast_df.at[index, 'send']
        if send < sstart:
            temp = sstart
            blast_df.at[index, 'sstart'] = send
            blast_df.at[index, 'send'] = temp
        qstart = blast_df.at[index, 'qstart']
        qend = blast_df.at[index, 'qend']
        if qend < qstart:
            temp = qstart
            blast_df.at[index, 'qstart'] = qend
            blast_df.at[index, 'qend'] = temp
    return blast_df


def read_fasta_dict(fasta_file):
    seqs = dict()
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqs[str(record.id)] = str(record.seq)
    handle.close()
    return seqs


def write_fasta_dict(seqs, fasta_file):
    with open(fasta_file, "w") as handle:
        for id in seqs:
            handle.write(">{}\n{}\n".format(id, seqs[id]))
    handle.close()


def verify_init(logger, database_dir):
    mob_init_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'mob_init.py')
    status_file = os.path.join(database_dir, 'status.txt')
    if not os.path.isfile(status_file):
        logger.info('MOB-databases need to be initialized, this will take some time')
        p = Popen(['python', mob_init_path, '-d', database_dir],
                  stdout=PIPE,
                  stderr=PIPE,
                  shell=False)

        stdout,stderr = p.communicate()
        print(stdout.decode(), stderr.decode()) #After completion print both streams
        return_code = p.returncode
        logger.info("Return code {}".format(return_code))

        #Verify if no errors were captured during the mob_init script run. Otherwise abort
        if len(re.findall("error",stderr.decode())) != 0 or return_code != 0:
            logger.error(   "Something went wrong with database download or unpacking"
                            "Check MOB-Suite databases directory and your Internet connection.")
            exit(-1)


def filter_overlaping_records(blast_df, overlap_threshold,contig_id_col,contig_start_col,contig_end_col,bitscore_col):
    prev_contig_id = ''
    prev_index = -1
    prev_contig_start = -1
    prev_contig_end = -1
    prev_score = -1
    filter_indexes = list()
    exclude_filter = dict()


    for index, row in blast_df.iterrows():
        contig_id = row['sseqid']
        contig_start = row['sstart']
        contig_end = row['send']
        score = row['bitscore']

        if prev_contig_id == '':
            prev_index = index
            prev_contig_id = contig_id
            prev_contig_start = contig_start
            prev_contig_end = contig_end
            prev_score = score
            continue

        if contig_id != prev_contig_id:
            prev_index = index
            prev_contig_id = contig_id
            prev_contig_start = contig_start
            prev_contig_end = contig_end
            prev_score = score
            continue

        if (contig_start >= prev_contig_start and contig_start <= prev_contig_end) or (contig_end >= prev_contig_start and contig_end <= prev_contig_end):
            overlap = abs(contig_start - prev_contig_end)
            if overlap > overlap_threshold:
                if prev_score > score:
                    filter_indexes.append(index)
                else:
                    filter_indexes.append(prev_index)


        prev_index = index
        prev_contig_id = contig_id
        prev_contig_start = contig_start
        prev_contig_end = contig_end
        prev_score = score

    for index in exclude_filter:
        filter_indexes.append(index)
    indexes = dict()
    for i in blast_df.iterrows():
        indexes[i[0]] = ''

    blast_df.drop(filter_indexes, inplace=True)

    return blast_df.reset_index(drop=True)





def replicon_blast(input_fasta, ref_db, min_ident, min_cov, evalue, tmp_dir,blast_results_file,overlap=5,num_threads=1):
    blast_runner = BlastRunner(input_fasta, tmp_dir)
    blast_runner.makeblastdb(ref_db, 'nucl')
    blast_runner.run_blast(query_fasta_path=input_fasta, blast_task='megablast', db_path=ref_db,
                             db_type='nucl', min_cov=min_cov, min_ident=min_ident, evalue=evalue,
                             blast_outfile=blast_results_file,
                              num_threads=num_threads)
    if os.path.getsize(blast_results_file) == 0:
        return dict()
    blast_df = BlastReader(blast_results_file).df
    blast_df = blast_df.loc[blast_df['pident'] >= min_ident]
    blast_df = blast_df.loc[blast_df['qcovhsp'] >= min_cov]
    blast_df = fixStart(blast_df)
    blast_df = blast_df.sort_values(['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])
    blast_df = blast_df.reset_index(drop=True)
    size = str(len(blast_df))
    blast_df = filter_overlaping_records(blast_df, overlap, 'sseqid', 'sstart', 'send', 'bitscore')
    prev_size = 0
    while size != prev_size:
        blast_df = filter_overlaping_records(blast_df, overlap, 'sseqid', 'sstart', 'send', 'bitscore')
        prev_size = size
        size = str(len(blast_df))

    return blast_df


def mob_blast(input_fasta, ref_db, min_ident, min_cov, evalue, tmp_dir,blast_results_file,overlap=5,num_threads=1):
    num_threads=1
    blast_runner = BlastRunner(input_fasta, tmp_dir)
    blast_runner.makeblastdb(ref_db, 'nucl')
    blast_runner.run_tblastn(query_fasta_path=input_fasta, blast_task='megablast', db_path=ref_db,
                             db_type='nucl', min_cov=min_cov, min_ident=min_ident, evalue=evalue,
                             blast_outfile=blast_results_file,
                              num_threads=num_threads)
    if os.path.getsize(blast_results_file) == 0:
        return dict()
    blast_df = BlastReader(blast_results_file).df
    blast_df = blast_df.loc[blast_df['pident'] >= min_ident]
    blast_df = blast_df.loc[blast_df['qcovhsp'] >= min_cov]
    blast_df = fixStart(blast_df)
    blast_df = blast_df.sort_values(['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])
    blast_df = blast_df.reset_index(drop=True)
    blast_df = filter_overlaping_records(blast_df, overlap, 'sseqid', 'sstart', 'send', 'bitscore')
    prev_size = 0
    size = str(len(blast_df))
    while size != prev_size:
        blast_df = filter_overlaping_records(blast_df, overlap, 'sseqid', 'sstart', 'send', 'bitscore')
        prev_size = size
        size = str(len(blast_df))
    #print(blast_df)
    return blast_df



def repetitive_blast(input_fasta, ref_db, min_ident, min_cov, evalue, min_length, tmp_dir, blast_results_file,num_threads=1):
    blast_runner = BlastRunner(input_fasta, tmp_dir)
    #blast_runner.makeblastdb(ref_db, 'nucl')
    blast_runner.run_blast(query_fasta_path=input_fasta, blast_task='megablast', db_path=ref_db,
                           db_type='nucl', min_cov=min_cov, min_ident=min_ident, evalue=evalue,
                           blast_outfile=blast_results_file,
                              num_threads=num_threads)
    if os.path.getsize(blast_results_file) == 0:
        return dict()

    blast_df = BlastReader(blast_results_file).df
    blast_df = blast_df.loc[blast_df['length'] >= min_length]
    blast_df = blast_df.loc[blast_df['pident'] >= min_ident]
    blast_df = blast_df.loc[blast_df['qcovs'] >= min_cov]
    blast_df = fixStart(blast_df)
    blast_df = blast_df.sort_values(['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])
    blast_df = blast_df.reset_index(drop=True)

    contig_list = dict()
    for index, row in blast_df.iterrows():
        if not row['qseqid'] in contig_list:
            contig_list[row['qseqid']] = {'id': row['sseqid'], 'score': row['bitscore'], 'contig_start': row['sstart'],
                                          'contig_end': row['send']}
        else:
            if contig_list[row['qseqid']]['score'] > row['bitscore']:
                contig_list[row['qseqid']] = {'id': row['sseqid'], 'score': row['bitscore'],
                                              'contig_start': row['sstart'], 'contig_end': row['send']}

    return contig_list


def getRepliconContigs(blast_df):
    contigs = dict()
    if isinstance(blast_df,dict) or blast_df is None:
        return contigs
    for index, row in blast_df.iterrows():
        contig_id = row['sseqid']
        ident = row['pident']
        start = row['sstart']
        end = row['send']
        hit_id = row['qseqid']
        coverage = row['qcovs']
        if not contig_id in contigs:
            contigs[contig_id] = dict()
        contigs[contig_id][hit_id] = {'id': hit_id, 'ident': ident, 'start': start, 'end': end, 'coverage': coverage,
                                      'length': abs(end - start)}
    return contigs


def fix_fasta_header(in_fasta, out_fasta):
    in_basename = os.path.basename(in_fasta)
    fh = open(out_fasta, 'w')
    with open(in_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            fh.write(">" + str(record.description).replace(' ', '_') + "\n" + str(
                record.seq) + "\n")
    handle.close()
    fh.close()


def getMashBestHit(mash_results):
    score = 1
    matches = 0
    top_hit = ''
    top_hit_size = 0
    seqid = ''
    mash_clustid = ''

    for line in mash_results:
        row = line.strip("\n").split("\t")

        if float(score) > float(row[2]):
            seqid, mash_clustid = row[0].split('|')
            score = row[2]
            matches = row[4]

    return {
        'top_hit': seqid,
        'mash_hit_score': score,
        'top_hit_size': top_hit_size,
        'clustid': mash_clustid
    }

''''
    Accepts fasta file and returns size, number of sequence records and gc %
'''


def calcFastaStats(fasta):
    num_seqs = 0;
    seq = ''
    for record in SeqIO.parse(fasta, "fasta"):
        num_seqs += 1
        seq = seq + record.seq
    genome_size = len(seq)
    gc = GC(seq)

    return {
        'num_seq': num_seqs,
        'size': genome_size,
        'gc_content': gc
    }


def init_console_logger(lvl=2):

    LOG_FORMAT = '%(asctime)s %(name)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'

    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    report_lvl = logging_levels[min(lvl, 3)]

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)
    return logging.getLogger(__name__)



