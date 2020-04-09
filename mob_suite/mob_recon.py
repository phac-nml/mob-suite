#!/usr/bin/env python3
from mob_suite.version import __version__
from collections import OrderedDict
import logging, os, shutil, sys, operator,re
import pandas as pd
from argparse import (ArgumentParser, FileType)
from mob_suite.blast import BlastRunner
from mob_suite.blast import BlastReader
from mob_suite.wrappers import mash
from mob_suite.wrappers import detectCircularity

import glob

from mob_suite.constants import  \
    MOB_CLUSTER_INFO_HEADER,\
    MOB_RECON_INFO_HEADER,\
    MOB_TYPER_REPORT_HEADER, \
    ETE3_LOCK_FILE, \
    ETE3DBTAXAFILE, \
    NCBI_PLASMID_TAXONOMY_FILE, \
    NCBI_PLASMID_TAXONOMY_HEADER, \
    default_database_dir,\
    LOG_FORMAT, \
    LIT_PLASMID_TAXONOMY_FILE, \
    LIT_PLASMID_TAXONOMY_HEADER

from mob_suite.mob_typer import \
    hostrange, \
    getAssocValues, \
    filter_invalid_taxids



from mob_suite.utils import \
    read_fasta_dict, \
    write_fasta_dict, \
    filter_overlaping_records, \
    fix_fasta_header, \
    getMashBestHit, \
    verify_init, \
    check_dependencies, \
    read_sequence_info, \
    fixStart, \
    calc_md5, \
    sort_biomarkers, \
    GC, \
    recursive_filter_overlap_records, \
    determine_mpf_type, \
    ETE3_db_status_check, \
    writeReport, \
    dict_from_alt_key_list, \
    read_file_to_dict

def parse_args():
    "Parse the input arguments, use '-h' for help"

    parser = ArgumentParser(
        description="MOB-Recon: Typing and reconstruction of plasmids from draft and complete assemblies version: {}".format(
            __version__))
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output Directory to put results')
    parser.add_argument('-i', '--infile', type=str, required=True, help='Input assembly fasta file to process')
    parser.add_argument('-n', '--num_threads', type=int, required=False, help='Number of threads to be used', default=1)
    parser.add_argument('-s', '--sample_id', type=str, required=False, help='Sample Prefix for reports')
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')
    parser.add_argument('-b', '--filter_db', type=str, required=False, help='Path to fasta file to mask sequences')
    parser.add_argument('-g', '--genome_filter_db_prefix', type=str, required=False,
                        help='Prefix of mash sketch and blastdb of closed chromosomes to use for auto selection of close genomes for filtering')
    parser.add_argument('--mash_genome_neighbor_threshold', type=int, required=False,
                        help='Mash distance selecting valid closed genomes to filter', default=0.002)

    parser.add_argument('--max_contig_size', type=int, required=False,
                        help='Maximum size of a contig to be considered a plasmid',
                        default=310000)
    parser.add_argument('--max_plasmid_size', type=int, required=False,
                        help='Maximum size of a reconstructed plasmid',
                        default=350000)
    parser.add_argument('--min_rep_evalue', type=str, required=False,
                        help='Minimum evalue threshold for replicon blastn',
                        default=0.00001)
    parser.add_argument('--min_mob_evalue', type=str, required=False,
                        help='Minimum evalue threshold for relaxase tblastn',
                        default=0.00001)
    parser.add_argument('--min_con_evalue', type=str, required=False, help='Minimum evalue threshold for contig blastn',
                        default=0.00001)
    parser.add_argument('--min_rpp_evalue', type=str, required=False,
                        help='Minimum evalue threshold for repetitve elements blastn',
                        default=0.00001)

    parser.add_argument('--min_length', type=str, required=False, help='Minimum length of contigs to classify',
                        default=1000)

    parser.add_argument('--min_rep_ident', type=int, required=False, help='Minimum sequence identity for replicons',
                        default=80)
    parser.add_argument('--min_mob_ident', type=int, required=False, help='Minimum sequence identity for relaxases',
                        default=80)
    parser.add_argument('--min_con_ident', type=int, required=False, help='Minimum sequence identity for contigs',
                        default=80)
    parser.add_argument('--min_rpp_ident', type=int, required=False,
                        help='Minimum sequence identity for repetitive elements', default=80)

    parser.add_argument('--min_rep_cov', type=int, required=False,
                        help='Minimum percentage coverage of replicon query by input assembly',
                        default=80)

    parser.add_argument('--min_mob_cov', type=int, required=False,
                        help='Minimum percentage coverage of relaxase query by input assembly',
                        default=80)

    parser.add_argument('--min_con_cov', type=int, required=False,
                        help='Minimum percentage coverage of assembly contig by the plasmid reference database to be considered',
                        default=70)

    parser.add_argument('--min_rpp_cov', type=int, required=False,
                        help='Minimum percentage coverage of contigs by repetitive elements',
                        default=80)

    parser.add_argument('--min_overlap', type=int, required=False,
                        help='Minimum overlap of fragments',
                        default=10)

    parser.add_argument('-u', '--unicycler_contigs', required=False,
                        help='Check for circularity flag generated by unicycler in fasta headers', action='store_true')

    parser.add_argument('-c', '--run_overhang', required=False,
                        help='Detect circular contigs with assembly overhangs', action='store_true')

    parser.add_argument('-k', '--keep_tmp', required=False, help='Do not delete temporary file directory',
                        action='store_true')

    parser.add_argument('-t', '--run_typer', required=False,
                        help='Automatically run Mob-typer on the identified plasmids',
                        action='store_true')

    parser.add_argument('--debug', required=False, help='Show debug information', action='store_true')

    parser.add_argument('--plasmid_db', type=str, required=False, help='Reference Database of complete plasmids',
                        default=os.path.join(default_database_dir,
                                             'ncbi_plasmid_full_seqs.fas'))
    parser.add_argument('--plasmid_mash_db', type=str, required=False,
                        help='Companion Mash database of reference database',
                        default=os.path.join(default_database_dir,
                                             'ncbi_plasmid_full_seqs.fas.msh'))
    parser.add_argument('-m','--plasmid_meta', type=str, required=False,
                        help='MOB-cluster plasmid cluster formatted file matched to the reference plasmid db',
                        default=os.path.join(default_database_dir,
                                             'clusters.txt'))
    parser.add_argument('--plasmid_db_type', type=str, required=False, help='Blast database type of reference database',
                        default='blastn')
    parser.add_argument('--plasmid_replicons', type=str, required=False, help='Fasta of plasmid replicons',
                        default=os.path.join(default_database_dir,
                                             'rep.dna.fas'))
    parser.add_argument('--repetitive_mask', type=str, required=False, help='Fasta of known repetitive elements',
                        default=os.path.join(default_database_dir,
                                             'repetitive.dna.fas'))
    parser.add_argument('--plasmid_mob', type=str, required=False, help='Fasta of plasmid relaxases',
                        default=os.path.join(default_database_dir,
                                             'mob.proteins.faa'))
    parser.add_argument('--plasmid_mpf', type=str, required=False, help='Fasta of known plasmid mate-pair proteins',
                        default=os.path.join(default_database_dir,
                                             'mpf.proteins.faa'))
    parser.add_argument('--plasmid_orit', type=str, required=False, help='Fasta of known plasmid oriT dna sequences',
                        default=os.path.join(default_database_dir,
                                             'orit.fas'))
    parser.add_argument('-d', '--database_directory',
                        default=default_database_dir,
                        required=False,
                        help='Directory you want to use for your databases. If the databases are not already '
                             'downloaded, they will be downloaded automatically. Defaults to {}'.format(default_database_dir))
    parser.add_argument('--primary_cluster_dist', type=int, required=False, help='Mash distance for assigning primary cluster id 0 - 1', default=0.06)
    parser.add_argument('--secondary_cluster_dist', type=int, required=False, help='Mash distance for assigning primary cluster id 0 - 1',
                        default=0.025)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__ )


    return parser.parse_args()


def validate_args(args,logger):

    if not args.outdir:
        logger.error('Error, no output directory specified, please specify one')
        sys.exit(-1)

    if not args.infile:
        logger.error('Error, no fasta specified, please specify one')
        sys.exit(-1)

    if not os.path.isfile(args.infile):
        logger.error('Error, input fasta file does not exist: "{}"'.format(args.infile))
        sys.exit(-1)

        # Input numeric params
        min_rep_ident = float(args.min_rep_ident)
        min_mob_ident = float(args.min_mob_ident)
        min_con_ident = float(args.min_con_ident)
        min_rpp_ident = float(args.min_rpp_ident)

        idents = {'min_rep_ident': min_rep_ident, 'min_mob_ident': min_mob_ident, 'min_con_ident': min_con_ident,
                  'min_rpp_ident': min_rpp_ident}

        for param in idents:
            value = float(idents[param])
            if value < 60:
                logger.error("Error: {} is too low, please specify an integer between 70 - 100".format(param))
                sys.exit(-1)
            if value > 100:
                logger.error("Error: {} is too high, please specify an integer between 70 - 100".format(param))
                sys.exit(-1)

        min_rep_cov = float(args.min_rep_cov)
        min_mob_cov = float(args.min_mob_cov)
        min_con_cov = float(args.min_con_cov)
        min_rpp_cov = float(args.min_rpp_cov)

        covs = {'min_rep_cov': min_rep_cov, 'min_mob_cov': min_mob_cov, 'min_con_cov': min_con_cov,
                'min_rpp_cov': min_rpp_cov}

        for param in covs:
            value = float(covs[param])
            if value < 50:
                logger.error("Error: {} is too low, please specify an integer between 50 - 100".format(param))
                sys.exit(-1)
            if value > 100:
                logger.error("Error: {} is too high, please specify an integer between 50 - 100".format(param))
                sys.exit(-1)

        min_rep_evalue = float(args.min_rep_evalue)
        min_mob_evalue = float(args.min_mob_evalue)
        min_con_evalue = float(args.min_con_evalue)
        min_rpp_evalue = float(args.min_rpp_evalue)

        evalues = {'min_rep_evalue': min_rep_evalue, 'min_mob_evalue': min_mob_evalue, 'min_con_evalue': min_con_evalue,
                   'min_rpp_evalue': min_rpp_evalue}

        for param in evalues:
            value = float(evalues[param])
            if value > 1:
                logger.error("Error: {} is too high, please specify an float evalue between 0 to 1".format(param))
                sys.exit(-1)


        # Input numeric params
        min_rep_ident = float(args.min_rep_ident)
        min_mob_ident = float(args.min_mob_ident)
        min_con_ident = float(args.min_con_ident)
        min_rpp_ident = float(args.min_rpp_ident)

        idents = {'min_rep_ident': min_rep_ident, 'min_mob_ident': min_mob_ident, 'min_con_ident': min_con_ident,
                  'min_rpp_ident': min_rpp_ident}

        for param in idents:
            value = idents[param]
            if value < 70:
                logger.error("Error: {} is too low, please specify an integer between 70 - 100".format(param))
                sys.exit(-1)
            if value > 100:
                logger.error("Error: {} is too high, please specify an integer between 70 - 100".format(param))
                sys.exit(-1)

        min_rep_cov = float(args.min_rep_cov)
        min_mob_cov = float(args.min_mob_cov)
        min_con_cov = float(args.min_con_cov)
        min_rpp_cov = float(args.min_rpp_cov)

        covs = {'min_rep_cov': min_rep_cov, 'min_mob_cov': min_mob_cov, 'min_con_cov': min_con_cov,
                'min_rpp_cov': min_rpp_cov}

        for param in covs:
            value = covs[param]
            if value < 50:
                logger.error("Error: {} is too low, please specify an integer between 50 - 100".format(param))
                sys.exit(-1)
            if value > 100:
                logger.error("Error: {} is too high, please specify an integer between 50 - 100".format(param))
                sys.exit(-1)

        min_rep_evalue = float(args.min_rep_evalue)
        min_mob_evalue = float(args.min_mob_evalue)
        min_con_evalue = float(args.min_con_evalue)
        min_rpp_evalue = float(args.min_rpp_evalue)

        evalues = {'min_rep_evalue': min_rep_evalue, 'min_mob_evalue': min_mob_evalue, 'min_con_evalue': min_con_evalue,
                   'min_rpp_evalue': min_rpp_evalue}

        for param in evalues:
            value = evalues[param]
            if value > 1:
                logger.error("Error: {} is too high, please specify an float evalue between 0 to 1".format(param))
                sys.exit(-1)


def init_console_logger(lvl=2):
    root = logging.getLogger()

    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]

    report_lvl = logging_levels[lvl]
    root.setLevel(report_lvl)  # set root logger level

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)

    return logging.getLogger(__name__)


def circularize(input_fasta, outdir,logging):
    c = detectCircularity()
    return c.run(input_fasta, outdir,logging)


def filter_blastdf_by_seqs(blast_df,include_seqs,column_name):
    blast_df = blast_df[blast_df[column_name].isin(include_seqs)]
    blast_df = blast_df.reset_index(drop=True)
    return blast_df


def blastn(input_fasta, blastdb, min_ident, min_cov, evalue, min_length, out_dir, blast_results_file,logging,seq_filterfile=None,num_threads=1,max_length=400000):

    blast_runner = BlastRunner(input_fasta, out_dir)
    blast_runner.run_blast(query_fasta_path=input_fasta, blast_task='megablast', db_path=blastdb,
                           db_type='nucl', min_cov=min_cov, min_ident=min_ident, evalue=evalue,
                           blast_outfile=blast_results_file,logging=logging, num_threads=num_threads, word_size=11,seq_id_file=seq_filterfile)

    if os.path.getsize(blast_results_file) == 0:
        os.remove(blast_results_file)
        return False

    blast_df = BlastReader(blast_results_file,logging).df
    blast_df = blast_df.loc[blast_df['length'] >= min_length]
    blast_df = blast_df.loc[blast_df['qlen'] <= max_length]
    blast_df = blast_df.loc[blast_df['qcovs'] >= min_cov]
    blast_df = blast_df.loc[blast_df['evalue'] <= evalue]
    blast_df = blast_df.loc[blast_df['pident'] >= min_ident]

    blast_df = blast_df.reset_index(drop=True)
    blast_df = fixStart(blast_df)
    blast_df.to_csv(blast_results_file, sep='\t', header=True, line_terminator='\n', index=False)

    return True


def tblastn(input_fasta, blastdb, min_ident, min_covs, evalue, out_dir,blast_results_file,logging,num_threads=1,min_covhsp=25,seq_id_file=None):
    blast_runner = BlastRunner(input_fasta, out_dir)

    blast_runner.run_tblastn(query_fasta_path=input_fasta, blast_task='megablast', db_path=blastdb,
                             db_type='protein', min_cov=min_covs, min_ident=min_ident, evalue=evalue,
                             blast_outfile=blast_results_file,
                              num_threads=num_threads,seq_id_file=seq_id_file,logging=logging)

    if os.path.getsize(blast_results_file) == 0:
        os.remove(blast_results_file)
        return False

    blast_df = BlastReader(blast_results_file,logging).df
    blast_df = blast_df.loc[blast_df['pident'] >= min_ident]
    blast_df = blast_df.loc[blast_df['qcovs'] >= min_covs]
    blast_df = blast_df.loc[blast_df['qcovhsp'] >= min_covhsp]
    blast_df = blast_df.loc[blast_df['evalue'] <= evalue]
    blast_df = fixStart(blast_df)
    blast_df = blast_df.sort_values(['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])
    blast_df = blast_df.reset_index(drop=True)
    blast_df.to_csv(blast_results_file, sep='\t', header=True, line_terminator='\n', index=False)

    return True

def parseMashScreen(mash_results):
    mash_results = mash_results.decode('utf-8').split("\n")

    hits = {}

    for line in mash_results:
        if len(line) < 4:
            continue
        row = line.strip("\n").split("\t")

        ref_id = str(row[4])
        score = float(row[0])

        if ref_id not in hits:
            hits[ref_id] = { }

        hits[ref_id] = score

    return hits

def parseMash(mash_results):

    mash_results = mash_results.decode('utf-8').split("\n")
    hits = {}

    for line in mash_results:
        row = line.strip("\n").split("\t")
        if len(row) < 4:
            continue
        ref_id = str(row[0])
        query_id = str(row[1])
        score = float(row[2])

        if query_id not in hits:
            hits[query_id] = { }

        hits[query_id][ref_id] = score

    return hits

def filter_sequences(input_fasta,blastdb,min_ident, min_cov, evalue, min_length, out_dir, blast_results_file,seq_filterfile=None,num_threads=1,max_length=400000):
    blastn(input_fasta, blastdb, min_ident, min_cov, evalue, min_length, out_dir, blast_results_file,
           seq_filterfile, num_threads, max_length)

    if os.path.getsize(blast_results_file) == 0:
        os.remove(blast_results_file)
        return pd.DataFrame()

    blast_df = BlastReader(blast_results_file).df

    if seq_filterfile:
        blast_df = filter_blastdf_by_seqs(blast_df,seq_filterfile)
        blast_df = blast_df.reset_index(drop=True)

    return  blast_df


def find_mash_genomes(reference_mash_sketch,fasta_query,outfile,cutoff_distance,num_threads=1):
    cutoff_distance = float(cutoff_distance)

    if (not os.path.isfile(fasta_query)):
        sys.exit(-1)

    if (not os.path.isfile(reference_mash_sketch)):
        sys.exit(-1)
    m = mash()
    distances = parseMash(m.run_mash(reference_mash_sketch,fasta_query,num_threads=num_threads))


    genomes = []

    for query in distances:
        for ref in distances[query]:
            score = distances[query][ref]
            if score < cutoff_distance:
                genomes.append(ref)

    return genomes


def add_biomarker_results(biomarker_df,df_column_name_biomarker,df_column_name_seqid,contig_info,contig_info_type_key,contig_info_acs_key,type_col_num=1,type_acs_num=0,delimeter='|'):
    results = {}
    for index,row in biomarker_df.iterrows():

        biomarker = row[df_column_name_biomarker].split(delimeter)

        if  type_col_num +1  > len(biomarker):
            logging.error("Specified column for biomarker type: {} does not exist in biomarker id field: {}".format(type_col_num,row[df_column_name_biomarker]))
            continue

        if type_acs_num+1 > len(biomarker):
            logging.error("Specified column for biomarker type: {} does not exist in biomarker id field: {}".format(
                type_acs_num, row[df_column_name_biomarker]))
            continue



        seqid = row[df_column_name_seqid]

        if not seqid in results:
            results[seqid] = {'types':[],'acs':[]}

        results[seqid]['types'].append(biomarker[type_col_num])
        results[seqid]['acs'].append(biomarker[type_acs_num])

    results = sort_biomarkers(results)

    for seqid in results:
        if not seqid in contig_info:
            continue

        if contig_info_type_key in contig_info[seqid]:
            contig_info[seqid][contig_info_type_key] = ','.join(results[seqid]['types'])
        else:
            print(
                "Error: {} contig_info_type_key not found in contig_info, check name and fix".format(
                    contig_info_type_key))

        if contig_info_acs_key in contig_info[seqid]:
            contig_info[seqid][contig_info_acs_key] = ','.join(results[seqid]['acs'])
        else:
            print(
                "Error: {} contig_info_acs_key not found in contig_info, check name and fix".format(
                    contig_info_type_key))


def calc_hit_coverage(blast_df,overlap_threshold,reference_sequence_meta):
    blast_df = blast_df.sort_values(['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])
    hit_scores = {}
    size = str(len(blast_df))
    prev_size = 0
    while size != prev_size:
        blast_df = filter_overlaping_records(blast_df, overlap_threshold, 'sseqid', 'sstart', 'send', 'bitscore')
        prev_size = size
        size = str(len(blast_df))

    blast_df['qseqid'].apply(str)
    blast_df['sseqid'].apply(str)

    for index, row in blast_df.iterrows():
        query = str(row['qseqid'])
        pID = str(row['sseqid'])
        score = float(row['bitscore'])
        aln_length = int(row['length'])
        total_len = int(row['slen'])
        if pID not in reference_sequence_meta:
            logging.warning("Seqid {} in blast results but not cluster file".format(pID))
            continue
        else:
            clust_id = reference_sequence_meta[pID]['primary_cluster_id']

        if pID not in hit_scores:
            hit_scores[pID] = {'score':0,'length':total_len,'covered_bases':0,'clust_id':clust_id,'contigs':[]}

        hit_scores[pID]['covered_bases'] += aln_length
        hit_scores[pID]['score'] += score
        hit_scores[pID]['contigs'].append(query)
        hit_scores[pID]['contigs'] = list(set(hit_scores[pID]['contigs']))

    return hit_scores


def calc_contig_reference_cov(blast_df,overlap_threshold,reference_sequence_meta):
    blast_df = blast_df.sort_values(['qseqid', 'sseqid','qstart', 'qend', 'bitscore'], ascending=[True, True, True, True, False])
    contig_scores = {}
    size = str(len(blast_df))
    prev_size = 0

    while size != prev_size:
        blast_df = filter_overlaping_records(blast_df, overlap_threshold, 'sseqid', 'sstart', 'send', 'bitscore')
        prev_size = size
        size = str(len(blast_df))

    blast_df['qseqid'].apply(str)
    blast_df['sseqid'].apply(str)

    for index, row in blast_df.iterrows():
        query = str(row['qseqid'])
        pID = str(row['sseqid'])
        score = float(row['bitscore'])

        if pID not in reference_sequence_meta:
            logging.warning("Seqid {} in blast results but not cluster file".format(pID))
            continue
        else:
            clust_id = reference_sequence_meta[pID]['primary_cluster_id']

        if query not in contig_scores:
            contig_scores[query] = {}

        if not pID in contig_scores[query]:
            contig_scores[query][pID] = 0

        contig_scores[query][pID] += score

    for contig_id in contig_scores:
        contig_scores[contig_id] = OrderedDict(
            sorted(iter(list(contig_scores[contig_id].items())), key=lambda x: x[1], reverse=True))

    return contig_scores


def calc_cluster_scores(reference_hit_coverage):
    cluster_scores = {}
    for ref_id in reference_hit_coverage:
        clust_id = reference_hit_coverage[ref_id]['clust_id']
        if not clust_id in cluster_scores:
            cluster_scores[clust_id] = 0
        score = reference_hit_coverage[ref_id]['score']
        cluster_scores[clust_id]+= score

    return OrderedDict(sorted(iter(list(cluster_scores.items())), key=lambda x: x[1], reverse=True))

def assign_contigs_to_clusters(contig_blast_df,reference_sequence_meta,contig_info,out_dir,contig_seqs,mash_db,primary_distance,secondary_distance,num_threads=1):
    print(contig_info)
    reference_feature_associations = calc_feature_associations(reference_sequence_meta)

    #Individual reference sequence coverage and overall score along with contig associations
    reference_hit_coverage = calc_hit_coverage(contig_blast_df, 1000, reference_sequence_meta)
    contig_reference_coverage = calc_contig_reference_cov(contig_blast_df,1000,reference_sequence_meta)

    group_membership = {}
    replicon_contigs = {}
    relaxase_contigs = {}
    circular_contigs = {}
    repetitive_contigs = {}

    for contig_id in contig_info:

        repetitive = contig_info[contig_id]['repetitive_dna_id']
        if len(repetitive) > 0:
            repetitive_contigs[contig_id] = repetitive

        #Skip contigs which were flagged to be filtered
        if contig_info[contig_id]['filtering_reason'] != 'none':
            continue

        replicon = contig_info[contig_id]['rep_type_accession(s)']
        relaxase = contig_info[contig_id]['relaxase_type_accession(s)']
        circular_status = contig_info[contig_id]['circularity_status']



        if not contig_id in replicon_contigs:
            if replicon != '-' and replicon != '':
                replicon_contigs[contig_id] = replicon.split(',')

        if not contig_id in relaxase_contigs:

            if relaxase != '-' and relaxase != '':
                relaxase_contigs[contig_id] = relaxase.split(',')

        if circular_status == 'circular' and (contig_id in replicon_contigs or contig_id in relaxase_contigs):
            circular_contigs[contig_id] = "circular"


    #Use circular contigs with replicon or relaxase as the initial seed for group memberships
    for contig_id in circular_contigs:

        if contig_id in replicon_contigs or contig_id in relaxase_contigs:

            if contig_id not in group_membership:
                group_membership[contig_id] = {
                    'clust_id': None,
                    'score': 0,
                    'is_circular': True,
                    'contains_replicon': False,
                    'contains_relaxase': False,
                    'rep_type': '',
                    'mob_type': ''
                }

            if contig_id in replicon_contigs:
                group_membership[contig_id]['rep_type'] = replicon_contigs[contig_id]
                group_membership[contig_id]['contains_replicon'] = True

            if contig_id in relaxase_contigs:
                group_membership[contig_id]['mob_type'] = relaxase_contigs[contig_id]
                group_membership[contig_id]['contains_relaxase'] = True

    contig_blast_df = contig_blast_df[contig_blast_df.qseqid.isin(list(group_membership.keys()))]
    contig_blast_df.reset_index(drop=True)

    cluster_contig_links = get_seq_links(reference_hit_coverage)
    cluster_scores = calc_cluster_scores(reference_hit_coverage)

    contig_link_counts = {}
    contig_clust_assoc = {}
    for clust_id in cluster_contig_links:
        contigs = cluster_contig_links[clust_id]
        for contig_id in contigs:
            if not contig_id in contig_link_counts:
                contig_link_counts[contig_id] = 0
                contig_clust_assoc[contig_id] = {}
            contig_link_counts[contig_id]+=1
            contig_clust_assoc[contig_id][clust_id] = cluster_scores[clust_id]

    for contig_id in contig_clust_assoc:
        contig_clust_assoc[contig_id] = OrderedDict(sorted(iter(contig_clust_assoc[contig_id].items()), key=lambda x: x[1], reverse=True))

    contig_link_counts = OrderedDict(sorted(iter(contig_link_counts.items()), key=lambda x: x[1], reverse=False))




    for contig_id in contig_link_counts:
        clust_id = next(iter(contig_clust_assoc[contig_id]))
        for c_id in cluster_contig_links[clust_id]:
            print("{}\t{}\t{}\t{}".format(contig_id,c_id,clust_id,cluster_scores[clust_id]))
            if contig_info[c_id]['filtering_reason'] in ['chromosome', 'user filter']:
                continue
            contig_clust_id = contig_info[c_id]['primary_cluster_id']
            if contig_clust_id != '':
                continue
            contig_info[c_id]['primary_cluster_id'] = clust_id
            contig_info[c_id]['molecule_type'] = 'plasmid'

    cluster_prioritzation = {}
    print('-------->')
    for clust_id in cluster_contig_links:
        if not clust_id in cluster_prioritzation:
            cluster_prioritzation[clust_id] = {
                'score': 0,
                'is_circular':False,
                'num_contigs':0,
                'contains_replicon': False,
                'contains_relaxase': False,
                'contigs': [],
            }

        cluster_prioritzation[clust_id]['contigs'] = cluster_contig_links[clust_id]
        cluster_prioritzation[clust_id]['num_contigs'] = len(cluster_contig_links[clust_id])
        cluster_prioritzation[clust_id]['score'] = cluster_scores[clust_id]

        if cluster_prioritzation[clust_id]['num_contigs'] == 1:
            contig_id = next(iter(cluster_prioritzation[clust_id]['contigs']))

            if contig_id in circular_contigs:
                cluster_prioritzation[clust_id]['is_circular'] = True

        if contig_id in replicon_contigs:
            cluster_prioritzation[clust_id]['contains_replicon'] = True

        if contig_id in relaxase_contigs:
            cluster_prioritzation[clust_id]['contains_relaxase'] = True

    for clust_id in cluster_scores:
        print("{}\t{}".format(clust_id,cluster_prioritzation[clust_id]))

    unassigned_contigs = {}

    for contig_id in contig_info:
        if contig_info[contig_id]['filtering_reason'] in ['chromosome', 'user filter']:
            continue

        if contig_info[contig_id]['primary_cluster_id']  != '':
            unassigned_contigs[contig_id] = ''

    #prioritize replicon/relaxase containing clusters
    for clust_id in cluster_scores:
        if not cluster_prioritzation[clust_id]['contains_replicon'] and not cluster_prioritzation[clust_id]['contains_relaxase']:
            continue

        contigs = cluster_contig_links[clust_id]
        # Skip contigs which were flagged to be filtered
        if contig_info[contig_id]['filtering_reason'] in ['chromosome', 'user filter']:
            continue

        for contig_id in contigs:

            # Skip contigs which were flagged to be filtered
            if contig_info[contig_id]['filtering_reason'] in ['chromosome','user filter'] :
                continue

            contig_clust_id = contig_info[contig_id]['primary_cluster_id']
            if len(contig_clust_id) > 0:
                continue
            contig_info[contig_id]['primary_cluster_id'] = clust_id
            contig_info[contig_id]['molecule_type'] = 'plasmid'


    for clust_id in cluster_scores:
        if clust_id in cluster_contig_links:
            contigs = cluster_contig_links[clust_id]
            print("{}\t{}".format(clust_id,contigs))
            for contig_id in contigs:

                # Skip contigs which were flagged to be filtered
                if contig_info[contig_id]['filtering_reason'] in ['chromosome','user filter'] :
                    continue

                contig_clust_id = contig_info[contig_id]['primary_cluster_id']
                if len(contig_clust_id) > 0:
                    continue
                contig_info[contig_id]['primary_cluster_id'] = clust_id
                contig_info[contig_id]['molecule_type'] = 'plasmid'

    cluster_links = {}
    print(contig_info)
    for contig_id in contig_info:
        clust_id = contig_info[contig_id]['primary_cluster_id']
        if len(clust_id) == 0:
            continue
        if not clust_id in cluster_links:
            cluster_links[clust_id] = []

        cluster_links[clust_id].append(contig_id)

    recon_cluster_dists = get_reconstructed_cluster_dists(mash_db,0.1,cluster_links,out_dir,contig_seqs,num_threads)

    print(recon_cluster_dists)

    #get lowest distance cluster
    counter = 0
    increment = False
    for clust_id in recon_cluster_dists:
        fail = False
        for top_ref_id in recon_cluster_dists[clust_id]:
            lowest_dist = recon_cluster_dists[clust_id][top_ref_id]
            if top_ref_id not in reference_sequence_meta:
                fail = True
                continue
            else:
                fail = False
                break

        if fail:
            continue

        contained_replicons = list(set(list(replicon_contigs.keys())) & set(cluster_links[clust_id]))
        contained_relaxases = list(set(list(relaxase_contigs.keys())) & set(cluster_links[clust_id]))
        contained_repettive = list(set(list(repetitive_contigs.keys())) & set(cluster_links[clust_id]))

        #if lowest_dist > primary_distance:

         #   if len(contained_replicons) == 0 and len(contained_relaxases) == 0:
          #      return


        if increment:
            counter+=1
            increment = False

        for contig_id in cluster_links[clust_id]:
            #skip clusters which are just repetitive elemenets
            if len(contained_repettive) == len(cluster_links[clust_id]):
                contig_info[contig_id]['primary_cluster_id'] = ''
                contig_info[contig_id]['molecule_type'] = 'chromosome'
                continue

            if lowest_dist <= primary_distance:
                contig_info[contig_id]['primary_cluster_id'] = clust_id
                contig_info[contig_id]['molecule_type'] = 'plasmid'
                contig_info[contig_id]['mash_nearest_neighbor'] = top_ref_id
                contig_info[contig_id]['mash_neighbor_distance'] = lowest_dist
                contig_info[contig_id]['mash_neighbor_identification'] = reference_sequence_meta[top_ref_id]['organism']

                if lowest_dist <= secondary_distance:
                    contig_info[contig_id]['secondary_cluster_id'] = reference_sequence_meta[top_ref_id]['secondary_cluster_id']
            else:
                if (len(contained_replicons) > 0 or len(contained_relaxases) > 0) and len(contained_repettive) != len(cluster_links[clust_id]):
                    contig_info[contig_id]['primary_cluster_id'] = "novel_{}".format(counter)
                    increment = True
                    contig_info[contig_id]['molecule_type'] = 'plasmid'
                    contig_info[contig_id]['mash_nearest_neighbor'] = top_ref_id
                    contig_info[contig_id]['mash_neighbor_distance'] = lowest_dist
                    contig_info[contig_id]['mash_neighbor_identification'] = reference_sequence_meta[top_ref_id]['organism']
                else:
                    contig_info[contig_id]['primary_cluster_id'] = ''
                    contig_info[contig_id]['molecule_type'] = 'chromosome'
    print(contig_info)
    return contig_info


def get_reconstructed_cluster_dists(mash_db,mash_distance,cluster_contig_links,out_dir,contig_seqs,num_threads=1):
    m = mash()
    cluster_dists = {}
    for clust_id in cluster_contig_links:
        contigs = cluster_contig_links[clust_id]
        seq_dict = {}
        tmp_fasta = os.path.join(out_dir,"clust_{}.fasta".format(clust_id))

        for contig_id in contigs:
            if contig_id in contig_seqs:
                seq_dict[contig_id] = contig_seqs[contig_id]
        write_fasta_dict(seq_dict, tmp_fasta)

        distances = parseMash(m.run_mash(reference_db=mash_db, input_fasta=tmp_fasta,  table=False,num_threads=num_threads))
        os.remove(tmp_fasta)



        for query in distances:
            for ref in distances[query]:
                score = float(distances[query][ref])
                if score <= mash_distance:
                    if clust_id not in cluster_dists:
                        cluster_dists[clust_id] = {}
                    cluster_dists[clust_id][ref] = score

        for clust_id in cluster_dists:
            cluster_dists[clust_id] =  OrderedDict(sorted(iter(list(cluster_dists[clust_id].items())), key=lambda x: x[1], reverse=False))

    return cluster_dists


def get_seq_links(reference_hit_coverage):
    reference_clust_members = {}
    for ref_id in reference_hit_coverage:
        contigs = reference_hit_coverage[ref_id]['contigs']
        clust_id = reference_hit_coverage[ref_id]['clust_id']

        if not clust_id in reference_clust_members:
            reference_clust_members[clust_id] = {}

        for contig_id in contigs:
            if contig_id not in reference_clust_members[clust_id]:
                reference_clust_members[clust_id][contig_id] = 0

            reference_clust_members[clust_id][contig_id]+=1

    return reference_clust_members



def update_group_members(target_contigs,group_membership,contig_reference_coverage,reference_sequence_meta,group_membership_key,reference_seq_key):

    for contig_id in target_contigs:
        types = target_contigs[contig_id]


        if not contig_id in group_membership:
            group_membership[contig_id] = {
                'clust_id': None,
                'score': 0,
                'is_circular': False,
                'contains_replicon': False,
                'contains_relaxase': False,
                'rep_type': '',
                'mob_type': ''
            }

        group_membership[contig_id][group_membership_key] = target_contigs[contig_id]

        contig_hit_scores = contig_reference_coverage[contig_id]

        if group_membership[contig_id]['clust_id'] is not None:
            continue

        for hsp in contig_hit_scores:
            hsp_score = contig_hit_scores[hsp]

            if hsp not in reference_sequence_meta:
                continue

            clust_id = reference_sequence_meta[hsp]['primary_cluster_id']

            if reference_sequence_meta[hsp][reference_seq_key] == '-':
                hsp_rep_types = []
            else:
                hsp_rep_types = reference_sequence_meta[hsp][reference_seq_key].split(",")

            for r in types:
                if r in hsp_rep_types:
                    group_membership[contig_id]['clust_id'] = clust_id
                    break

    return group_membership

def filter_contig_df_by_index(indicies,contig_blast_df,reference_hit_coverage):
    for index in indicies:

        row = contig_blast_df.iloc[index]
        query = str(row['qseqid'])
        pID = str(row['sseqid'])
        score = float(row['bitscore'])
        aln_length = int(row['length'])
        total_len = int(row['slen'])

        if pID not in reference_hit_coverage:
            logging.warning("Seqid {} in blast results but not cluster file".format(pID))
            continue

        if pID in reference_hit_coverage:

            reference_hit_coverage[pID]['score'] -= score
            reference_hit_coverage[pID]['covered_bases'] -= aln_length
        else:
            print("{} not found".format(pID))

    return reference_hit_coverage


def calc_feature_associations(reference_sequence_meta):
    replicon_relaxase = {}
    cluster_replicon = {}
    cluster_relaxase = {}

    replicon_counts = {}
    relaxase_counts = {}
    cluster_counts = {}

    for id in reference_sequence_meta:
        if 'rep_type_accession(s)' in reference_sequence_meta[id]:
            rep_type = reference_sequence_meta[id]['rep_type_accession(s)']
            if rep_type == '' or rep_type == '-':
                rep_type = []
            else:
                rep_type = rep_type.split(",")
        else:
            rep_type = []

        if 'relaxase_type_accession(s)' in reference_sequence_meta[id]:
            mob_type = reference_sequence_meta[id]['relaxase_type_accession(s)']
            if mob_type == '' or mob_type == '-':
                mob_type = []
            else:
                mob_type = mob_type.split(",")
        else:
            mob_type = []

        if 'primary_cluster_id' in reference_sequence_meta[id]:
            clust_id = reference_sequence_meta[id]['primary_cluster_id']
            if clust_id == '' or clust_id == '-':
                clust_id = None

        else:
            clust_id = None

        if clust_id is None:
            continue

        if clust_id not in cluster_counts:
            cluster_counts[clust_id] = 0

        cluster_counts[clust_id]+= 1

        if clust_id not in cluster_replicon:
            cluster_replicon[clust_id] = {}

        if clust_id not in cluster_relaxase:
            cluster_relaxase[clust_id] = {}

        for r in rep_type:
            if r not in replicon_counts:
                replicon_counts[r] = 0

            replicon_counts[r]+= 1

            if r not in cluster_replicon[clust_id]:
                cluster_replicon[clust_id][r] = 0

            cluster_replicon[clust_id][r] += 1

            if len(mob_type) > 0:
                if r not in replicon_relaxase:
                    replicon_relaxase[r] = {}
                for m in mob_type:
                    if m not in replicon_relaxase[r]:
                        replicon_relaxase[r][m] = 0
                    replicon_relaxase[r][m] += 1

        for m in mob_type:
            if m not in relaxase_counts:
                relaxase_counts[m] = 0
            relaxase_counts[m] += 1
            if m not in cluster_relaxase[clust_id]:
                cluster_relaxase[clust_id][m] = 0
            cluster_relaxase[clust_id][m] += 1

    #sort the associations by probability of co-occurance
    for r in replicon_relaxase:
        total_rep = replicon_counts[r]
        for m in replicon_relaxase[r]:
            total_mob = relaxase_counts[m]
            count = replicon_relaxase[r][m]
            replicon_relaxase[r][m] = count/(total_mob + total_rep)
        replicon_relaxase[r] = OrderedDict(sorted(iter(list(replicon_relaxase[r].items())), key=lambda x: x[1], reverse=True))

    for clust_id in cluster_replicon:
        total_cluster = cluster_counts[clust_id]
        for r in cluster_replicon[clust_id]:
            total_rep = replicon_counts[r]
            count = cluster_replicon[clust_id][r]
            cluster_replicon[clust_id][r] = count / (total_cluster + total_rep)
        cluster_replicon[clust_id][r] = OrderedDict(sorted(iter(list(cluster_replicon[clust_id].items())), key=lambda x: x[1], reverse=True))

    for clust_id in cluster_relaxase:
        total_cluster = cluster_counts[clust_id]
        for m in cluster_relaxase[clust_id]:
            total_mob = relaxase_counts[m]
            count = cluster_relaxase[clust_id][m]
            cluster_relaxase[clust_id][m] = count / (total_cluster + total_mob)
        cluster_relaxase[clust_id] = OrderedDict(
            sorted(iter(list(cluster_relaxase[clust_id].items())), key=lambda x: x[1], reverse=True))

    return {'cluster_relaxase':cluster_relaxase, 'cluster_replicon':cluster_replicon, 'replicon_relaxase':replicon_relaxase}


def build_mobtyper_report(plasmid_contig_info,out_dir,outfile,seq_dict,ncbi,lit):
    mob_typer_results = {}
    for clust_id in plasmid_contig_info:

        cluster_file = open(os.path.join(out_dir,"plasmid_{}.fasta".format(clust_id)),'w')
        logging.info("Writting plasmid sequences to {}".format(os.path.join(out_dir,"plasmid_{}.fasta".format(clust_id))))
        if clust_id not in mob_typer_results:
            mob_typer_results[clust_id] = {}
        for field in MOB_TYPER_REPORT_HEADER:
            if not field in mob_typer_results[clust_id]:
                mob_typer_results[clust_id][field] = []

        data = plasmid_contig_info[clust_id]

        #Put contig report data into MOB-typer report header
        #aggregating the data by cluster id
        cluster_seq = []
        for contig_id in data:
            if contig_id in seq_dict:
                cluster_file.write(">{}\n{}\n".format(contig_id,seq_dict[contig_id]))
                cluster_seq.append(seq_dict[contig_id])
            for field in MOB_TYPER_REPORT_HEADER:
                if field in data[contig_id]:
                    if isinstance(data[contig_id][field],list) and len(data[contig_id][field]) == 0:
                        continue
                    if data[contig_id][field] == '':
                        continue
                    mob_typer_results[clust_id][field].append(data[contig_id][field])


        #overwrite individual seq stat calculations with the overall
        seq = "".join(cluster_seq)
        mob_typer_results[clust_id]['md5'] = [calc_md5(seq)]
        mob_typer_results[clust_id]['gc'] = [GC(seq)]
        mob_typer_results[clust_id]['size'] = [len(seq)]
        mob_typer_results[clust_id]['num_contigs'] = len(cluster_seq)

        #Sort MOB-typer biomarker results
        replicon = sort_biomarkers({'rep':{'types':mob_typer_results[clust_id]['rep_type(s)'],'acs':mob_typer_results[clust_id]['rep_type_accession(s)']}})
        mob_typer_results[clust_id]['rep_type(s)'] = ",".join(replicon['rep']['types'])
        mob_typer_results[clust_id]['rep_type_accession(s)'] = ",".join(replicon['rep']['acs'])

        relaxase = sort_biomarkers({'mob':{'types':mob_typer_results[clust_id]['relaxase_type(s)'],'acs':mob_typer_results[clust_id]['relaxase_type_accession(s)']}})
        mob_typer_results[clust_id]['relaxase_type(s)'] = ",".join(relaxase['mob']['types'])
        mob_typer_results[clust_id]['relaxase_type_accession(s)'] = ",".join(relaxase['mob']['acs'])
        if len(mob_typer_results[clust_id]['mpf_type']) > 0:
            mob_typer_results[clust_id]['mpf_type'] = determine_mpf_type(mob_typer_results[clust_id]['mpf_type'])
        else:
            mob_typer_results[clust_id]['mpf_type'] = ''


        mob_typer_results[clust_id]['mpf_type_accession(s)'] = ",".join(mob_typer_results[clust_id]['mpf_type_accession(s)'])

        mob_typer_results[clust_id]['orit_type(s)'] = ",".join(mob_typer_results[clust_id]['orit_type(s)'])
        mob_typer_results[clust_id]['orit_accession(s)'] = ",".join(mob_typer_results[clust_id]['orit_accession(s)'])

        #Assign mobility
        mob_typer_results[clust_id]['predicted_mobility'] = 'non-mobilizable'
        if len(mob_typer_results[clust_id]['relaxase_type(s)']) > 0 and len(mob_typer_results[clust_id]['mpf_type']) > 0:
            mob_typer_results[clust_id]['predicted_mobility'] = 'conjugative'
        elif(len(mob_typer_results[clust_id]['relaxase_type(s)']) > 0 or len(mob_typer_results[clust_id]['orit_type(s)']) > 0):
            mob_typer_results[clust_id]['predicted_mobility'] = 'mobilizable'

        if isinstance(mob_typer_results[clust_id]['primary_cluster_id'],list) and len(mob_typer_results[clust_id]['primary_cluster_id']) >0:
            mob_typer_results[clust_id]['primary_cluster_id'] = mob_typer_results[clust_id]['primary_cluster_id'][0]
            mob_typer_results[clust_id]['mash_nearest_neighbor'] = mob_typer_results[clust_id]['mash_nearest_neighbor'][0]
            mob_typer_results[clust_id]['mash_neighbor_distance'] = mob_typer_results[clust_id]['mash_neighbor_distance'][0]
            mob_typer_results[clust_id]['mash_neighbor_identification'] = mob_typer_results[clust_id]['mash_neighbor_identification'][0]



        if isinstance(mob_typer_results[clust_id]['sample_id'],list) and len(mob_typer_results[clust_id]['sample_id']) > 0:
            mob_typer_results[clust_id]['sample_id'] = "{}:{}".format(mob_typer_results[clust_id]['sample_id'][0],mob_typer_results[clust_id]['primary_cluster_id'])

        if isinstance(mob_typer_results[clust_id]['secondary_cluster_id'], list) and len(mob_typer_results[clust_id]['secondary_cluster_id']) > 0:
            mob_typer_results[clust_id]['secondary_cluster_id'] = mob_typer_results[clust_id]['secondary_cluster_id'][0]


        if mob_typer_results[clust_id]['num_contigs'] > 1:
            mob_typer_results[clust_id]['circularity_status'] = 'incomplete'

        if len(replicon['rep']['types']) > 0:
            rep_types = replicon['rep']['types']
        else:
            rep_types = []

        if len(relaxase['mob']['acs']) > 0:
            relaxase_types = relaxase['mob']['acs']
        else:
            relaxase_types = []

        if mob_typer_results[clust_id]['primary_cluster_id'] != '' and not 'novel_' in mob_typer_results[clust_id]['primary_cluster_id']:
            mob_cluster_id = mob_typer_results[clust_id]['primary_cluster_id']
        else:
            mob_cluster_id = '-'

        host_range = hostrange(rep_types, relaxase_types, mob_cluster_id, ncbi, lit)

        for field in host_range:
            mob_typer_results[clust_id][field] = host_range[field]

        for element in MOB_TYPER_REPORT_HEADER:
            if element in mob_typer_results[clust_id]:
                if isinstance(mob_typer_results[clust_id][element],list):
                    mob_typer_results[clust_id][element] = ', '.join(str(x) for x in mob_typer_results[clust_id][element])

    results = []
    for clust_id in mob_typer_results:
        results.append( mob_typer_results[clust_id])

    writeReport(results,MOB_TYPER_REPORT_HEADER,outfile)
    return




def main():
    args = parse_args()

    if args.debug:
        logger = init_console_logger(3)
    else:
        logger = init_console_logger(2)

    logger.info("MOB-recon version {} ".format(__version__))
    logger.debug("Debug log reporting set on successfully")

    check_dependencies(logger)
    validate_args(args,logger)

    keep_tmp = args.keep_tmp
    plasmid_files = []
    input_fasta = args.infile
    out_dir = args.outdir
    num_threads = args.num_threads
    tmp_dir = os.path.join(out_dir, '__tmp')
    file_id = os.path.basename(input_fasta)
    fixed_fasta = os.path.join(tmp_dir, 'fixed.input.fasta')
    chromosome_file = os.path.join(out_dir, 'chromosome.fasta')
    replicon_blast_results = os.path.join(tmp_dir, 'replicon_blast_results.txt')
    mob_blast_results = os.path.join(tmp_dir, 'mob_blast_results.txt')
    mpf_blast_results = os.path.join(tmp_dir, 'mpf_blast_results.txt')
    orit_blast_results = os.path.join(tmp_dir, 'orit_blast_results.txt')
    repetitive_blast_results = os.path.join(tmp_dir, 'repetitive_blast_results.txt')
    contig_blast_results = os.path.join(tmp_dir, 'contig_blast_results.txt')
    contig_report = os.path.join(out_dir, 'contig_report.txt')

    logger.info('Processing fasta file {}'.format(args.infile))
    logger.info('Analysis directory {}'.format(args.outdir))

    database_dir = os.path.abspath(args.database_directory)

    if database_dir == default_database_dir:
        plasmid_ref_db = args.plasmid_db
        mob_ref = args.plasmid_mob
        mash_db = args.plasmid_mash_db
        replicon_ref = args.plasmid_replicons
        plasmid_meta = args.plasmid_meta
        repetitive_mask_file = args.repetitive_mask
        mpf_ref = args.plasmid_mpf
        plasmid_orit = args.plasmid_orit
    else:
        plasmid_ref_db = os.path.join(database_dir, 'ncbi_plasmid_full_seqs.fas')
        mob_ref = os.path.join(database_dir, 'mob.proteins.faa')
        mash_db = os.path.join(database_dir, 'ncbi_plasmid_full_seqs.fas.msh')
        replicon_ref = os.path.join(database_dir, 'rep.dna.fas')
        plasmid_meta = os.path.join(database_dir, 'clusters.txt')
        repetitive_mask_file = os.path.join(database_dir, 'repetitive.dna.fas')
        mpf_ref = os.path.join(database_dir, 'mpf.proteins.faa')
        plasmid_orit = os.path.join(database_dir, 'orit.fas')

    if args.sample_id is None:
        sample_id = re.sub(r"\.(fasta|fas|fa){1,1}", "", os.path.basename(args.infile))
    else:
        sample_id = args.sample_id

    run_overhang = args.run_overhang
    unicycler_contigs = args.unicycler_contigs

    #initialize analysis directory
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir, 0o755)

    elif not args.force:
        logger.error("Error output directory exists, please specify a new directory or use --force to overwrite")
        sys.exit(-1)
    else:
        shutil.rmtree(args.outdir)
        os.mkdir(args.outdir, 0o755)

    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir, 0o755)
    else:
        shutil.rmtree(tmp_dir)
        os.mkdir(tmp_dir, 0o755)

    #Initialize clustering distance thresholds
    if not (args.primary_cluster_dist >= 0 and args.primary_cluster_dist <= 1):
        logging.error('Error distance thresholds must be between 0 - 1: {}'.format(args.primary_cluster_dist))
        sys.exit()
    else:
        primary_distance = float(args.primary_cluster_dist)

    if not (args.secondary_cluster_dist >= 0 and args.secondary_cluster_dist <= 1):
        logging.error('Error distance thresholds must be between 0 - 1: {}'.format(args.secondary_cluster_dist))
        sys.exit()
    else:
        secondary_distance = float(args.secondary_cluster_dist)

    # Input numeric params

    min_overlapp = int(args.min_overlap)
    min_length = int(args.min_length)

    min_rep_ident = float(args.min_rep_ident)
    min_mob_ident = float(args.min_mob_ident)
    min_con_ident = float(args.min_con_ident)
    min_rpp_ident = float(args.min_rpp_ident)

    min_rep_cov = float(args.min_rep_cov)
    min_mob_cov = float(args.min_mob_cov)
    min_con_cov = float(args.min_con_cov)
    min_rpp_cov = float(args.min_rpp_cov)

    min_rep_evalue = float(args.min_rep_evalue)
    min_mob_evalue = float(args.min_mob_evalue)
    min_con_evalue = float(args.min_con_evalue)
    min_rpp_evalue = float(args.min_rpp_evalue)

    #Test that ETE3 db is ok and lock process check
    dbstatus = ETE3_db_status_check(1, ETE3_LOCK_FILE, ETE3DBTAXAFILE, logging)
    if dbstatus == False:
        logging.error("Exiting due to lock file not removed: {}".format(ETE3_LOCK_FILE))
        sys.exit(-1)

    #Parse reference cluster information
    reference_sequence_meta = read_sequence_info(plasmid_meta,MOB_CLUSTER_INFO_HEADER)


    #process input fasta
    logger.info('Writing cleaned header input fasta file from {} to {}'.format(input_fasta, fixed_fasta))
    fix_fasta_header(input_fasta, fixed_fasta)
    contig_seqs = read_fasta_dict(fixed_fasta)
    contig_info = {}
    br = BlastRunner(fixed_fasta, tmp_dir)
    br.makeblastdb(fixed_fasta, dbtype='nucl',logging=logging)
    del(br)



    #Detect circular sequences
    circular_contigs = {}
    if run_overhang:
        logger.info('Running internal circular contig detection on {}'.format(fixed_fasta))
        circular_contigs = circularize(fixed_fasta, tmp_dir,logging)


    for id in contig_seqs :
        seq = contig_seqs[id]
        contig_info[id] = {}
        for feature in MOB_RECON_INFO_HEADER:
            contig_info[id][feature] = ''
        contig_info[id]['md5'] = calc_md5(seq)
        contig_info[id]['gc'] = GC(seq)
        contig_info[id]['size'] = len(seq)
        contig_info[id]['contig_id'] = id
        contig_info[id]['sample_id'] = sample_id
        contig_info[id]['molecule_type'] = 'chromosome'
        contig_info[id]['filtering_reason'] = 'none'

        if run_overhang:
            if id in circular_contigs:
                contig_info[id]['circularity_status'] = 'circular'
            else:
                contig_info[id]['circularity_status'] = 'incomplete'

        if unicycler_contigs:
            if 'circular=true' in id or '_circ' in id:
                contig_info[id]['circularity_status'] = 'circular'
            elif id not in circular_contigs:
                contig_info[id]['circularity_status'] = 'incomplete'

        if contig_info[id]['circularity_status'] == '':
            contig_info[id]['circularity_status'] = 'not tested'




    #Blast reference databases

    #blast replicon database
    logging.info("Blasting replicon sequences {} against {}".format(replicon_ref,fixed_fasta))
    blastn(input_fasta=replicon_ref,blastdb=fixed_fasta,min_ident=min_rep_ident,min_cov=min_rep_cov,evalue=min_rep_evalue,min_length=25,out_dir=tmp_dir,
           blast_results_file=replicon_blast_results,num_threads=num_threads,logging=logging)

    logging.info("Filtering replicon blast results {} ".format(replicon_blast_results))
    rep_blast_df = BlastReader(replicon_blast_results, logging=logging).df
    if len(rep_blast_df) > 0:
        rep_blast_df = rep_blast_df.drop(0)
        rep_blast_df = fixStart(rep_blast_df).sort_values(['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])
        rep_blast_df = recursive_filter_overlap_records(rep_blast_df, 5, 'sseqid', 'sstart', 'send',
                                  'bitscore')



        add_biomarker_results(biomarker_df=rep_blast_df, df_column_name_biomarker='qseqid', df_column_name_seqid='sseqid', contig_info=contig_info,
                              contig_info_type_key='rep_type(s)', contig_info_acs_key='rep_type_accession(s)', delimeter='|')


    del(rep_blast_df)

    #blast relaxase database
    logging.info("Blasting relaxase sequences {} against {}".format(mob_ref, fixed_fasta))
    tblastn(input_fasta=mob_ref, blastdb=fixed_fasta, min_ident=min_mob_ident, min_covs=min_mob_cov, evalue=min_mob_evalue, out_dir=tmp_dir,logging=logging,
            blast_results_file=mob_blast_results, num_threads=num_threads)

    logging.info("Filtering relaxase blast results {} ".format(mob_blast_results))

    mob_blast_df = BlastReader(mob_blast_results,logging).df
    if len(mob_blast_df) > 0:
        mob_blast_df = fixStart(mob_blast_df.drop(0).sort_values(['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False]))
        mob_blast_df = recursive_filter_overlap_records(mob_blast_df, 5, 'sseqid', 'sstart', 'send',
                                  'bitscore')

        add_biomarker_results(biomarker_df=mob_blast_df, df_column_name_biomarker='qseqid', df_column_name_seqid='sseqid', contig_info=contig_info,
                              contig_info_type_key='relaxase_type(s)', contig_info_acs_key='relaxase_type_accession(s)', delimeter='|')


    del (mob_blast_df)


    #blast mpf database
    logging.info("Blasting MPF sequences {} against {}".format(mpf_ref, fixed_fasta))
    tblastn(input_fasta=mpf_ref, blastdb=fixed_fasta, min_ident=min_mob_ident, min_covs=min_mob_cov, evalue=min_mob_evalue, out_dir=tmp_dir,
            blast_results_file=mpf_blast_results, num_threads=num_threads,logging=logging)

    mpf_blast_df = BlastReader(mpf_blast_results, logging).df

    if len(mpf_blast_df) > 0:
        mpf_blast_df = fixStart(mpf_blast_df.drop(0)).sort_values(['sseqid', 'sstart', 'send', 'bitscore'],
                                                          ascending=[True, True, True, False])
        mpf_blast_df = recursive_filter_overlap_records(mpf_blast_df, 5, 'sseqid', 'sstart', 'send',
                                  'bitscore')

        logging.info("Filtering MPF blast results {} ".format(mpf_blast_results))

        add_biomarker_results(biomarker_df=mpf_blast_df, df_column_name_biomarker='qseqid', df_column_name_seqid='sseqid', contig_info=contig_info,
                              contig_info_type_key='mpf_type', contig_info_acs_key='mpf_type_accession(s)', delimeter='|')
    del(mpf_blast_results)

    #Assign overall mpf type
    for contig_id in contig_info:
        mpf_type = contig_info[contig_id]['mpf_type'].split(",")
        if len(mpf_type) > 0:
            contig_info[contig_id]['mpf_type'] = determine_mpf_type(mpf_type)

    #blast orit database
    logging.info("Blasting orit sequences {} against {}".format(plasmid_orit, fixed_fasta))
    blastn(input_fasta=plasmid_orit,blastdb=fixed_fasta,min_ident=min_rep_ident,min_cov=min_rep_cov,evalue=min_rep_evalue,min_length=min_length,out_dir=tmp_dir,
           blast_results_file=orit_blast_results,num_threads=num_threads,logging=logging)

    logging.info("Filtering orit blast results {} ".format(orit_blast_results))

    orit_blast_df = BlastReader(orit_blast_results,logging).df
    if len(orit_blast_df) > 0:
        orit_blast_df = recursive_filter_overlap_records(fixStart(orit_blast_df.drop(0).sort_values(['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])), 5, 'sseqid', 'sstart', 'send',
                              'bitscore')
        add_biomarker_results(biomarker_df=orit_blast_df, df_column_name_biomarker='qseqid', df_column_name_seqid='sseqid', contig_info=contig_info,
                          contig_info_type_key='orit_type(s)', contig_info_acs_key='orit_accession(s)', delimeter='|')

    del(orit_blast_df)


    #blast repetitive database
    logging.info("Blasting contigs against repetitive sequences db: {}".format(repetitive_mask_file))
    blastn(input_fasta=fixed_fasta,blastdb=repetitive_mask_file,min_ident=min_rpp_ident,min_cov=min_rpp_cov,evalue=min_rpp_evalue,min_length=min_length,out_dir=tmp_dir,
           blast_results_file=repetitive_blast_results,num_threads=num_threads,logging=logging)
    logging.info("Filtering repetitive blast results {} ".format(repetitive_blast_results))

    repetitive_blast_df = BlastReader(repetitive_blast_results,logging).df
    if len(repetitive_blast_df) > 0:
        repetitive_blast_df = recursive_filter_overlap_records(fixStart(repetitive_blast_df.drop(0).sort_values(['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])), 5, 'qseqid', 'qstart', 'qend',
                                  'bitscore')

        repetitive_list = repetitive_blast_df['qseqid'].tolist()


        #add filtering flag to contigs which are primarially a repetitive element
        for contig_id in repetitive_list:
            if contig_id in contig_info:
                logger.info('Filtering contig: {} due to repetitive sequence'.format(contig_id))
                contig_info[contig_id]['filtering_reason'] = 'repetitve element'
            else:
                logger.error('Contig: {} not found in contig_df this is likely an error'.format(contig_id))


        add_biomarker_results(biomarker_df=repetitive_blast_df, df_column_name_biomarker='sseqid', df_column_name_seqid='qseqid', contig_info=contig_info,
                              contig_info_type_key='repetitive_dna_type', contig_info_acs_key='repetitive_dna_id', delimeter='|',type_col_num=2,type_acs_num=1)


    del(repetitive_blast_df)



    #Filtering contigs against chromosome database

    chrom_filter = False

    if args.genome_filter_db_prefix:
        chrom_filter = True
        genome_filter_db_prefix = args.genome_filter_db_prefix
        logger.info('Genome filter sequences provided: {}'.format(genome_filter_db_prefix))
        matched = (glob.glob(genome_filter_db_prefix+"*"))
        extensions = ['nsq','nin','nhr','nal']
        found = [0,0,0,0]
        for f in matched:
            for i in range(0,len(extensions)):
                e = extensions[i]
                if e in f:
                    found[i]+=1

        for i in found:
            if i == 0:
                logger.error('Error blast database not found with prefix: {}'.format(genome_filter_db_prefix))
                sys.exit()
        if not os.path.isfile(genome_filter_db_prefix + '.msh') :
            logger.error('Error mash sketch not found with prefix: {}'.format(genome_filter_db_prefix))
            sys.exit()


    if chrom_filter:
        cutoff_distance = float(args.mash_genome_neighbor_threshold)
        chr_mash_dists = os.path.join(tmp_dir, 'mash_chr_dists.txt')
        chr_blast_filter = os.path.join(tmp_dir, 'chr_contig_filter_report.txt')
        chr_mash_sketch = genome_filter_db_prefix + ".msh"
        close_genome_reps = find_mash_genomes(chr_mash_sketch , fixed_fasta, chr_mash_dists, cutoff_distance, num_threads)

        if len(close_genome_reps) > 0:
            logger.info('Found close genome matches: {}'.format(",".join(close_genome_reps)))
            seq_id_file = os.path.join(tmp_dir,"seqids.txt")
            sf = open(seq_id_file,'w')
            for s in close_genome_reps:
                sf.write("{}\n".format(s))
            sf.close()

            #fix labels to match the seq id format parsed by makeblastdb
            for i in range(0,len(close_genome_reps)):
                close_genome_reps[i] = "ref|{}|".format(close_genome_reps[i])

            blastn(input_fasta=fixed_fasta, blastdb=genome_filter_db_prefix, min_ident=min_con_ident, min_cov=min_con_cov,
                   evalue=min_con_evalue, min_length=min_length, out_dir=tmp_dir,
                   blast_results_file=chr_blast_filter, num_threads=num_threads, logging=logging,seq_filterfile=seq_id_file)

            chromosome_filter_seqs = BlastReader(chr_blast_filter,logging).df.drop(0)['qseqid'].tolist()

            for contig_id in chromosome_filter_seqs:
                if contig_id in contig_info :
                    if contig_info[contig_id]['filtering_reason'] == 'none':
                        contig_info[contig_id]['filtering_reason'] = 'chromosome'
                        logger.info('Filtering contig: {} due to inclusion in genome filter {}'.format(contig_id,genome_filter_db_prefix))
                else:
                    logger.error('Contig: {} not found in contig_df this is likely an error'.format(contig_id))

            del(chromosome_filter_seqs)

        else:
            logger.info('No close genome matches found')


    # Filter out sequences based on user filter
    if args.filter_db:
        filter_db = args.filter_db
        logger.info('Filter sequences provided: {}'.format(filter_db))
        if not os.path.isfile(filter_db + '.nsq') or \
            os.path.isfile(filter_db + '.nin') or \
            os.path.isfile(filter_db + '.nhr'):
            br = BlastRunner(filter_db,os.path.dirname(filter_db))
            br.makeblastdb(filter_db,dbtype='nucl')

        run_filter = True
    else:
        run_filter = False


    if run_filter:
        logger.info('Blasting input fasta {} against filter db {}'.format(input_fasta, filter_db))
        blast_filter = os.path.join(tmp_dir, 'contig_filter_report.txt')
        blastn(input_fasta=fixed_fasta, blastdb=filter_db, min_ident=min_con_ident, min_cov=min_con_cov,
               evalue=min_con_evalue, min_length=min_length, out_dir=tmp_dir,
               blast_results_file=blast_filter, num_threads=num_threads, logging=logging)

        user_filter_seqs = BlastReader(blast_filter, logging).df
        if len(user_filter_seqs) > 0:
            user_filter_seqs = user_filter_seqs.drop(0)['qseqid'].tolist()
        else:
            user_filter_seqs = []

        for contig_id in user_filter_seqs:
            if contig_id in contig_info:
                contig_info[contig_id]['filtering_reason'] = 'user filter'
                logger.info(
                    'Filtering contig: {} due to inclusion in genome filter {}'.format(contig_id, filter_db))
            else:
                logger.error('Contig: {} not found in contig_df this is likely an error'.format(contig_id))

        del (user_filter_seqs)

    #Identify plasmids likely contained in the sample
    logging.info("Identifying candidate plasmids based on mash screen {}".format(plasmid_ref_db))
    m = mash()
    mash_screen_results = parseMashScreen(m.run_mash_screen(mash_db, fixed_fasta, winner_take_all=True, num_threads=num_threads))
    candidate_plasmids = {}
    logging.info("Filtering plasmid candidates based on distance 0.8")
    for id in mash_screen_results:
        if mash_screen_results[id] > 0.7:
            candidate_plasmids[id] =  mash_screen_results[id]

    sorted(iter(list(candidate_plasmids.items())), key=lambda x: x[1], reverse=True)

    #write seq_filter
    seq_filterfile = os.path.join(tmp_dir,"plasmid_candidates.txt")
    seq_filter = open(seq_filterfile,'w')
    for id in candidate_plasmids:
        seq_filter.write("{}\n".format(id))
    seq_filter.close()


    #blast plasmid database
    logging.info("Blasting contigs against reference sequence db: {}".format(plasmid_ref_db))
    blastn(input_fasta=fixed_fasta,blastdb=plasmid_ref_db,min_ident=min_con_ident,min_cov=min_con_cov,evalue=min_con_evalue,min_length=min_length,out_dir=tmp_dir,
           blast_results_file=contig_blast_results,num_threads=num_threads,logging=logging,seq_filterfile=seq_filterfile)

    logging.info("Filtering contig blast results: {}".format(contig_blast_results))
    contig_blast_df = BlastReader(contig_blast_results,logging).df

    if len(contig_blast_df) > 0:
        contig_blast_df = filter_overlaping_records(fixStart(contig_blast_df.drop(0)), 500, 'qseqid', 'qstart', 'qend',
                                  'bitscore')
        contig_blast_df.reset_index(drop=True)
        #remove blast formatting of seq id
        for index,row in contig_blast_df.iterrows():
            line = row['sseqid'].split('|')
            if len(line) >= 2:
                contig_blast_df.at[index, 'sseqid'] = line[1]

        contig_blast_df = contig_blast_df[contig_blast_df.sseqid.isin(list(reference_sequence_meta.keys()))]

        contig_blast_df.reset_index(drop=True)
        contig_info = assign_contigs_to_clusters(contig_blast_df, reference_sequence_meta, contig_info,tmp_dir,contig_seqs,mash_db,primary_distance,secondary_distance,num_threads)

    results = []
    contig_memberships = {'chromosome':{},'plasmid':{}}
    for contig_id in contig_info:
        data = contig_info[contig_id]

        if data['primary_cluster_id'] != '' and data['mash_nearest_neighbor'] == '':
            data['primary_cluster_id'] = ''
            data['molecule_type'] = 'chromosome'
        if data['filtering_reason'] == 'none':
            data['filtering_reason'] = ''

        if data['molecule_type'] == 'chromosome':
            contig_memberships['chromosome'][contig_id] = ''
        else:
            clust_id = data['primary_cluster_id']
            if not clust_id in contig_memberships['plasmid']:
                contig_memberships['plasmid'][clust_id] = {}
            contig_memberships['plasmid'][clust_id][contig_id] = data

        results.append(data)

    #Write contig report
    logging.info("Writting contig results to {}".format(contig_report))
    if len(results) > 0:
        writeReport(results, MOB_RECON_INFO_HEADER, contig_report)

        logging.info("Writting chromosome sequences to {}".format(chromosome_file))
        chr_fh = open(chromosome_file,'w')
        for contig_id in contig_memberships['chromosome']:
            if contig_id in contig_seqs:
                chr_fh.write(">{}\n{}\n".format(contig_id,contig_seqs[contig_id]))
        chr_fh.close()
        if len(contig_memberships['plasmid']) > 0:

            ncbi = dict_from_alt_key_list(
                read_file_to_dict(NCBI_PLASMID_TAXONOMY_FILE, NCBI_PLASMID_TAXONOMY_HEADER, separater="\t"),
                "sample_id")
            lit = dict_from_alt_key_list(
                read_file_to_dict(LIT_PLASMID_TAXONOMY_FILE, LIT_PLASMID_TAXONOMY_HEADER, separater="\t"), "sample_id")

            build_mobtyper_report(contig_memberships['plasmid'],out_dir,os.path.join(out_dir,"mobtyper_results.txt"),contig_seqs,ncbi,lit)


    if not keep_tmp:
        logging.info("Cleaning up temporary files {}".format(tmp_dir))
        shutil.rmtree(tmp_dir)
    logger.info("Run completed")


# call main function
if __name__ == '__main__':
    main()
