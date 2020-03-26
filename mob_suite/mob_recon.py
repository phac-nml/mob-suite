#!/usr/bin/env python3
from mob_suite.version import __version__
from collections import OrderedDict
import logging, os, shutil, sys, operator,re
from argparse import (ArgumentParser, FileType)
from mob_suite.blast import BlastRunner
from mob_suite.blast import BlastReader
from mob_suite.wrappers import mash
from mob_suite.wrappers import detectCircularity
import glob

from mob_suite.mob_typer import  MOB_TYPER_REPORT_HEADER

from mob_suite.utils import \
    read_fasta_dict, \
    write_fasta_dict, \
    filter_overlaping_records, \
    replicon_blast, \
    mob_blast, \
    repetitive_blast, \
    getRepliconContigs, \
    fix_fasta_header, \
    getMashBestHit, \
    verify_init, \
    check_dependencies, \
    run_mob_typer, \
    read_sequence_info


LOG_FORMAT = '%(asctime)s %(name)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'




def parse_args():
    "Parse the input arguments, use '-h' for help"
    default_database_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'databases')
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
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/ncbi_plasmid_full_seqs.fas'))
    parser.add_argument('--plasmid_mash_db', type=str, required=False,
                        help='Companion Mash database of reference database',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/ncbi_plasmid_full_seqs.fas.msh'))
    parser.add_argument('-m','--plasmid_meta', type=str, required=False,
                        help='MOB-cluster plasmid cluster formatted file matched to the reference plasmid db',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/clusters.txt'))
    parser.add_argument('--plasmid_db_type', type=str, required=False, help='Blast database type of reference database',
                        default='blastn')
    parser.add_argument('--plasmid_replicons', type=str, required=False, help='Fasta of plasmid replicons',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/rep.dna.fas'))
    parser.add_argument('--repetitive_mask', type=str, required=False, help='Fasta of known repetitive elements',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/repetitive.dna.fas'))
    parser.add_argument('--plasmid_mob', type=str, required=False, help='Fasta of plasmid relaxases',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/mob.proteins.faa'))
    parser.add_argument('--plasmid_mpf', type=str, required=False, help='Fasta of known plasmid mate-pair proteins',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/mpf.proteins.faa'))
    parser.add_argument('--plasmid_orit', type=str, required=False, help='Fasta of known plasmid oriT dna sequences',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/orit.fas'))
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


def init_console_logger(lvl=2):
    root = logging.getLogger()

    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]

    report_lvl = logging_levels[lvl]
    root.setLevel(report_lvl)  # set root logger level

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)

    return logging.getLogger(__name__)


def find_mash_genomes(reference_mash_sketch,fasta_query,outfile,cutoff_distance,num_threads=1):
    cutoff_distance = float(cutoff_distance)
    if (not os.path.isfile(fasta_query)):
        logger.error('Error query fasta file missing "{}"'.format(fasta_query))
        sys.exit(-1)

    if (not os.path.isfile(reference_mash_sketch)):
        logger.error('Error genome mash sketch file missing "{}"'.format(reference_mash_sketch))
        sys.exit(-1)

    output_fh = open(outfile,'w',encoding="utf-8")
    m = mash()
    m.run_mash(reference_mash_sketch,fasta_query,output_fh,num_threads=num_threads)
    output_fh.close()
    mash_results = m.read_mash(outfile)

    genomes = []

    for line in mash_results:
        row = line.strip("\n").split("\t")
        seqid = row[0]
        score = float(row[2])
        matches = row[4]
        if score < cutoff_distance:
            genomes.append(seqid)

    return genomes


def blast_closed_genomes(input_fasta,filter_db,out_dir,blast_filter,min_con_ident, min_con_cov, min_con_evalue, include_list=[], seq_id_file=None, num_threads=1):
    blast_runner = BlastRunner(input_fasta, out_dir)
    blast_runner.run_blast(query_fasta_path=input_fasta, blast_task='megablast', db_path=filter_db,
                           db_type='nucl', min_cov=min_con_cov, min_ident=min_con_ident, evalue=min_con_evalue,
                           blast_outfile=blast_filter, num_threads=num_threads, word_size=11,seq_id_file=seq_id_file )


    if os.path.getsize(blast_filter) == 0:
        return []

    blast_df = BlastReader(blast_filter).df
    blast_df = blast_df.loc[blast_df['pident'] >= min_con_ident]
    blast_df = blast_df.loc[blast_df['evalue'] <= min_con_evalue]
    blast_df = blast_df.loc[blast_df['qcovs'] >= min_con_cov]
    blast_df = blast_df[blast_df['sseqid'].isin(include_list)]
    blast_df = blast_df.reset_index(drop=True)
    blast_df.to_csv(blast_filter, sep='\t', header=True, line_terminator='\n', index=False)

    return list(set(blast_df['qseqid'].tolist()))






def contig_blast(input_fasta, plasmid_db, min_ident, min_cov, evalue, min_length, tmp_dir, blast_results_file,
                 num_threads=1):
    blast_runner = None
    filtered_blast = os.path.join(tmp_dir, 'filtered_blast.txt')
    blast_runner = BlastRunner(input_fasta, tmp_dir)
    blast_runner.run_blast(query_fasta_path=input_fasta, blast_task='megablast', db_path=plasmid_db,
                           db_type='nucl', min_cov=min_cov, min_ident=min_ident, evalue=evalue,
                           blast_outfile=blast_results_file, num_threads=num_threads, word_size=11)
    if os.path.getsize(blast_results_file) == 0:
        fh = open(filtered_blast, 'w', encoding="utf-8")
        fh.write('')
        fh.close()
        return dict()
    blast_df = BlastReader(blast_results_file).df
    blast_df = blast_df.loc[blast_df['length'] >= min_length]
    blast_df = blast_df.loc[blast_df['qlen'] <= 400000]
    blast_df = blast_df.loc[blast_df['qlen'] >= min_length]
    blast_df = blast_df.loc[blast_df['qcovs'] >= min_cov]
    blast_df = blast_df.reset_index(drop=True)
    blast_df.to_csv(filtered_blast, sep='\t', header=False, line_terminator='\n', index=False)


def membership_voting(reference_sequence_hits,contig_hit_scores):
    filtered_reference_hits = {}
    filtered_cluster_ids = {}

    for ref_id in reference_sequence_hits:
        if float(reference_sequence_hits[ref_id]['covered_bases'])/reference_sequence_hits[ref_id]['length'] < 0.8:
            continue
        filtered_reference_hits[ref_id] = reference_sequence_hits[ref_id]
        filtered_cluster_ids[reference_sequence_hits[ref_id]['clust_id']] = reference_sequence_hits[ref_id]['score']
    memberships = {}

    for contig_id in contig_hit_scores:
        top_hit_score = 0
        top_hit_cluster_ids = []
        top_hit_ids = []

        for hit_id in contig_hit_scores[contig_id]:
            score = contig_hit_scores[contig_id][hit_id]['score']
            if score >= top_hit_score:
                top_hit_ids.append(hit_id)
                top_hit_cluster_ids.append(reference_sequence_hits[hit_id]['clust_id'])
                top_hit_score = score


        for hit_id in top_hit_ids:
            if hit_id in filtered_reference_hits:
                memberships[contig_id] = {'ref_id': hit_id, 'clust_id': reference_sequence_hits[hit_id]['clust_id'],'score':reference_sequence_hits[hit_id]['score']}
                break

        if contig_id not in memberships:
            for i in range(0,len(top_hit_cluster_ids)):
                clust_id = top_hit_cluster_ids[i]
                hit_id = top_hit_ids[i]
                if clust_id in filtered_cluster_ids:
                    memberships[contig_id] = {'ref_id': hit_id, 'clust_id': reference_sequence_hits[hit_id]['clust_id'],'score':reference_sequence_hits[hit_id]['score']}
                    break


    return memberships

def calc_hit_coverage(blast_df,overlap_threshold,reference_sequence_meta):
    blast_df = blast_df.sort_values(['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])
    hit_scores = {}
    size = str(len(blast_df))
    prev_size = 0
    while size != prev_size:
        blast_df = filter_overlaping_records(blast_df, overlap_threshold, 'qseqid', 'qstart', 'qend', 'bitscore')
        prev_size = size
        size = str(len(blast_df))

    for index, row in blast_df.iterrows():
        query = row['qseqid']
        pID = row['sseqid']
        score = row['bitscore']

        if pID not in reference_sequence_meta:
            print("->{}\t{}".format(query,pID))
            continue
        else:
            clust_id = reference_sequence_meta[pID]['primary_cluster_id']


        if query not in hit_scores:
            hit_scores[query] = {}

        if clust_id not in hit_scores[query]:
            hit_scores[query][clust_id] = 0

        if hit_scores[query][clust_id] < score:
            hit_scores[query][clust_id] = score

    cluster_scores = {}

    for q in hit_scores:
        for clust_id in hit_scores[q]:
            if clust_id not in cluster_scores:
                cluster_scores[clust_id] = 0
            cluster_scores[clust_id] += hit_scores[q][clust_id]

    sorted_cluster_scores = OrderedDict(sorted(iter(list(cluster_scores.items())), key=lambda x: x[1], reverse=True))

    return sorted_cluster_scores






def contig_blast_group(blast_results_file, overlap_threshold,reference_sequence_meta):
    if os.path.getsize(blast_results_file) == 0:
        return ({},{},{})
    blast_df = BlastReader(blast_results_file).df
    blast_df = blast_df.sort_values(['qseqid', 'sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, True, False])


    blast_df = filter_overlaping_records(blast_df, overlap_threshold, 'sseqid', 'sstart', 'send', 'bitscore')
    blast_df = blast_df.sort_values(['qseqid', 'sseqid', 'qstart', 'qend', 'bitscore'], ascending=[True, True, True, True, False])
    blast_df = filter_overlaping_records(blast_df, overlap_threshold, 'sseqid', 'qstart', 'qend', 'bitscore')

    size = str(len(blast_df))
    prev_size = 0
    while size != prev_size:
        blast_df = filter_overlaping_records(blast_df, overlap_threshold, 'sseqid', 'sstart', 'send', 'bitscore')
        prev_size = size
        size = str(len(blast_df))

    sorted_cluster_scores = calc_hit_coverage(blast_df, overlap_threshold,reference_sequence_meta)


    cluster_scores = dict()
    groups = dict()
    hits = dict()
    contigs = dict()
    query_hit_scores = dict()


    for index, row in blast_df.iterrows():
        query = row['qseqid']
        pID = row['sseqid']
        if pID not in reference_sequence_meta:
            continue
        else:
            clust_id = reference_sequence_meta[pID]['primary_cluster_id']

        score = row['bitscore']
        pLen = row['slen']
        contig_id = row['qseqid']
        mlength = row['length']

        if not query in query_hit_scores:
            query_hit_scores[query] = {}

        if not pID in query_hit_scores[query]:
            query_hit_scores[query][pID] = {'score': 0, 'length': pLen, 'covered_bases': 0, 'clust_id': clust_id}

        query_hit_scores[query][pID]['score']+= score
        query_hit_scores[query][pID]['covered_bases'] += mlength

        if not pID in hits:
            hits[pID] = {'score': 0, 'length': pLen, 'covered_bases': 0, 'clust_id': clust_id}

        if not clust_id in cluster_scores:
            cluster_scores[clust_id] = score

        elif score > cluster_scores[clust_id]:
            cluster_scores[clust_id] = score

        if not clust_id in groups:
            groups[clust_id] = dict()

        if not query in groups[clust_id]:
            groups[clust_id][query] = dict()

        if not contig_id in contigs:
            contigs[contig_id] = dict()

        if not clust_id in contigs[contig_id]:
            contigs[contig_id][clust_id] = 0

        if contigs[contig_id][clust_id] < score:
            contigs[contig_id][clust_id] = score

        groups[clust_id][query][contig_id] = score

        hits[pID]['score'] += score
        hits[pID]['covered_bases'] += mlength



    #sorted_cluster_scores = OrderedDict(sorted(iter(list(cluster_scores.items())), key=lambda x: x[1], reverse=True))

    unassigned_contigs = contigs
    print(len(unassigned_contigs))
    while len(unassigned_contigs) > 0 and len(sorted_cluster_scores) > 0:
        clust_id = next(iter(sorted_cluster_scores))
        (unassigned_contigs, assigned_contigs, sorted_cluster_scores) = assign_contigs(clust_id, contigs, query_hit_scores, sorted_cluster_scores)
        for contig_id in assigned_contigs:
            contigs[contig_id] = assigned_contigs[contig_id]
        print(len(unassigned_contigs))
        print(sorted_cluster_scores)

        remove_contig_list = []

        for contig_id in unassigned_contigs:
            clust_found = False
            for contig_clust_id in unassigned_contigs[contig_id]:
                if contig_clust_id in sorted_cluster_scores:
                    clust_found = True
            if not clust_found:
                remove_contig_list.append(contig_id)

        for contig_id in remove_contig_list:
            if contig_id in unassigned_contigs:
                del(unassigned_contigs[contig_id])



    return (contigs,hits,query_hit_scores)


def assign_contigs(clust_id,contigs,query_hit_scores,cluster_scores):
    print(cluster_scores)
    unassigned_contigs = {}
    assigned_contigs = {}
    for contig_id in contigs:
        if len(contigs[contig_id]) <= 1:
            print(contigs[contig_id])
            continue
        print("---->{}\t{}".format(contig_id,clust_id))
        print(contigs[contig_id])
        if clust_id in contigs[contig_id]:
            assigned_contigs[contig_id] = {clust_id: contigs[contig_id][clust_id]}
            for hit_id in query_hit_scores[contig_id]:
                score = query_hit_scores[contig_id][hit_id]['score']
                hit_clust_id = query_hit_scores[contig_id][hit_id]['clust_id']

                if hit_clust_id in cluster_scores:
                    print("==>Found")
                    print(cluster_scores[hit_clust_id])
                    cluster_scores[hit_clust_id] -= score
                    print(cluster_scores[hit_clust_id])
                    if cluster_scores[hit_clust_id] <= 0:
                        del(cluster_scores[hit_clust_id])

        else:
            print("{}\t{}-unasigned".format(clust_id,contig_id))
            unassigned_contigs[contig_id] = contigs[contig_id]
    if(clust_id in cluster_scores):
        del(cluster_scores[clust_id])
    cluster_scores = OrderedDict(sorted(iter(list(cluster_scores.items())), key=lambda x: x[1], reverse=True))

    return(unassigned_contigs,assigned_contigs,cluster_scores)





def circularize(input_fasta, outdir):
    c = detectCircularity()
    return c.run(input_fasta, outdir)

def filter_user_sequences(input_fasta,filter_db,out_dir,blast_filter,min_con_ident, min_con_cov, min_con_evalue, num_threads=1):

    blast_runner = BlastRunner(input_fasta, out_dir)
    blast_runner.run_blast(query_fasta_path=input_fasta, blast_task='megablast', db_path=filter_db,
                           db_type='nucl', min_cov=min_con_cov, min_ident=min_con_ident, evalue=min_con_evalue,
                           blast_outfile=blast_filter, num_threads=num_threads, word_size=11)
    if os.path.getsize(blast_filter) == 0:
        return []

    blast_df = BlastReader(blast_filter).df
    blast_df = blast_df.loc[blast_df['pident'] >= min_con_ident]
    blast_df = blast_df.loc[blast_df['evalue'] <= min_con_evalue]
    blast_df = blast_df.loc[blast_df['qcovs'] >= min_con_cov]
    blast_df = blast_df.reset_index(drop=True)
    blast_df.to_csv(blast_filter, sep='\t', header=False, line_terminator='\n', index=False)

    return list(set(blast_df['qseqid'].tolist()))




def main():

    args = parse_args()


    if args.debug:
        logger = init_console_logger(3)
    else:
        logger = init_console_logger(2)

    logger.info("MOB-recon version {} ".format(__version__))
    logger.debug("Debug log reporting set on successfully")


    if not args.outdir:
        logger.error('Error, no output directory specified, please specify one')
        sys.exit(-1)

    if not args.infile:
        logger.error('Error, no fasta specified, please specify one')
        sys.exit(-1)

    if not os.path.isfile(args.infile):
        logger.error('Error, input fasta file does not exist: "{}"'.format(args.infile))
        sys.exit(-1)

    logger.info('Processing fasta file {}'.format(args.infile))
    logger.info('Analysis directory {}'.format(args.outdir))

    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir, 0o755)
    elif not args.force:
        logger.error("Error output directory exisits, please specify a new directory or use --force to overwrite")
        sys.exit(-1)
    else:
        shutil.rmtree(args.outdir)
        os.mkdir(args.outdir, 0o755)

    max_contig_size = args.max_contig_size
    max_plasmid_size = args.max_plasmid_size

    # Check that the needed databases have been initialized
    database_dir = os.path.abspath(args.database_directory)
    verify_init(logger, database_dir)

    if args.sample_id is None:
        sample_id = re.sub(r"\.(fasta|fas|fa){1,1}", "", os.path.basename(args.infile))
    else:
        sample_id = args.sample_id



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

    plasmid_files = []
    input_fasta = args.infile
    out_dir = args.outdir
    num_threads = args.num_threads
    tmp_dir = os.path.join(out_dir, '__tmp')
    file_id = os.path.basename(input_fasta)
    fixed_fasta = os.path.join(tmp_dir, 'fixed.input.fasta')
    chromosome_file = os.path.join(out_dir, 'chromosome.fasta')
    replicon_blast_results = os.path.join(tmp_dir, 'replicon_blast_results.txt')
    mob_blast_results = os.path.join(tmp_dir, 'mobrecon_blast_results.txt')
    repetitive_blast_results = os.path.join(tmp_dir, 'repetitive_blast_results.txt')
    contig_blast_results = os.path.join(tmp_dir, 'contig_blast_results.txt')



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

    min_overlapp = int(args.min_overlap)
    min_length = int(args.min_length)



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


    # Input Databases
    default_database_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'databases')


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
        mpf_ref = args.os.path.join(database_dir, 'mpf.proteins.faa')
        plasmid_orit = args.os.path.join(database_dir, 'orit.fas')

    reference_sequence_meta = read_sequence_info(plasmid_meta)


    check_dependencies(logger)

    needed_dbs = [plasmid_ref_db, replicon_ref, mob_ref, mash_db, repetitive_mask_file,"{}.nin".format(repetitive_mask_file)]

    for db in needed_dbs:
        if (not os.path.isfile(db)):
            logger.error('Error needed database missing "{}"'.format(db))
            sys.exit(-1)

    contig_report_file = os.path.join(out_dir, 'contig_report.txt')
    filtered_blast = os.path.join(tmp_dir, 'filtered_blast.txt')
    repetitive_blast_report = os.path.join(out_dir, 'repetitive_blast_report.txt')
    mobtyper_results_file = os.path.join(out_dir, 'mobtyper_aggregate_report.txt')
    keep_tmp = args.keep_tmp

    run_overhang = args.run_overhang
    unicycler_contigs = args.unicycler_contigs

    if not isinstance(args.num_threads, int):
        logger.info('Error number of threads must be an integer, you specified "{}"'.format(args.num_threads))

    logger.info('Creating tmp working directory {}'.format(tmp_dir))
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir, 0o755)

    #Filtering against chromosome database

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

    logger.info('Writing cleaned header input fasta file from {} to {}'.format(input_fasta, fixed_fasta))
    fix_fasta_header(input_fasta, fixed_fasta)
    contig_seqs = read_fasta_dict(fixed_fasta)


    seq_filter = []

    if chrom_filter:
        cutoff_distance = float(args.mash_genome_neighbor_threshold)
        chr_blast_filter = os.path.join(out_dir, 'contig_filter_report.txt')
        chr_mash_sketch = genome_filter_db_prefix + ".msh"
        close_genome_reps = find_mash_genomes(chr_mash_sketch , fixed_fasta, chr_blast_filter, cutoff_distance, num_threads=1)
        logger.info('Found close genome matches: {}'.format(",".join(close_genome_reps)))
        if len(close_genome_reps) > 0:
            seq_id_file = os.path.join(tmp_dir,"seqids.txt")
            sf = open(seq_id_file,'w')
            for s in close_genome_reps:
                sf.write("{}\n".format(s))
            sf.close()

            #fix labels to match the seq id format parsed by makeblastdb
            for i in range(0,len(close_genome_reps)):
                close_genome_reps[i] = "ref|{}|".format(close_genome_reps[i])

            seq_filter = blast_closed_genomes(fixed_fasta, genome_filter_db_prefix, out_dir, chr_blast_filter, min_con_ident, min_con_cov,
                                 min_con_evalue, close_genome_reps,seq_id_file, num_threads)



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




    #Filter out sequences based on user filter
    if run_filter:
        logger.info('Blasting input fasta {} against filter db {}'.format(input_fasta, filter_db))
        blast_filter = os.path.join(out_dir, 'contig_filter_report.txt')
        seq_filter = seq_filter  + filter_user_sequences(fixed_fasta,filter_db,out_dir,blast_filter,min_con_ident, min_con_cov, min_con_evalue, num_threads)

    if run_filter or chrom_filter:
        filtered_seqs = {}
        for seq_id in contig_seqs:
            if seq_id not in seq_filter:
                filtered_seqs[seq_id] = contig_seqs[seq_id]

        write_fasta_dict(filtered_seqs,fixed_fasta)
        del(filtered_seqs)

        if os.path.getsize(fixed_fasta) == 0:
            logger.error('Error no sequences are present in fasta file {} after filtering against filter db {}'.format(input_fasta, filter_db))
            sys.exit()


    logger.info('Running replicon blast on {}'.format(replicon_ref))
    replicon_contigs = getRepliconContigs(
        replicon_blast(replicon_ref, fixed_fasta, min_rep_ident, min_rep_cov, min_rep_evalue, tmp_dir,
                       replicon_blast_results,
                       num_threads=num_threads))

    logger.info('Running relaxase blast on {}'.format(mob_ref))
    mob_contigs = getRepliconContigs(
        mob_blast(mob_ref, fixed_fasta, min_mob_ident, min_mob_cov, min_mob_evalue, tmp_dir, mob_blast_results,
                  num_threads=num_threads))

    logger.info('Running contig blast on {}'.format(plasmid_ref_db))
    contig_blast(fixed_fasta, plasmid_ref_db, min_con_ident, min_con_cov, min_con_evalue, min_length,
                 tmp_dir, contig_blast_results)

    (pcl_clusters,reference_sequence_hits,contig_hit_scores) = contig_blast_group(filtered_blast, min_overlapp,reference_sequence_meta)
    contig_cluster_membership = membership_voting(reference_sequence_hits,contig_hit_scores)

    for contig_id in contig_cluster_membership:
        if contig_id in pcl_clusters:
            pcl_clusters[contig_id] = {contig_cluster_membership[contig_id]['clust_id']:contig_cluster_membership[contig_id]['score']}


    logger.info('Running repetitive contig masking blast on {}'.format(mob_ref))
    repetitive_contigs = repetitive_blast(fixed_fasta, repetitive_mask_file, min_rpp_ident, min_rpp_cov, min_rpp_evalue,
                                          min_length, tmp_dir,
                                          repetitive_blast_results, num_threads=num_threads)

    circular_contigs = dict()


    if run_overhang:
        logger.info('Running internal circular contig detection on {}'.format(fixed_fasta))
        circular_contigs = circularize(fixed_fasta, tmp_dir)

    if unicycler_contigs:
        for seqid in contig_seqs:
            if 'circular=true' in seqid:
                circular_contigs[seqid] = ''

    repetitive_dna = dict()
    results_fh = open(repetitive_blast_report, 'w', encoding="utf-8")
    results_fh.write("contig_id\tmatch_id\tmatch_type\tscore\tcontig_match_start\tcontig_match_end\n")

    for contig_id in repetitive_contigs:
        match_info = repetitive_contigs[contig_id]['id'].split('|')
        repetitive_dna[contig_id] = "{}\t{}\t{}\t{}\t{}".format(
            match_info[1],
            match_info[len(match_info) - 1],
            repetitive_contigs[contig_id]['score'],
            repetitive_contigs[contig_id]['contig_start'],
            repetitive_contigs[contig_id]['contig_end'])
        results_fh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(contig_id,
                                                           match_info[1],
                                                           match_info[len(match_info) - 1],
                                                           repetitive_contigs[contig_id]['score'],
                                                           repetitive_contigs[contig_id]['contig_start'],
                                                           repetitive_contigs[contig_id]['contig_end']))

    results_fh.close()

    seq_clusters = dict()
    cluster_bitscores = dict()
    for seqid in pcl_clusters:
        cluster_id = list(pcl_clusters[seqid].keys())[0]
        bitscore = pcl_clusters[seqid][cluster_id]
        cluster_bitscores[cluster_id] = bitscore

    sorted_cluster_bitscores = sorted(list(cluster_bitscores.items()), key=operator.itemgetter(1))
    sorted_cluster_bitscores.reverse()
    contigs_assigned = dict()
    for cluster_id, bitscore in sorted_cluster_bitscores:

        if not cluster_id in seq_clusters:
            seq_clusters[cluster_id] = dict()
        for seqid in pcl_clusters:
            if not cluster_id in pcl_clusters[seqid]:
                continue
            if seqid in contig_seqs and seqid not in contigs_assigned:
                seq_clusters[cluster_id][seqid] = contig_seqs[seqid]
                contigs_assigned[seqid] = cluster_id

    # Add sequences with known replicons regardless of whether they belong to a mcl cluster
    clust_id = 0
    refined_clusters = dict()
    for contig_id in mob_contigs:
        if not contig_id in pcl_clusters:
            if contig_id in contig_seqs:
                if not clust_id in seq_clusters:
                    seq_clusters["Novel_" + str(clust_id)] = dict()
                    if not contig_id in pcl_clusters:
                        pcl_clusters[contig_id] = dict()

                    pcl_clusters[contig_id]["Novel_" + str(clust_id)] = 0
                seq_clusters["Novel_" + str(clust_id)][contig_id] = contig_seqs[contig_id]
            clust_id += 1

    # Add sequences with known relaxases regardless of whether they belong to a mcl cluster

    for contig_id in replicon_contigs:
        if not contig_id in pcl_clusters:
            if contig_id in contig_seqs:
                if not clust_id in seq_clusters:
                    seq_clusters["Novel_" + str(clust_id)] = dict()
                    if not contig_id in pcl_clusters:
                        pcl_clusters[contig_id] = dict()

                    pcl_clusters[contig_id]["Novel_" + str(clust_id)] = dict()
                seq_clusters["Novel_" + str(clust_id)][contig_id] = contig_seqs[contig_id]
            clust_id += 1

    refined_clusters = dict()

    # split out circular sequences from each other

    replicon_clusters = dict()
    for contig_id in replicon_contigs:

        for hit_id in replicon_contigs[contig_id]:
            id, rep_type = hit_id.split('|')

            cluster = list(pcl_clusters[contig_id].keys())[0]
            if not cluster in replicon_clusters:
                replicon_clusters[cluster] = 0
            replicon_clusters[cluster] += 1

    for id in seq_clusters:
        cluster = seq_clusters[id]

        if not id in refined_clusters:
            refined_clusters[id] = dict()

        for contig_id in cluster:
            if contig_id in circular_contigs and len(cluster) > 1 and (
                    id in replicon_clusters and replicon_clusters[id] > 1):
                if not clust_id in refined_clusters:
                    refined_clusters["Novel_" + str(clust_id)] = dict()
                refined_clusters["Novel_" + str(clust_id)][contig_id] = cluster[contig_id]
                clust_id += 1
                continue

            refined_clusters[id][contig_id] = cluster[contig_id]


    seq_clusters = refined_clusters
    m = mash()

    results_fh = open(contig_report_file, 'w', encoding="utf-8")
    results_fh.write("sample_id\tprimary_cluster_id\tsecondary_cluster_id\tcontig_id\tcontig_length\tcircularity_status\trep_type\t" \
                     "rep_type_accession\trelaxase_type\trelaxase_type_accession\tmash_nearest_neighbor\t"
                     " mash_neighbor_distance\trepetitive_dna_id\tmatch_type\tscore\tcontig_match_start\tcontig_match_end\n")

    filter_list = dict()
    counter = 0

    for cluster in seq_clusters:
        clusters = seq_clusters[cluster]
        total_cluster_length = 0

        count_seqs = len(clusters)
        count_rep = 0
        count_small = 0
        temp = dict()

        for contig_id in clusters:

            if contig_id in repetitive_contigs:
                count_rep += 1
            length = len(clusters[contig_id])
            total_cluster_length += length
            if length > max_contig_size:
                continue
            if length < 3000:
                count_small += 1
            temp[contig_id] = ''

        if count_rep == count_seqs or (float(
                count_rep) / count_seqs * 100 > 50 and
                        count_small == count_seqs) or \
                        total_cluster_length < 1400 or \
                        total_cluster_length > max_plasmid_size:
            continue

        for contig_id in temp:
            filter_list[contig_id] = ''

        cluster_file = os.path.join(tmp_dir, 'clust_' + str(cluster) + '.fasta')
        mash_file = os.path.join(tmp_dir, 'clust_' + str(cluster) + '.txt')
        write_fasta_dict(clusters, cluster_file)

        mashfile_handle = open(mash_file, 'w',encoding="utf-8")
        m.run_mash(mash_db, cluster_file, mashfile_handle)

        mash_results = m.read_mash(mash_file)
        mash_top_hit = getMashBestHit(mash_results)

        plasmid_primary_acs = "{}_novel".format(counter)
        plasmid_secondary_acs = "{}_novel".format(counter)

        if mash_top_hit['top_hit'] in reference_sequence_meta:
            if mash_top_hit['mash_hit_score'] <= primary_distance:
                plasmid_primary_acs = reference_sequence_meta[mash_top_hit['top_hit']]['primary_cluster_id']
            if mash_top_hit['mash_hit_score'] <= secondary_distance:
                plasmid_secondary_acs = reference_sequence_meta[mash_top_hit['top_hit']]['secondary_cluster_id']

        # delete low scoring clusters
        if float(mash_top_hit['mash_hit_score']) > primary_distance:
            skip = True
            for contig_id in clusters:
                if contig_id in replicon_contigs:
                    skip = False
                    break
                if contig_id in circular_contigs:
                    skip = False
                    break
                if contig_id in mob_contigs:
                    skip = False
                    break
            if skip:
                for contig_id in clusters:
                    del (filter_list[contig_id])
                continue

        new_clust_file = None
        if os.path.isfile(cluster_file):
            if float(mash_top_hit['mash_hit_score']) < primary_distance:
                new_clust_file = os.path.join(out_dir, "{}_plasmid_{}.fasta".format(sample_id, plasmid_primary_acs))

            else:
                new_clust_file = os.path.join(out_dir, "{}_plasmid_{}.fasta".format(sample_id, plasmid_primary_acs))
                counter += 1

            if os.path.isfile(new_clust_file):
                temp_fh = open(cluster_file, 'r', encoding="utf-8")

                data = temp_fh.read()

                temp_fh.close()
                temp_fh = open(new_clust_file, 'a', encoding="utf-8")
                temp_fh.write(data)
                temp_fh.close()
                mash_file = os.path.join(tmp_dir, "clust_{}.txt".format(plasmid_primary_acs))
                mashfile_handle = open(mash_file, 'w', encoding="utf-8")
                m.run_mash(mash_db, cluster_file, mashfile_handle)
                mash_results = m.read_mash(mash_file)
                mash_top_hit = getMashBestHit(mash_results)

            else:
                os.rename(cluster_file, new_clust_file)

        if new_clust_file is not None:
            plasmid_files.append(new_clust_file)


        for contig_id in clusters:
            found_replicon_string = ''
            found_replicon_id_string = ''
            found_mob_string = ''
            found_mob_id_string = ''
            if run_overhang or unicycler_contigs:
                contig_status = 'Incomplete'
            else:
                contig_status = 'Not Tested'

            if contig_id in circular_contigs:
                contig_status = 'Circular'

            if contig_id in replicon_contigs:
                rep_ids = dict()
                rep_hit_ids = dict()

                for hit_id in replicon_contigs[contig_id]:
                    id, rep_type = hit_id.split('|')
                    rep_ids[rep_type] = ''
                    rep_hit_ids[id] = ''

                found_replicon_string = ','.join(list(rep_ids.keys()))
                found_replicon_id_string = ','.join(list(rep_hit_ids.keys()))

            if contig_id in mob_contigs:
                mob_ids = dict()
                mob_hit_ids = dict()

                for hit_id in mob_contigs[contig_id]:
                    id, mob_type = hit_id.split('|')
                    mob_ids[mob_type] = ''
                    mob_hit_ids[id] = ''

                found_mob_string = ','.join(list(mob_ids.keys()))
                found_mob_id_string = ','.join(list(mob_hit_ids.keys()))

            rep_dna_info = "\t\t\t\t"
            if contig_id in repetitive_dna:
                rep_dna_info = repetitive_dna[contig_id]

            results_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample_id,
                                                                                       plasmid_primary_acs,plasmid_secondary_acs, contig_id,
                                                                                       len(clusters[contig_id]),
                                                                                       contig_status,
                                                                                       found_replicon_string,
                                                                                       found_replicon_id_string,
                                                                                       found_mob_string,
                                                                                       found_mob_id_string,
                                                                                       mash_top_hit['top_hit'],
                                                                                       mash_top_hit['mash_hit_score'],
                                                                                       rep_dna_info))
    chr_contigs = dict()

    for contig_id in contig_seqs:
        if contig_id not in filter_list:
            chr_contigs[contig_id] = contig_seqs[contig_id]
            rep_dna_info = "\t\t\t\t"
            if contig_id in repetitive_dna:
                rep_dna_info = repetitive_dna[contig_id]

            found_replicon_string = ''
            found_replicon_id_string = ''
            found_mob_string = ''
            found_mob_id_string = ''

            if run_overhang or unicycler_contigs:
                contig_status = 'Incomplete'
            else:
                contig_status = 'Not Tested'
            if contig_id in circular_contigs:
                contig_status = 'Circular'

            if contig_id in replicon_contigs:
                rep_ids = dict()
                rep_hit_ids = dict()

                for hit_id in replicon_contigs[contig_id]:
                    id, rep_type = hit_id.split('|')
                    rep_ids[rep_type] = ''
                    rep_hit_ids[id] = ''

                found_replicon_string = ','.join(list(rep_ids.keys()))
                found_replicon_id_string = ','.join(list(rep_hit_ids.keys()))

            if contig_id in mob_contigs:
                mob_ids = dict()
                mob_hit_ids = dict()

                for hit_id in mob_contigs[contig_id]:
                    id, mob_type = hit_id.split('|')
                    mob_ids[mob_type] = ''
                    mob_hit_ids[id] = ''

                found_mob_string = ','.join(list(mob_ids.keys()))
                found_mob_id_string = ','.join(list(mob_hit_ids.keys()))

            results_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample_id,
                                                                                       'chromosome',
                                                                                       '-', contig_id,
                                                                                       len(contig_seqs[contig_id]),
                                                                                       contig_status,
                                                                                       found_replicon_string,
                                                                                       found_replicon_id_string,
                                                                                       found_mob_string,
                                                                                       found_mob_id_string,
                                                                                       '-',
                                                                                       '-',
                                                                                       rep_dna_info))
    results_fh.close()
    write_fasta_dict(chr_contigs, chromosome_file)

    if args.run_typer:
        mobtyper_results = "{}\n".format("\t".join(MOB_TYPER_REPORT_HEADER))
        for plasmid_file_abs_path in plasmid_files:
            mob_out_file = os.path.join(out_dir,"mob_typer_{}_report.txt".format(os.path.splitext(os.path.basename(plasmid_file_abs_path) )[0]))

            mobtyper_results = mobtyper_results + "{}\n".format(
                run_mob_typer(plasmid_file_abs_path=plasmid_file_abs_path,
                          out_file=mob_out_file,
                          mash_db=mash_db,
                          replicon_ref=replicon_ref,
                          plasmid_meta=plasmid_meta,
                          mob_ref=mob_ref,
                          mpf_ref=mpf_ref,
                          plasmid_orit=plasmid_orit,
                          num_threads=str(num_threads)))

        fh = open(mobtyper_results_file, 'w', encoding="utf-8")
        fh.write(mobtyper_results)
        fh.close()

    if not keep_tmp:
        shutil.rmtree(tmp_dir)


# call main function
if __name__ == '__main__':
    main()

#TODO
#Resolve empty aggregated reports bug
#Plasmid size inflation with the same contigs