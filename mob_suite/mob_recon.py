#!/usr/bin/env python3
from mob_suite.version import __version__
from collections import OrderedDict
import logging, os, shutil, sys, re
import pandas as pd
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter,RawDescriptionHelpFormatter)
from mob_suite.blast import BlastRunner
from mob_suite.blast import BlastReader
from mob_suite.wrappers import mash
from mob_suite.wrappers import detectCircularity

import glob

from mob_suite.constants import \
    MOB_CLUSTER_INFO_HEADER, \
    MOB_RECON_INFO_HEADER, \
    default_database_dir, \
    LOG_FORMAT, \
    LIT_PLASMID_TAXONOMY_HEADER

from mob_suite.utils import \
    read_fasta_dict, \
    write_fasta_dict, \
    filter_overlaping_records, \
    fix_fasta_header, \
    verify_init, \
    check_dependencies, \
    read_sequence_info, \
    fixStart, \
    calc_md5, \
    gc_fraction, \
    ETE3_db_status_check, \
    writeReport, \
    dict_from_alt_key_list, \
    read_file_to_dict, \
    blastn, \
    identify_biomarkers, \
    build_mobtyper_report, \
    parseMash, \
    blast_mge, \
    writeMGEresults, \
    create_biomarker_dataframe


def parse_args():
    "Parse the input arguments, use '-h' for help"
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass
    parser = ArgumentParser(
        description="MOB-Recon: Typing and reconstruction of plasmids from draft and complete assemblies version: {}".format(
            __version__), formatter_class=CustomFormatter)
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output Directory to put results')
    parser.add_argument('-i', '--infile', type=str, required=True, help='Input assembly fasta file to process')
    parser.add_argument('-n', '--num_threads', type=int, required=False, help='Number of threads to be used', default=1)
    parser.add_argument('-s', '--sample_id', type=str, required=False, help='Sample Prefix for reports')
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')
    parser.add_argument('-b', '--filter_db', type=str, required=False, help='Path to fasta file to mask sequences')
    parser.add_argument('-g', '--genome_filter_db_prefix', type=str, required=False,
                        help='Prefix of mash sketch and blastdb of closed chromosomes to use for auto selection of close genomes for filtering')
    parser.add_argument('-p', '--prefix', type=str, required=False, help='Prefix to append to result files')
    parser.add_argument('--mash_genome_neighbor_threshold', type=float, required=False,
                        help='Mash distance selecting valid closed genomes to filter', default=0.002)

    parser.add_argument('--max_contig_size', type=int, required=False,
                        help='Maximum size of a contig to be considered a plasmid',
                        default=450000)
    parser.add_argument('--max_plasmid_size', type=int, required=False,
                        help='Maximum size of a reconstructed plasmid',
                        default=450000)
    parser.add_argument('--min_rep_evalue', type=float, required=False,
                        help='Minimum evalue threshold for replicon blastn',
                        default=0.00001)
    parser.add_argument('--min_mob_evalue', type=float, required=False,
                        help='Minimum evalue threshold for relaxase tblastn',
                        default=0.00001)
    parser.add_argument('--min_con_evalue', type=float, required=False, help='Minimum evalue threshold for contig blastn',
                        default=0.00001)
    parser.add_argument('--min_rpp_evalue', type=float, required=False,
                        help='Minimum evalue threshold for repetitve elements blastn',
                        default=0.00001)

    parser.add_argument('--min_length', type=int, required=False, help='Minimum length of contigs to classify',
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
                        default=60)

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

    parser.add_argument('--debug', required=False, help='Show debug information', action='store_true')

    parser.add_argument('--plasmid_db', type=str, required=False, help='Reference Database of complete plasmids',
                        default=os.path.join(default_database_dir,
                                             'ncbi_plasmid_full_seqs.fas'))
    parser.add_argument('--plasmid_mash_db', type=str, required=False,
                        help='Companion Mash database of reference database',
                        default=os.path.join(default_database_dir,
                                             'ncbi_plasmid_full_seqs.fas.msh'))
    parser.add_argument('-m', '--plasmid_meta', type=str, required=False,
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
                             'downloaded, they will be downloaded automatically. Defaults to {}'.format(
                            default_database_dir))
    parser.add_argument('--primary_cluster_dist', type=float, required=False,
                        help='Mash distance for assigning primary cluster id 0 - 1', default=0.06)
    parser.add_argument('--secondary_cluster_dist', type=float, required=False,
                        help='Mash distance for assigning primary cluster id 0 - 1',
                        default=0.025)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()


def validate_args(args, logger):
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


def circularize(input_fasta, outdir, logging):
    c = detectCircularity()
    return c.run(input_fasta, outdir, logging)


def filter_blastdf_by_seqs(blast_df, include_seqs, column_name):
    blast_df = blast_df[blast_df[column_name].isin(include_seqs)]
    blast_df = blast_df.reset_index(drop=True)
    return blast_df


def filter_sequences(input_fasta, blastdb, min_ident, min_cov, evalue, min_length, out_dir, blast_results_file,
                     seq_filterfile=None, num_threads=1, max_length=400000):
    blastn(input_fasta, blastdb, min_ident, min_cov, evalue, min_length, out_dir, blast_results_file,
           seq_filterfile, num_threads, max_length)

    if os.path.getsize(blast_results_file) == 0:
        os.remove(blast_results_file)
        return pd.DataFrame()

    blast_df = BlastReader(blast_results_file).df

    if seq_filterfile:
        blast_df = filter_blastdf_by_seqs(blast_df, seq_filterfile)
        blast_df = blast_df.reset_index(drop=True)

    return blast_df


def find_mash_genomes(reference_mash_sketch, fasta_query, outfile, cutoff_distance, num_threads=1):
    cutoff_distance = float(cutoff_distance)

    if (not os.path.isfile(fasta_query)):
        sys.exit(-1)

    if (not os.path.isfile(reference_mash_sketch)):
        sys.exit(-1)
    m = mash()
    distances = parseMash(m.run_mash(reference_mash_sketch, fasta_query, num_threads=num_threads))

    genomes = []

    for query in distances:
        for ref in distances[query]:
            score = distances[query][ref]
            if score < cutoff_distance:
                genomes.append(ref)

    return genomes


def calc_hit_coverage(blast_df, overlap_threshold, reference_sequence_meta):
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
            hit_scores[pID] = {'score': 0, 'length': total_len, 'covered_bases': 0, 'clust_id': clust_id, 'contigs': []}

        hit_scores[pID]['covered_bases'] += aln_length
        hit_scores[pID]['score'] += score
        hit_scores[pID]['contigs'].append(query)
        hit_scores[pID]['contigs'] = list(set(hit_scores[pID]['contigs']))

    return hit_scores


def calc_contig_reference_cov(blast_df, overlap_threshold, reference_sequence_meta):
    blast_df = blast_df.sort_values(['qseqid', 'sseqid', 'qstart', 'qend', 'bitscore'],
                                    ascending=[True, True, True, True, False])
    contig_scores = {}
    size = str(len(blast_df))
    prev_size = 0

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
        cluster_scores[clust_id] += score

    return OrderedDict(sorted(iter(list(cluster_scores.items())), key=lambda x: x[1], reverse=True))


def get_contig_link_counts(reference_hit_coverage, cluster_contig_links):
    cluster_scores = calc_cluster_scores(reference_hit_coverage)

    contig_link_counts = {}
    contig_clust_assoc = {}
    for clust_id in cluster_contig_links:
        contigs = cluster_contig_links[clust_id]
        for contig_id in contigs:
            if not contig_id in contig_link_counts:
                contig_link_counts[contig_id] = 0
                contig_clust_assoc[contig_id] = {}
            contig_link_counts[contig_id] += 1
            contig_clust_assoc[contig_id][clust_id] = cluster_scores[clust_id]

    for contig_id in contig_clust_assoc:
        contig_clust_assoc[contig_id] = OrderedDict(
            sorted(iter(contig_clust_assoc[contig_id].items()), key=lambda x: x[1], reverse=True))

    return OrderedDict(sorted(iter(contig_link_counts.items()), key=lambda x: x[1], reverse=False))


def get_contigs_by_key_value(contig_info, column_key, reasons):
    filtered = []
    for contig_id in contig_info:
        if contig_info[contig_id][column_key] in reasons:
            filtered.append(contig_id)
    return list(set(filtered))


def get_contigs_with_value_set(contig_info, column_key):
    filtered = []
    for contig_id in contig_info:
        if contig_info[contig_id][column_key] is not None and contig_info[contig_id][column_key] != '':
            filtered.append(contig_id)
    return list(set(filtered))


def assign_contigs_to_clusters(contig_blast_df, reference_sequence_meta, contig_info, out_dir, contig_seqs, mash_db,
                               primary_distance, secondary_distance, num_threads=1):
    # Individual reference sequence coverage and overall score along with contig associations
    reference_hit_coverage = calc_hit_coverage(contig_blast_df, 1000, reference_sequence_meta)
    contig_reference_coverage = calc_contig_reference_cov(contig_blast_df, 1000, reference_sequence_meta)
    filtered_contigs = get_contigs_by_key_value(contig_info, 'filtering_reason', ['user', 'chromosome'])
    repetitive_contigs = get_contigs_by_key_value(contig_info, 'filtering_reason', ['repetitve element'])
    circular_contigs = get_contigs_by_key_value(contig_info, 'circularity_status', ['circular'])
    replicon_contigs = get_contigs_with_value_set(contig_info, 'rep_type(s)')
    relaxase_contigs = get_contigs_with_value_set(contig_info, 'relaxase_type(s)')

    contig_list = list(contig_reference_coverage.keys())
    for contig_id in filtered_contigs:
        if contig_id in contig_reference_coverage:
            del(contig_reference_coverage[contig_id])
        if contig_id in contig_reference_coverage:
            del (contig_reference_coverage[contig_id])


    cluster_contig_links = get_seq_links(contig_reference_coverage, reference_sequence_meta)
    contig_link_counts = get_contig_link_counts(reference_hit_coverage, cluster_contig_links)
    cluster_scores = calc_cluster_scores(reference_hit_coverage)

    contig_cluster_scores = {}

    for clust_id in cluster_contig_links:
        for contig_id in cluster_contig_links[clust_id]:
            if contig_id in filtered_contigs:
                continue
            score = cluster_scores[clust_id]
            if not contig_id in contig_cluster_scores:
                contig_cluster_scores[contig_id] = {}
            contig_cluster_scores[contig_id][clust_id] = float(score)

    for contig_id in contig_cluster_scores:
        contig_cluster_scores[contig_id] = OrderedDict(
            sorted(iter(contig_cluster_scores[contig_id].items()), key=lambda x: x[1], reverse=True))

    black_list_clusters = {}
    group_membership = {}
    # assign circular contigs with replicon or relaxase first
    for contig_id in circular_contigs:
        if contig_id in repetitive_contigs or contig_id in filtered_contigs:
            continue
        if contig_id in replicon_contigs or contig_id in relaxase_contigs:
            if contig_id not in contig_cluster_scores:
                continue
            for clust_id in contig_cluster_scores[contig_id]:

                if contig_id not in group_membership:
                    black_list_clusters[clust_id] = contig_id
                    group_membership[contig_id] = clust_id
                    # update cluster scores to remove contigs already assigned
                    if contig_id in contig_reference_coverage:
                        for ref_hit_id in contig_reference_coverage[contig_id]:
                            if ref_hit_id in reference_hit_coverage:
                                reference_hit_coverage[ref_hit_id]['score'] -= contig_reference_coverage[contig_id][
                                    ref_hit_id]
                        del (contig_reference_coverage[contig_id])
                        break

    cluster_scores = calc_cluster_scores(reference_hit_coverage)

    # find plasmids well covered by contigs
    high_confidence_references = {}
    for ref_id in reference_hit_coverage:
        data = reference_hit_coverage[ref_id]
        score = data['score']
        coverage = data['covered_bases'] / data['length']

        if coverage > 0.8:
            high_confidence_references[ref_id] = score

    # Assign contigs according to highly coverged plasmids
    high_confidence_references = OrderedDict(
        sorted(iter(high_confidence_references.items()), key=lambda x: x[1], reverse=True))

    for ref_id in high_confidence_references:
        if not ref_id in reference_sequence_meta:
            continue
        clust_id = reference_sequence_meta[ref_id]['primary_cluster_id']
        if clust_id in black_list_clusters:
            continue
        data = reference_hit_coverage[ref_id]
        contigs = data['contigs']
        for contig_id in contigs:
            if contig_id not in group_membership:
                group_membership[contig_id] = clust_id
                # update cluster scores to remove contigs already assigned
                if contig_id in contig_reference_coverage:
                    for ref_hit_id in contig_reference_coverage[contig_id]:
                        if ref_hit_id in reference_hit_coverage:
                            reference_hit_coverage[ref_hit_id]['score'] -= contig_reference_coverage[contig_id][
                                ref_hit_id]
                    del (contig_reference_coverage[contig_id])
    cluster_scores = calc_cluster_scores(reference_hit_coverage)

    # Assign low linkage contigs first
    for c_id in contig_link_counts:
        count = contig_link_counts[c_id]
        if count > 5:
            break
        scores = contig_cluster_scores[c_id]
        for clust_id in scores:
            score = scores[clust_id]
            if clust_id in black_list_clusters:
                continue

            for contig_id in cluster_contig_links[clust_id]:
                if contig_id not in group_membership:
                    group_membership[contig_id] = clust_id
                    # update cluster scores to remove contigs already assigned
                    if contig_id in contig_reference_coverage:
                        for ref_hit_id in contig_reference_coverage[contig_id]:
                            if ref_hit_id in reference_hit_coverage:
                                reference_hit_coverage[ref_hit_id]['score'] -= contig_reference_coverage[contig_id][
                                    ref_hit_id]
                        del (contig_reference_coverage[contig_id])
            break

    clusters_with_biomarkers = {}
    for clust_id in cluster_contig_links:
        contigs = cluster_contig_links[clust_id]
        for contig_id in contigs:
            if contig_id in relaxase_contigs or contig_id in replicon_contigs:
                if clust_id not in clusters_with_biomarkers:
                    clusters_with_biomarkers[clust_id] = []
                clusters_with_biomarkers[clust_id].append(contig_id)

    max = len(cluster_scores)
    iteration = 0
    while iteration < max:
        iteration += 1

        clust_id = next(iter(cluster_scores))
        if clust_id not in cluster_contig_links:
            del (cluster_scores[clust_id])
            continue

        contigs = cluster_contig_links[clust_id]

        # assign contigs to clusters
        for contig_id in contigs:
            if contig_id not in group_membership:
                group_membership[contig_id] = clust_id
                # update cluster scores to remove contigs already assigned
                if contig_id in contig_reference_coverage:
                    for ref_hit_id in contig_reference_coverage[contig_id]:
                        if ref_hit_id in reference_hit_coverage:
                            reference_hit_coverage[ref_hit_id]['score'] -= contig_reference_coverage[contig_id][
                                ref_hit_id]
                    del (contig_reference_coverage[contig_id])
        cluster_scores = calc_cluster_scores(reference_hit_coverage)
    cluster_links = {}
    for contig_id in group_membership:
        clust_id = group_membership[contig_id]
        if clust_id not in cluster_links:
            cluster_links[clust_id] = []
        cluster_links[clust_id].append(contig_id)
    recon_cluster_dists = get_reconstructed_cluster_dists(mash_db, 0.1, cluster_links, out_dir, contig_seqs,
                                                          num_threads)
    cluster_md5 = {}
    for clust_id in cluster_contig_links:
        contigs = cluster_contig_links[clust_id]

        seq = []
        for contig_id in contigs:
            if contig_id in contig_seqs:
                seq.append(contig_seqs[contig_id])

        seq.sort(key=len)
        cluster_md5[clust_id] = calc_md5(''.join(seq))



    # get lowest distance cluster
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

        if increment:
            counter += 1
            increment = False

        contained_repettive = list(set(repetitive_contigs) & set(cluster_links[clust_id]))

        for contig_id in cluster_links[clust_id]:

            # skip clusters which are just repetitive elemenets
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
                    contig_info[contig_id]['secondary_cluster_id'] = reference_sequence_meta[top_ref_id][
                        'secondary_cluster_id']
            else:
                if len(contained_repettive) != len(cluster_links[clust_id]) and clust_id in clusters_with_biomarkers:
                    contig_info[contig_id]['primary_cluster_id'] = "novel_{}".format(cluster_md5[clust_id])
                    increment = True
                    contig_info[contig_id]['molecule_type'] = 'plasmid'
                    contig_info[contig_id]['mash_nearest_neighbor'] = top_ref_id
                    contig_info[contig_id]['mash_neighbor_distance'] = lowest_dist
                    contig_info[contig_id]['mash_neighbor_identification'] = reference_sequence_meta[top_ref_id][
                        'organism']
                else:
                    contig_info[contig_id]['primary_cluster_id'] = ''
                    contig_info[contig_id]['molecule_type'] = 'chromosome'

    ########working##########
    cluster_links = {}
    for contig_id in contig_info:
        data = contig_info[contig_id]
        if contig_id in filtered_contigs:
            data['molecule_type'] = 'chromosome'
            data['primary_cluster_id'] = ''
            data['secondary_cluster_id'] = ''
            data['mash_nearest_neighbor	'] = ''
            data['mash_neighbor_distance'] = ''
            data['mash_neighbor_identification'] = ''

        if data['molecule_type'] == 'chromosome' or data['primary_cluster_id'] == '':
            continue

        clust_id = data['primary_cluster_id']
        if clust_id not in cluster_links:
            cluster_links[clust_id] = []
        cluster_links[clust_id].append(contig_id)

    recon_cluster_dists = get_reconstructed_cluster_dists(mash_db, 0.1, cluster_links, out_dir, contig_seqs,
                                                          num_threads)
    clusters_with_biomarkers = {}
    for clust_id in cluster_links:
        contigs = cluster_links[clust_id]
        for contig_id in contigs:
            if contig_id in relaxase_contigs or contig_id in replicon_contigs:
                if clust_id not in clusters_with_biomarkers:
                    clusters_with_biomarkers[clust_id] = []
                clusters_with_biomarkers[clust_id].append(contig_id)

    cluster_md5 = {}
    for clust_id in cluster_links:
        contigs = cluster_links[clust_id]
        seq = []
        for contig_id in contigs:
            if contig_id in contig_seqs:
                seq.append(contig_seqs[contig_id])

        seq.sort(key=len)
        cluster_md5[clust_id] = calc_md5(''.join(seq))

    for clust_id in recon_cluster_dists:
        fail = False
        for top_ref_id in recon_cluster_dists[clust_id]:
            lowest_dist = float(recon_cluster_dists[clust_id][top_ref_id])
            if top_ref_id not in reference_sequence_meta:
                fail = True
                continue
            else:
                fail = False
                break

        if fail:
            continue

        contained_repettive = list(set(repetitive_contigs) & set(cluster_links[clust_id]))
        top_ref_cluster_id = reference_sequence_meta[top_ref_id]['primary_cluster_id']

        for contig_id in cluster_links[clust_id]:

            # skip clusters which are just repetitive elemenets
            if len(contained_repettive) == len(cluster_links[clust_id]):

                contig_info[contig_id]['primary_cluster_id'] = ''
                contig_info[contig_id]['molecule_type'] = 'chromosome'
                continue

            if lowest_dist <= primary_distance:
                contig_info[contig_id]['primary_cluster_id'] = top_ref_cluster_id
                contig_info[contig_id]['molecule_type'] = 'plasmid'
                contig_info[contig_id]['mash_nearest_neighbor'] = top_ref_id
                contig_info[contig_id]['mash_neighbor_distance'] = lowest_dist
                contig_info[contig_id]['mash_neighbor_identification'] = reference_sequence_meta[top_ref_id]['organism']

                if lowest_dist <= secondary_distance:
                    contig_info[contig_id]['secondary_cluster_id'] = reference_sequence_meta[top_ref_id][
                        'secondary_cluster_id']
            else:
                if len(contained_repettive) != len(cluster_links[clust_id]) and clust_id in clusters_with_biomarkers:
                    contig_info[contig_id]['primary_cluster_id'] = "novel_{}".format(cluster_md5[clust_id])
                    contig_info[contig_id]['molecule_type'] = 'plasmid'
                    contig_info[contig_id]['mash_nearest_neighbor'] = top_ref_id
                    contig_info[contig_id]['mash_neighbor_distance'] = lowest_dist
                    contig_info[contig_id]['mash_neighbor_identification'] = reference_sequence_meta[top_ref_id][
                        'organism']
                else:
                    contig_info[contig_id]['primary_cluster_id'] = ''
                    contig_info[contig_id]['molecule_type'] = 'chromosome'


    #Fix MD5 assignment for plasmids which had changes
    cluster_membership = {}

    for contig_id in contig_info:
        data = contig_info[contig_id]
        if data['molecule_type'] == 'chromosome' or data['primary_cluster_id'] == '':
            continue
        cluster_id = data['primary_cluster_id']
        if not cluster_id in cluster_membership:
            cluster_membership[cluster_id] = []
        cluster_membership[cluster_id].append(contig_id)


    for clust_id in cluster_membership:
        contigs = cluster_membership[clust_id]
        seq = []

        for contig_id in contigs:
            if contig_id in contig_seqs:
                seq.append(contig_seqs[contig_id])
        if 'novel' in clust_id:
            clust_id = "novel_{}".format(calc_md5(''.join(sorted(seq,key=len))))

        for contig_id in contigs:
            contig_info[contig_id]['primary_cluster_id'] = clust_id


    return evaluate_contig_assignments(contig_info, primary_distance, secondary_distance)


def evaluate_contig_assignments(contig_info, primary_distance, secondary_distance):
    cluster_membership = {}
    biomarker_clusters = {}
    circular_contigs = []
    for contig_id in contig_info:
        data = contig_info[contig_id]
        cluster_id = data['primary_cluster_id']

        if cluster_id == '':
            continue

        if data['circularity_status'] == 'circular':
            circular_contigs.append(contig_id)

        if data['rep_type(s)'] != '' or data['relaxase_type(s)']:
            biomarker_clusters[cluster_id] = ''

        if not cluster_id in cluster_membership:
            cluster_membership[cluster_id] = []

        cluster_membership[cluster_id].append(contig_id)

    for contig_id in contig_info:
        data = contig_info[contig_id]
        cluster_id = data['primary_cluster_id']

        if cluster_id == '':
            continue

        if cluster_id not in biomarker_clusters:
            if contig_id in circular_contigs:
                continue
            if data['mash_neighbor_distance'] > primary_distance:
                contig_info[contig_id]['primary_cluster_id'] = ''
                contig_info[contig_id]['secondary_cluster_id'] = ''
                contig_info[contig_id]['molecule_type'] = 'chromosome'
                contig_info[contig_id]['mash_nearest_neighbor'] = ''
                contig_info[contig_id]['mash_neighbor_distance'] = ''
                contig_info[contig_id]['mash_neighbor_identification'] = ''

    return contig_info


def get_reconstructed_cluster_dists(mash_db, mash_distance, cluster_contig_links, out_dir, contig_seqs, num_threads=1):
    m = mash()
    cluster_dists = {}
    for clust_id in cluster_contig_links:
        contigs = cluster_contig_links[clust_id]
        seq_dict = {}
        tmp_fasta = os.path.join(out_dir, "clust_{}.fasta".format(clust_id))

        for contig_id in contigs:
            if contig_id in contig_seqs:
                seq_dict[contig_id] = contig_seqs[contig_id]
        write_fasta_dict(seq_dict, tmp_fasta)

        distances = parseMash(
            m.run_mash(reference_db=mash_db, input_fasta=tmp_fasta, table=False, num_threads=num_threads))
        os.remove(tmp_fasta)

        for query in distances:
            for ref in distances[query]:
                score = float(distances[query][ref])
                if score <= mash_distance:
                    if clust_id not in cluster_dists:
                        cluster_dists[clust_id] = {}
                    cluster_dists[clust_id][ref] = score

        for clust_id in cluster_dists:
            cluster_dists[clust_id] = OrderedDict(
                sorted(iter(list(cluster_dists[clust_id].items())), key=lambda x: x[1], reverse=False))

    return cluster_dists


def get_seq_links(contig_reference_coverage, reference_sequence_meta):
    reference_clust_members = {}
    for contig_id in contig_reference_coverage:

        for ref_id in contig_reference_coverage[contig_id]:
            if ref_id in reference_sequence_meta:
                clust_id = reference_sequence_meta[ref_id]['primary_cluster_id']
                if not clust_id in reference_clust_members:
                    reference_clust_members[clust_id] = {}
                if not contig_id in reference_clust_members[clust_id]:
                    reference_clust_members[clust_id][contig_id] = 0
                reference_clust_members[clust_id][contig_id] += 1

    return reference_clust_members


def update_group_members(target_contigs, group_membership, contig_reference_coverage, reference_sequence_meta,
                         group_membership_key, reference_seq_key):
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


def filter_contig_df_by_index(indicies, contig_blast_df, reference_hit_coverage):
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
            logging.warning("{} not found".format(pID))

    return reference_hit_coverage


def main():
    args = parse_args()

    if args.debug:
        logger = init_console_logger(3)
    else:
        logger = init_console_logger(2)

    logger.info("MOB-recon version {} ".format(__version__))
    logger.debug("Debug log reporting set on successfully")

    check_dependencies(logger)
    validate_args(args, logger)
    max_contig_size = args.max_contig_size
    max_plasmid_size = args.max_plasmid_size
    keep_tmp = args.keep_tmp
    input_fasta = args.infile
    out_dir = args.outdir
    num_threads = args.num_threads
    tmp_dir = os.path.join(out_dir, '__tmp')
    fixed_fasta = os.path.join(tmp_dir, 'fixed.input.fasta')
    chromosome_file = os.path.join(out_dir, 'chromosome.fasta')
    replicon_blast_results = os.path.join(tmp_dir, 'replicon_blast_results.txt')
    mob_blast_results = os.path.join(tmp_dir, 'mob_blast_results.txt')
    mpf_blast_results = os.path.join(tmp_dir, 'mpf_blast_results.txt')
    orit_blast_results = os.path.join(tmp_dir, 'orit_blast_results.txt')
    repetitive_blast_results = os.path.join(tmp_dir, 'repetitive_blast_results.txt')
    contig_blast_results = os.path.join(tmp_dir, 'contig_blast_results.txt')
    prefix = None
    if args.prefix is not None:
        prefix = args.prefix
    contig_report = os.path.join(out_dir, 'contig_report.txt')
    if prefix is not None:
        contig_report = os.path.join(out_dir, "{}.contig_report.txt".format(prefix))
        chromosome_file = os.path.join(out_dir, "{}.chromosome.fasta".format(prefix))
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
        verify_init(logger, database_dir)
    else:
        plasmid_ref_db = os.path.join(database_dir, 'ncbi_plasmid_full_seqs.fas')
        mob_ref = os.path.join(database_dir, 'mob.proteins.faa')
        mash_db = os.path.join(database_dir, 'ncbi_plasmid_full_seqs.fas.msh')
        replicon_ref = os.path.join(database_dir, 'rep.dna.fas')
        plasmid_meta = os.path.join(database_dir, 'clusters.txt')
        repetitive_mask_file = os.path.join(database_dir, 'repetitive.dna.fas')
        mpf_ref = os.path.join(database_dir, 'mpf.proteins.faa')
        plasmid_orit = os.path.join(database_dir, 'orit.fas')
    ETE3DBTAXAFILE = os.path.abspath(database_dir + "/taxa.sqlite")

    LIT_PLASMID_TAXONOMY_FILE = os.path.join(database_dir, "host_range_literature_plasmidDB.txt")
    NCBI_PLASMID_TAXONOMY_FILE = plasmid_meta

    if args.sample_id is None:
        sample_id = re.sub(r"\.(fasta|fas|fa){1,1}", "", os.path.basename(args.infile))
    else:
        sample_id = args.sample_id


    run_overhang = args.run_overhang
    unicycler_contigs = args.unicycler_contigs

    # initialize analysis directory
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir, 0o755)

    elif not args.force:
        logger.error("Error output directory exists, please specify a new directory or use --force to overwrite")
        sys.exit(-1)
    else:
        shutil.rmtree(args.outdir)
        os.makedirs(args.outdir, 0o755)

    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir, 0o755)
    else:
        shutil.rmtree(tmp_dir)
        os.makedirs(tmp_dir, 0o755)

    # Initialize clustering distance thresholds
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
    min_mpf_ident = float(args.min_mob_ident)  # Left in for future if we decide to allow modification
    min_con_ident = float(args.min_con_ident)
    min_rpp_ident = float(args.min_rpp_ident)

    min_rep_cov = float(args.min_rep_cov)
    min_mob_cov = float(args.min_mob_cov)
    min_mpf_cov = float(args.min_mob_cov)  # Left in for future if we decide to allow modification
    min_con_cov = float(args.min_con_cov)
    min_rpp_cov = float(args.min_rpp_cov)

    min_rep_evalue = float(args.min_rep_evalue)
    min_mob_evalue = float(args.min_mob_evalue)
    min_mpf_evalue = float(args.min_mob_evalue)  # Left in for future if we decide to allow modification
    min_con_evalue = float(args.min_con_evalue)
    min_rpp_evalue = float(args.min_rpp_evalue)


    # Parse reference cluster information
    reference_sequence_meta = read_sequence_info(plasmid_meta, MOB_CLUSTER_INFO_HEADER)

    # process input fasta
    logger.info('Writing cleaned header input fasta file from {} to {}'.format(input_fasta, fixed_fasta))
    id_mapping = fix_fasta_header(input_fasta, fixed_fasta)
    contig_seqs = read_fasta_dict(fixed_fasta)
    contig_info = {}
    br = BlastRunner(fixed_fasta, tmp_dir)
    br.makeblastdb(fixed_fasta, dbtype='nucl', logging=logging)
    del (br)

    # Detect circular sequences
    circular_contigs = {}
    if run_overhang:
        logger.info('Running internal circular contig detection on {}'.format(fixed_fasta))
        circular_contigs = circularize(fixed_fasta, tmp_dir, logging)

    for id in contig_seqs:
        seq = contig_seqs[id]
        contig_info[id] = {}
        for feature in MOB_RECON_INFO_HEADER:
            contig_info[id][feature] = ''
        contig_info[id]['md5'] = calc_md5(seq)
        contig_info[id]['gc'] = gc_fraction(seq)
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
            if 'circular=true' in id or id in circular_contigs:
                contig_info[id]['circularity_status'] = 'circular'
            elif id not in circular_contigs:
                contig_info[id]['circularity_status'] = 'incomplete'

        if contig_info[id]['circularity_status'] == '':
            contig_info[id]['circularity_status'] = 'not tested'

    # Blast reference databases

    identify_biomarkers(contig_info, fixed_fasta, tmp_dir, min_length, logging, \
                        replicon_ref, min_rep_ident, min_rep_cov, min_rep_evalue, replicon_blast_results, \
                        mob_ref, min_mob_ident, min_mob_cov, min_mob_evalue, mob_blast_results, \
                        mpf_ref, min_mpf_ident, min_mpf_cov, min_mpf_evalue, mpf_blast_results, \
                        repetitive_mask_file, min_rpp_ident, min_rpp_cov, min_rpp_evalue, \
                        plasmid_orit, orit_blast_results, repetitive_blast_results, \
                        num_threads=num_threads)

    # Filtering contigs against chromosome database

    chrom_filter = False

    if args.genome_filter_db_prefix:
        chrom_filter = True
        genome_filter_db_prefix = args.genome_filter_db_prefix
        logger.info('Genome filter sequences provided: {}'.format(genome_filter_db_prefix))
        matched = (glob.glob(genome_filter_db_prefix + "*"))
        extensions = ['nsq', 'nin', 'nhr']
        found = [0, 0, 0]
        for f in matched:
            for i in range(0, len(extensions)):
                e = extensions[i]
                if e in f:
                    found[i] += 1

        for i in found:
            if i == 0:
                logger.error('Error blast database not found with prefix: {}'.format(genome_filter_db_prefix))
                sys.exit()
        if not os.path.isfile(genome_filter_db_prefix + '.msh'):
            logger.error('Error mash sketch not found with prefix: {}'.format(genome_filter_db_prefix))
            sys.exit()

    if chrom_filter:
        cutoff_distance = float(args.mash_genome_neighbor_threshold)
        chr_mash_dists = os.path.join(tmp_dir, 'mash_chr_dists.txt')
        chr_blast_filter = os.path.join(tmp_dir, 'chr_contig_filter_report.txt')
        chr_mash_sketch = genome_filter_db_prefix + ".msh"
        close_genome_reps = find_mash_genomes(chr_mash_sketch, fixed_fasta, chr_mash_dists, cutoff_distance,
                                              num_threads)

        if len(close_genome_reps) > 0:
            logger.info('Found close genome matches: {}'.format(",".join(close_genome_reps)))
            seq_id_file = os.path.join(tmp_dir, "seqids.txt")
            sf = open(seq_id_file, 'w')
            for s in close_genome_reps:
                sf.write("{}\n".format(s))
            sf.close()

            # fix labels to match the seq id format parsed by makeblastdb
            for i in range(0, len(close_genome_reps)):
                close_genome_reps[i] = "ref|{}|".format(close_genome_reps[i])

            blastn(input_fasta=fixed_fasta, blastdb=genome_filter_db_prefix, min_ident=min_con_ident,
                   min_cov=min_con_cov,
                   evalue=min_con_evalue, min_length=min_length, out_dir=tmp_dir,
                   blast_results_file=chr_blast_filter, num_threads=num_threads, logging=logging,
                   seq_filterfile=seq_id_file)

            chromosome_filter_seqs = BlastReader(chr_blast_filter, logging).df.drop(0)['qseqid'].tolist()

            for contig_id in chromosome_filter_seqs:
                if contig_id in contig_info:
                        contig_info[contig_id]['filtering_reason'] = 'chromosome'
                        logger.info('Filtering contig: {} due to inclusion in genome filter {}'.format(contig_id,
                                                                                                       genome_filter_db_prefix))
                else:
                    logger.error('Contig: {} not found in contig_df this is likely an error'.format(contig_id))

            del (chromosome_filter_seqs)

        else:
            logger.info('No close genome matches found')

    # Filter out sequences based on user filter
    if args.filter_db:
        filter_db = args.filter_db
        logger.info('Filter sequences provided: {}'.format(filter_db))
        if not os.path.isfile(filter_db + '.nsq') or \
                os.path.isfile(filter_db + '.nin') or \
                os.path.isfile(filter_db + '.nhr'):
            br = BlastRunner(filter_db, os.path.dirname(filter_db))
            br.makeblastdb(filter_db, dbtype='nucl', logging=logging)

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
                contig_info[contig_id]['filtering_reason'] = 'user'
                logger.info(
                    'Filtering contig: {} due to inclusion in genome filter {}'.format(contig_id, filter_db))
            else:
                logger.error('Contig: {} not found in contig_df this is likely an error'.format(contig_id))

        del (user_filter_seqs)


    # blast plasmid database
    logging.info("Blasting contigs against reference sequence db: {}".format(plasmid_ref_db))
    blastn(input_fasta=fixed_fasta, blastdb=plasmid_ref_db, min_ident=min_con_ident, min_cov=min_con_cov,
           evalue=min_con_evalue, min_length=min_length, out_dir=tmp_dir,
           blast_results_file=contig_blast_results, num_threads=num_threads, logging=logging, seq_filterfile=None)

    logging.info("Filtering contig blast results: {}".format(contig_blast_results))
    contig_blast_df = BlastReader(contig_blast_results, logging).df

    if len(contig_blast_df) > 0:
        contig_blast_df = fixStart(contig_blast_df.drop(0)).sort_values(
            ['sseqid', 'qseqid', 'sstart', 'send', 'bitscore'])
        contig_blast_df = filter_overlaping_records(contig_blast_df, 500, 'qseqid', 'qstart', 'qend', 'bitscore')
        contig_blast_df.reset_index(drop=True)
        # remove blast formatting of seq id
        for index, row in contig_blast_df.iterrows():
            line = row['sseqid'].split('|')
            if len(line) >= 2:
                contig_blast_df.at[index, 'sseqid'] = line[1]
        contig_blast_df = contig_blast_df[contig_blast_df.sseqid.isin(list(reference_sequence_meta.keys()))]
        contig_blast_df.reset_index(drop=True)
        logging.info("Assigning contigs to plasmid groups")
        contig_info = assign_contigs_to_clusters(contig_blast_df, reference_sequence_meta, contig_info, tmp_dir,
                                                 contig_seqs, mash_db, primary_distance, secondary_distance,
                                                 num_threads)

    # Triage Novel plasmids with biomarkers but not enough similaririty for valid contig hits
    for contig_id in contig_info:
        data = contig_info[contig_id]

        if data['circularity_status'] == 'circular' and data['repetitive_dna_id'] == '' and data['primary_cluster_id'] == '':
            rep_types = data['rep_type(s)']
            mob_types = data['relaxase_type(s)']

            if rep_types == '' and mob_types == '':
                continue
            md5 = data['md5']
            primary_clust_id = "{}".format(md5)
            data['primary_cluster_id'] = "novel_{}".format(primary_clust_id)
            data['secondary_cluster_id'] = ""
            data['molecule_type'] = 'plasmid'
            data['mash_nearest_neighbor'] = ''
            data['mash_neighbor_distance'] = ''
            data['mash_neighbor_identification'] = ''

            cluster_dists = get_reconstructed_cluster_dists(mash_db, 0.1, {primary_clust_id: [contig_id]}, out_dir,
                                                            contig_seqs, num_threads)

            lowest_dist = 1
            for clust_id in cluster_dists:
                fail = False

                lowest_dist = 1
                for top_ref_id in cluster_dists[clust_id]:
                    lowest_dist = cluster_dists[clust_id][top_ref_id]
                    if top_ref_id not in reference_sequence_meta:
                        fail = True
                        continue
                    else:
                        fail = False
                        break
                if fail:
                    continue

                if lowest_dist <= 0.1:
                    data['mash_nearest_neighbor'] = top_ref_id
                    data['mash_neighbor_distance'] = lowest_dist
                    data['mash_neighbor_identification'] = reference_sequence_meta[top_ref_id]['organism']

            contig_info[contig_id] = data

    #Assign sequences to chromosome which exceed length limit
    plasmid_sizes = {}
    for contig_id in contig_info:
        length = contig_info[contig_id]['size']
        if contig_info[contig_id]['primary_cluster_id'] != '' and contig_info[contig_id]['mash_nearest_neighbor'] == '':
            contig_info[contig_id]['primary_cluster_id'] = ''
            contig_info[contig_id]['molecule_type'] = 'chromosome'
            data['filtering_reason'] = "Contig length > {}bp".format(length)
        if contig_info[contig_id]['size'] > max_contig_size:
            contig_info[contig_id]['molecule_type'] == 'chromosome'
        if contig_info[contig_id]['molecule_type']  == 'plasmid':
            pID = contig_info[contig_id]['primary_cluster_id']

            if not pID in plasmid_sizes:
                plasmid_sizes[pID] = 0
            plasmid_sizes[pID]+= length

    for contig_id in contig_info:
        if contig_info[contig_id]['molecule_type'] == 'plasmid':
            pID = contig_info[contig_id]['primary_cluster_id']
            length = contig_info[contig_id]['size']
            if plasmid_sizes[pID] > max_plasmid_size:
                contig_info[contig_id]['primary_cluster_id'] = ''
                contig_info[contig_id]['secondary_cluster_id'] = ''
                contig_info[contig_id]['molecule_type'] = 'chromosome'
                data['filtering_reason'] = "Reconstructed plasmid length > {}bp".format(length)

    results = []
    contig_memberships = {'chromosome': {}, 'plasmid': {}}

    for contig_id in contig_info:
        data = contig_info[contig_id]
        if contig_id in id_mapping:
            original_id = id_mapping[contig_id]
        data['contig_id'] = original_id
        if data['primary_cluster_id'] != '' and data['mash_nearest_neighbor'] == '':
            data['primary_cluster_id'] = ''
            data['molecule_type'] = 'chromosome'
        if data['filtering_reason'] == 'none':
            data['filtering_reason'] = ''

        if data['molecule_type'] == 'chromosome':
            contig_memberships['chromosome'][contig_id] = data
        else:
            clust_id = data['primary_cluster_id']
            if not clust_id in contig_memberships['plasmid']:
                contig_memberships['plasmid'][clust_id] = {}
            contig_memberships['plasmid'][clust_id][contig_id] = data

        results.append(data)

    # Write contig report
    logging.info("Writting contig results to {}".format(contig_report))
    if len(results) > 0:

        if len(contig_memberships['plasmid']) > 0:
            ncbi = dict_from_alt_key_list(
                read_file_to_dict(NCBI_PLASMID_TAXONOMY_FILE, MOB_CLUSTER_INFO_HEADER, separater="\t"),
                "sample_id")
            lit = dict_from_alt_key_list(
                read_file_to_dict(LIT_PLASMID_TAXONOMY_FILE, LIT_PLASMID_TAXONOMY_HEADER, separater="\t"), "sample_id")

            mobtyper_report = os.path.join(out_dir, "mobtyper_results.txt")
            if prefix is not None:
                mobtyper_report = os.path.join(out_dir, "{}.mobtyper_results.txt".format(prefix))
            build_mobtyper_report(contig_memberships['plasmid'], out_dir, mobtyper_report,contig_seqs, ncbi, lit, ETE3DBTAXAFILE, database_dir)

        writeReport(results, MOB_RECON_INFO_HEADER, contig_report)

        logging.info("Writting chromosome sequences to {}".format(chromosome_file))
        chr_fh = open(chromosome_file, 'w')
        for contig_id in contig_memberships['chromosome']:
            if contig_id in contig_seqs:
                original_id = id_mapping[contig_id]
                chr_fh.write(">{}\n{}\n".format(original_id, contig_seqs[contig_id]))
        chr_fh.close()

    #fix plasmid fastas
    clusters = list(contig_memberships['plasmid'].keys())
    for cluster in clusters:
        file = os.path.join(out_dir, "plasmid_{}.fasta".format(cluster))
        if prefix is None:
            update = False
        else:
            update = True
        seqs = read_fasta_dict(file)
        fh = open(file,'w')
        for seq_id in seqs:
            seq = seqs[seq_id]
            if seq_id in id_mapping:
                seq_id = id_mapping[seq_id]
            fh.write(">{}\n{}\n".format(seq_id,seq))
        fh.close()
        if update:
            os.rename(file, os.path.join(out_dir, "{}.plasmid_{}.fasta".format(prefix,cluster)))

    #Peform MGE detection
    mge_results = blast_mge(fixed_fasta, repetitive_mask_file, tmp_dir, min_length,
                            logging, min_rpp_ident, min_rpp_cov, min_rpp_evalue,num_threads)
    mge_report_file = os.path.join(out_dir,"mge.report.txt")
    if prefix is not None:
        mge_report_file = os.path.join(out_dir, "{}.mge_report.txt".format(prefix))

    writeMGEresults(contig_memberships, mge_results, mge_report_file)

    biomarker_params = {
        'oriT': {
            'file':orit_blast_results,
            'min_length': 80,
            'min_cov':min_rep_cov,
            'min_hsp_cov': 25,
            'evalue':min_rep_evalue,
            'min_ident':min_rep_ident
        },
        'replicon': {
            'file':replicon_blast_results,
            'min_length': 80,
            'min_cov':min_rep_cov,
            'min_hsp_cov': 25,
            'evalue':min_rep_evalue,
            'min_ident':min_rep_ident
        },
        'relaxase': {
            'file':mob_blast_results,
            'min_length': 40,
            'min_cov':min_mob_cov,
            'min_hsp_cov': 25,
            'evalue':min_mob_evalue,
            'min_ident':min_mob_ident
        },        
        'mate-pair-formation': {
            'file':mpf_blast_results,
            'min_length': 40,
            'min_cov':min_mpf_cov,
            'min_hsp_cov': 25,
            'evalue':min_mpf_evalue,
            'min_ident':min_mpf_ident
        },

    }

    biomarker_df = create_biomarker_dataframe(biomarker_params,id_mapping,logging)
    if biomarker_df.empty == False:
        logging.info(f"Writting plasmid biomarkers report contaning {biomarker_df.shape[0]} biomarkers")
        biomarker_df.to_csv(os.path.join(out_dir,"biomarkers.blast.txt"),header=True,sep="\t",index=False)
    else:
        logging.warning("Plasmid biomarkers report is empty and will not be reported")
    
    if not keep_tmp:
        logging.info("Cleaning up temporary files {}".format(tmp_dir))
        shutil.rmtree(tmp_dir)
    logger.info("MOB-recon completed and results written to {}".format(out_dir))


# call main function
if __name__ == '__main__':
    main()
