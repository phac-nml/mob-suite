#!/usr/bin/env python3

import logging
import os, re, shutil, sys, tempfile
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from mob_suite.version import __version__
import mob_suite.mob_init
from collections import OrderedDict
from operator import itemgetter
from mob_suite.blast import BlastRunner
from mob_suite.wrappers import mash
from Bio import SeqIO

from mob_suite.utils import fix_fasta_header, \
    calcFastaStats, \
    verify_init, \
    check_dependencies, \
    read_sequence_info, \
    writeReport, \
    sort_biomarkers, \
    calc_md5, \
    gc_fraction, \
    read_fasta_dict, \
    identify_biomarkers, \
    parseMash, \
    determine_mpf_type, \
    hostrange, \
    dict_from_alt_key_list, \
    read_file_to_dict, \
    blast_mge, \
    writeMGEresults, \
    create_biomarker_dataframe

from mob_suite.constants import MOB_TYPER_REPORT_HEADER, \
    MOB_CLUSTER_INFO_HEADER, \
    default_database_dir, \
    LIT_PLASMID_TAXONOMY_HEADER


def init_console_logger(lvl=2):
    root = logging.getLogger()

    LOG_FORMAT = '%(asctime)s %(name)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]

    report_lvl = logging_levels[lvl]
    root.setLevel(report_lvl)  # set root logger level

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)
    return logging.getLogger(__name__)


def parse_args():
    "Parse the input arguments, use '-h' for help"

    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="MOB-Typer: Plasmid typing and mobility prediction: {}".format(
            __version__), formatter_class=CustomFormatter)
    parser.add_argument('-i', '--infile', type=str, required=True, help='Input assembly fasta file to process')
    parser.add_argument('-o', '--out_file', type=str, required=True, help='Output file to write results')
    parser.add_argument('--biomarker_report_file', type=str, required=False, help='Output file for biomarker blast results')
    parser.add_argument('-g', '--mge_report_file', type=str, required=False, help='Output file for MGE results')
    parser.add_argument('-a', '--analysis_dir', type=str, required=False,
                        help='Working directory for storing temporary results')
    parser.add_argument('-n', '--num_threads', type=int, required=False, help='Number of threads to be used', default=1)
    parser.add_argument('-s', '--sample_id', type=str, required=False, help='Sample Prefix for reports')
    parser.add_argument('-x', '--multi', required=False, help='Treat each sequence as an independant plasmid',
                        action='store_true')
    parser.add_argument('--min_rep_evalue', type=float, required=False,
                        help='Minimum evalue threshold for replicon blastn',
                        default=0.00001)
    parser.add_argument('--min_mob_evalue', type=float, required=False,
                        help='Minimum evalue threshold for relaxase tblastn',
                        default=0.00001)
    parser.add_argument('--min_con_evalue', type=float, required=False, help='Minimum evalue threshold for contig blastn',
                        default=0.00001)

    parser.add_argument('--min_length', type=str, required=False, help='Minimum length of blast hits',
                        default=500)
    parser.add_argument('--min_rep_ident', type=int, required=False, help='Minimum sequence identity for replicons',
                        default=80)
    parser.add_argument('--min_mob_ident', type=int, required=False, help='Minimum sequence identity for relaxases',
                        default=80)
    parser.add_argument('--min_con_ident', type=int, required=False, help='Minimum sequence identity for contigs',
                        default=80)
    parser.add_argument('--min_rpp_ident', type=int, required=False,
                        help='Minimum sequence identity for MGE', default=80)

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
                        help='Minimum percentage coverage of MGE',
                        default=80)
    parser.add_argument('--min_rpp_evalue', type=float, required=False,
                        help='Minimum evalue threshold for repetitve elements blastn',
                        default=0.00001)

    parser.add_argument('--min_overlap', type=int, required=False,
                        help='Minimum overlap of fragments',
                        default=10)

    parser.add_argument('-k', '--keep_tmp', required=False, help='Do not delete temporary file directory',
                        action='store_true')

    parser.add_argument('--debug', required=False, help='Show debug information', action='store_true')

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
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()


def initMOBTyperReportTemplate(header):
    data = {}
    for i in header:
        data[i] = '-'
    return data


def main():
    args = parse_args()

    if args.debug:
        logger = init_console_logger(3)
    else:
        logger = init_console_logger(2)

    logger.info('Running Mob-typer version {}'.format(__version__))

    logger.info('Processing fasta file {}'.format(args.infile))

    if not os.path.isfile(args.infile):
        logger.info('Error, fasta file does not exist {}'.format(args.infile))
        sys.exit()

    if not args.analysis_dir:
        tmp_dir = tempfile.TemporaryDirectory(dir=tempfile.gettempdir()).name
    else:
        tmp_dir = args.analysis_dir


    if not os.path.isdir(tmp_dir):
        logger.info('Creating analysis directory {}'.format(tmp_dir))
        os.makedirs(tmp_dir, 0o755)


    if not isinstance(args.num_threads, int):
        logger.info('Error number of threads must be an integer, you specified "{}"'.format(args.num_threads))

    database_dir = os.path.abspath(args.database_directory)

    if args.sample_id is None:
        sample_id = re.sub(r"\.(fasta|fa|fas){1,1}", "", os.path.basename(args.infile))
    else:
        sample_id = args.sample_id

    # Script arguments
    input_fasta = args.infile
    report_file = args.out_file
    mge_report_file = args.mge_report_file
    num_threads = int(args.num_threads)
    keep_tmp = args.keep_tmp
    biomarker_report_file = args.biomarker_report_file

    if args.multi:
        multi = True
    else:
        multi = False

    if not (args.primary_cluster_dist >= 0 and args.primary_cluster_dist <= 1):
        logging.error('Error distance thresholds must be between 0 - 1: {}'.format(args.primary_cluster_dist))
        sys.exit()
    else:
        primary_distance = float(args.primary_cluster_dist)

    min_length = int(args.min_length)
    ETE3DBTAXAFILE = os.path.abspath(database_dir + "/taxa.sqlite")

    if database_dir == default_database_dir:
        mob_ref = args.plasmid_mob
        mash_db = args.plasmid_mash_db
        replicon_ref = args.plasmid_replicons
        plasmid_meta = args.plasmid_meta
        mpf_ref = args.plasmid_mpf
        plasmid_orit = args.plasmid_orit
        repetitive_mask_file = args.repetitive_mask
        verify_init(logger, database_dir)
    else:
        mob_ref = os.path.join(database_dir, 'mob.proteins.faa')
        mash_db = os.path.join(database_dir, 'ncbi_plasmid_full_seqs.fas.msh')
        replicon_ref = os.path.join(database_dir, 'rep.dna.fas')
        plasmid_meta = os.path.join(database_dir, 'clusters.txt')
        mpf_ref = os.path.join(database_dir, 'mpf.proteins.faa')
        plasmid_orit = os.path.join(database_dir, 'orit.fas')
        repetitive_mask_file = os.path.join(database_dir, 'repetitive.dna.fas')


    LIT_PLASMID_TAXONOMY_FILE = os.path.join(database_dir, "host_range_literature_plasmidDB.txt")
    NCBI_PLASMID_TAXONOMY_FILE = plasmid_meta

    fixed_fasta = os.path.join(tmp_dir, 'fixed.input.fasta')
    replicon_blast_results = os.path.join(tmp_dir, 'replicon_blast_results.txt')
    mob_blast_results = os.path.join(tmp_dir, 'mobtyper_blast_results.txt')
    mpf_blast_results = os.path.join(tmp_dir, 'mpf_blast_results.txt')
    orit_blast_results = os.path.join(tmp_dir, 'orit_blast_results.txt')
    repetitive_blast_results = os.path.join(tmp_dir, 'repetitive_blast_results.txt')

    if os.path.isfile(mob_blast_results):
        os.remove(mob_blast_results)
    if os.path.isfile(mpf_blast_results):
        os.remove(mpf_blast_results)
    if os.path.isfile(orit_blast_results):
        os.remove(orit_blast_results)
    if os.path.isfile(replicon_blast_results):
        os.remove(replicon_blast_results)

    # Input numeric params

    min_rep_ident = float(args.min_rep_ident)
    min_mob_ident = float(args.min_mob_ident)
    min_ori_ident = float(args.min_rep_ident)
    min_mpf_ident = float(args.min_mob_ident)
    min_rpp_ident = float(args.min_rpp_ident)

    idents = {'min_rep_ident': min_rep_ident, 'min_mob_ident': min_mob_ident, 'min_ori_ident': min_ori_ident}

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
    min_ori_cov = float(args.min_rep_cov)
    min_mpf_cov = float(args.min_mob_cov)
    min_rpp_cov = float(args.min_rpp_cov)

    covs = {'min_rep_cov': min_rep_cov, 'min_mob_cov': min_mob_cov, 'min_con_cov': min_ori_cov,
            'min_rpp_cov': min_ori_cov}

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
    min_ori_evalue = float(args.min_rep_evalue)
    min_mpf_evalue = float(args.min_mob_evalue)
    min_rpp_evalue = float(args.min_rpp_evalue)

    evalues = {'min_rep_evalue': min_rep_evalue, 'min_mob_evalue': min_mob_evalue, 'min_con_evalue': min_ori_evalue}

    for param in evalues:

        value = float(evalues[param])

        if value > 1:
            logger.error("Error: {} is too high, please specify an float evalue between 0 to 1".format(param))
            sys.exit(-1)

    check_dependencies(logger)

    needed_dbs = [replicon_ref, mob_ref, mash_db, mpf_ref]

    for db in needed_dbs:
        if (not os.path.isfile(db)):
            logger.info('Warning! Needed database missing "{}"'.format(db))
            mob_suite.mob_init.main()

    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir, 0o755)

    # Get cluster information
    reference_sequence_meta = read_sequence_info(plasmid_meta, MOB_CLUSTER_INFO_HEADER)

    # initilize master record tracking
    id_mapping = fix_fasta_header(input_fasta, fixed_fasta)
    contig_info = {}
    with open(fixed_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id = str(record.id)
            contig_info[id] = {}
            for feature in MOB_TYPER_REPORT_HEADER:
                contig_info[id][feature] = ''
            seq = str(record.seq)
            contig_info[id]['md5'] = calc_md5(seq)
            contig_info[id]['gc'] = gc_fraction(seq)
            contig_info[id]['size'] = len(seq)
            contig_info[id]['contig_id'] = id
            contig_info[id]['sample_id'] = sample_id
    handle.close()

    # Makeblastdb
    blast_runner = BlastRunner(fixed_fasta, tmp_dir)
    build_success = blast_runner.makeblastdb(fixed_fasta, 'nucl', logging=logging)
    if build_success == False:
        logging.error("Could not build blast database, check error messages..cannot continue")
        sys.exit()

    # run individual marker blasts

    contig_info = identify_biomarkers(contig_info, fixed_fasta, tmp_dir, 25, logging, \
                                      replicon_ref, min_rep_ident, min_rep_cov, min_rep_evalue, replicon_blast_results, \
                                      mob_ref, min_mob_ident, min_mob_cov, min_mob_evalue, mob_blast_results, \
                                      mpf_ref, min_mpf_ident, min_mpf_cov, min_mpf_evalue, mpf_blast_results, \
                                      None, None, None, None, \
                                      plasmid_orit, orit_blast_results, repetitive_blast_results, \
                                      num_threads=num_threads)

    m = mash()
    mobtyper_results = []

    mash_input_fasta = fixed_fasta + '.msh'

    ncbi = dict_from_alt_key_list(
        read_file_to_dict(NCBI_PLASMID_TAXONOMY_FILE, MOB_CLUSTER_INFO_HEADER, separater="\t"),
        "sample_id")
    lit = dict_from_alt_key_list(
        read_file_to_dict(LIT_PLASMID_TAXONOMY_FILE, LIT_PLASMID_TAXONOMY_HEADER, separater="\t"), "sample_id")

    if multi:
        m.mashsketch(input_fasta=fixed_fasta, output_path=mash_input_fasta, sketch_ind=True, num_threads=num_threads)
        mash_results = parseMash(
            m.run_mash(reference_db=mash_db, input_fasta=mash_input_fasta, table=False, num_threads=num_threads))

        for seq_id in mash_results:
            record = {}
            for field in MOB_TYPER_REPORT_HEADER:
                if field in contig_info[seq_id]:
                    record[field] = contig_info[seq_id][field]
                else:
                    record[field] = ''
            record['sample_id'] = seq_id
            record['num_contigs'] = 1
            distances = OrderedDict(sorted(mash_results[seq_id].items(), key=itemgetter(1), reverse=False))

            for mash_neighbor_id in distances:
                dist = distances[mash_neighbor_id]
                if mash_neighbor_id not in reference_sequence_meta:
                    continue
                else:
                    record['mash_nearest_neighbor'] = mash_neighbor_id
                    record['mash_neighbor_distance'] = dist
                    record['primary_cluster_id'] = reference_sequence_meta[mash_neighbor_id]['primary_cluster_id']
                    record['secondary_cluster_id'] = reference_sequence_meta[mash_neighbor_id]['secondary_cluster_id']
                    record['mash_neighbor_identification'] = reference_sequence_meta[mash_neighbor_id]['organism']
                    break
            mobtyper_results.append(record)

    else:
        m.mashsketch(input_fasta=fixed_fasta, output_path=mash_input_fasta, sketch_ind=False, num_threads=num_threads)
        mash_results = parseMash(
            m.run_mash(reference_db=mash_db, input_fasta=mash_input_fasta, table=False, num_threads=num_threads))
        record = {}

        for field in MOB_TYPER_REPORT_HEADER:
            record[field] = ''

        record['sample_id'] = sample_id
        fastaSeqStats = calcFastaStats(fixed_fasta)
        record['md5'] = fastaSeqStats['md5']
        record['size'] = fastaSeqStats['size']
        record['num_contigs'] = fastaSeqStats['num_seq']
        record['gc'] = fastaSeqStats['gc_content']
        record['mash_nearest_neighbor'] = '-'
        record['mash_neighbor_distance'] = 1
        record['primary_cluster_id'] = '-'
        record['secondary_cluster_id'] = '-'
        record['mash_neighbor_identification'] = '-'

        for seq_id in mash_results:
            distances = OrderedDict(sorted(mash_results[seq_id].items(), key=itemgetter(1), reverse=False))
            mash_neighbor_id = next(iter(distances))
            dist = distances[mash_neighbor_id]
            if mash_neighbor_id not in reference_sequence_meta:
                continue
            record['mash_nearest_neighbor'] = mash_neighbor_id
            record['mash_neighbor_distance'] = dist
            record['primary_cluster_id'] = reference_sequence_meta[mash_neighbor_id]['primary_cluster_id']
            record['secondary_cluster_id'] = reference_sequence_meta[mash_neighbor_id]['secondary_cluster_id']
            record['mash_neighbor_identification'] = reference_sequence_meta[mash_neighbor_id]['organism']

        record['rep_type(s)'] = []
        record['rep_type_accession(s)'] = []
        record['relaxase_type(s)'] = []
        record['relaxase_type_accession(s)'] = []
        record['mpf_type'] = []
        record['mpf_type_accession(s)'] = []
        record['orit_type(s)'] = []
        record['orit_accession(s)'] = []

        for seq_id in contig_info:
            record['rep_type(s)'].append(contig_info[seq_id]['rep_type(s)'])
            record['rep_type_accession(s)'].append(contig_info[seq_id]['rep_type_accession(s)'])
            record['relaxase_type(s)'].append(contig_info[seq_id]['relaxase_type(s)'])
            record['relaxase_type_accession(s)'].append(contig_info[seq_id]['relaxase_type_accession(s)'])
            record['mpf_type'].append(contig_info[seq_id]['mpf_type'])
            record['mpf_type_accession(s)'].append(contig_info[seq_id]['mpf_type_accession(s)'])
            record['orit_type(s)'].append(contig_info[seq_id]['orit_type(s)'])
            record['orit_accession(s)'].append(contig_info[seq_id]['orit_accession(s)'])

        for field in record:
            tmp = []
            if record[field] == None:
                continue
            if isinstance(record[field], list):
                length = len(record[field])
                for i in range(0, length):
                    tmp += record[field][i].split(',')
            elif isinstance(record[field], str) and len(record[field]) > 0:
                tmp += record[field].split(',')
            if len(tmp) > 0:
                record[field] = []
                for d in tmp:
                    if len(d) > 0:
                        record[field].append(d)

        mobtyper_results.append(record)



    for i in range(0, len(mobtyper_results)):
        record = mobtyper_results[i]
        sample_id = record['sample_id']
        if isinstance(record['sample_id'], list):
            sample_id = record['sample_id'][0]
        if sample_id in id_mapping:
            original_id = id_mapping[sample_id]
            record['sample_id'] = original_id
        bio_markers = sort_biomarkers({0: {'types': record['rep_type(s)'], 'acs': record['rep_type_accession(s)']},
                                       1: {'types': record['relaxase_type(s)'],
                                           'acs': record['relaxase_type_accession(s)']},
                                       2: {'types': record['mpf_type'], 'acs': record['mpf_type_accession(s)']},
                                       3: {'types': record['orit_type(s)'], 'acs': record['orit_accession(s)']}, })

        record['rep_type(s)'] = bio_markers[0]['types']
        record['rep_type_accession(s)'] = bio_markers[0]['acs']
        record['relaxase_type(s)'] = bio_markers[1]['types']
        record['relaxase_type_accession(s)'] = bio_markers[1]['acs']
        record['mpf_type'] = bio_markers[2]['types']
        record['mpf_type_accession(s)'] = bio_markers[2]['acs']
        record['orit_type(s)'] = bio_markers[3]['types']
        record['orit_accession(s)'] = bio_markers[3]['acs']

        if (isinstance(record['mash_neighbor_distance'], float) or isinstance(record['mash_neighbor_distance'],
                                                                              int)) and record[
            'mash_neighbor_distance'] <= primary_distance:
            mob_cluster_id = record['primary_cluster_id']
        else:
            mob_cluster_id = None

        # Patches that sometimes results are concatonated into strings if contigs are merged into a single results
        if isinstance(record['rep_type(s)'], list):
            record['rep_type(s)'] = ",".join(record['rep_type(s)'])
        if isinstance(record['relaxase_type_accession(s)'], list):
            record['relaxase_type_accession(s)'] = ",".join(record['relaxase_type_accession(s)'])

        host_range = hostrange(record['rep_type(s)'].split(','), record['relaxase_type_accession(s)'].split(','),
                               mob_cluster_id, ncbi, lit, ETE3DBTAXAFILE, database_dir)

        for field in host_range:
            record[field] = host_range[field]

        if isinstance(record['mpf_type'], list):
            record['mpf_type'] = determine_mpf_type(record['mpf_type'])
        elif isinstance(record['mpf_type'], str):
            record['mpf_type'] = determine_mpf_type(record['mpf_type'].split(','))

        for field in record:
            if isinstance(record[field], list):
                record[field] = ",".join(record[field])

        record['predicted_mobility'] = 'non-mobilizable'
        if len(record['relaxase_type(s)']) > 0 and len(record['mpf_type']):
            record['predicted_mobility'] = 'conjugative'
        elif len(record['relaxase_type(s)']) > 0 or len(record['orit_type(s)']) > 0:
            record['predicted_mobility'] = 'mobilizable'

        mobtyper_results[i] = record

    writeReport(mobtyper_results, MOB_TYPER_REPORT_HEADER, report_file)

    # Peform MGE detection
    if mge_report_file is not None:
        mge_results = blast_mge(fixed_fasta, repetitive_mask_file, tmp_dir, min_length,
                                logging, min_rpp_ident, min_rpp_cov, min_rpp_evalue, num_threads)

        tmp = {}
        for contig_id in mge_results:
            if contig_id in id_mapping:
                label = id_mapping[contig_id]
            else:
                continue
            tmp[label] = mge_results[contig_id]

        mge_results = tmp
        contig_memberships = {'chromosome': {}, 'plasmid': {}}
        for i in range(0, len(mobtyper_results)):
            primary_cluster_id = mobtyper_results[i]['primary_cluster_id']
            if not primary_cluster_id in contig_memberships['plasmid']:
                contig_memberships['plasmid'][primary_cluster_id] = {}
            contig_id = mobtyper_results[i]['sample_id']
            mobtyper_results[i]['molecule_type'] = 'plasmid'
            mobtyper_results[i]['contig_id'] = contig_id
            contig_memberships['plasmid'][primary_cluster_id][contig_id] = mobtyper_results[i]

        if len(mge_results) > 0:
            writeMGEresults(contig_memberships, mge_results, mge_report_file)
            logger.info("MOB-typer MGE results written to {}".format(mge_report_file))
        else:
            logger.info("No MOB-typer MGE to write")

    if biomarker_report_file is not None:
            
        biomarker_params = {
            'oriT': {
                'file':orit_blast_results,
                'min_length': 80,
                'min_cov':min_rep_cov,
                'min_hsp_cov': 15,
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
        biomarker_df.to_csv(biomarker_report_file,header=True,sep="\t",index=False)


    if not keep_tmp:
        shutil.rmtree(tmp_dir)
    logger.info("MOB-typer completed and results written to {}".format(report_file))


# call main function
if __name__ == '__main__':
    main()
