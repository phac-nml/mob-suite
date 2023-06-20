import logging, os, shutil, sys, re, scipy
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter,RawDescriptionHelpFormatter)
from mob_suite.version import __version__
import pandas as pd
import scipy.cluster.hierarchy as sch
from Bio import SeqIO

from scipy.cluster.hierarchy import fcluster
from mob_suite.blast import BlastRunner
from mob_suite.utils import \
    check_dependencies, \
    read_sequence_info, \
    read_file_to_dict, \
    read_fasta_dict, \
    NamesToTaxIDs

from mob_suite.wrappers import mash

from mob_suite.constants import LOG_FORMAT, ACS_LETTER_VALUES, ACS_FORMAT_VALUES, ACS_VALUES_TO_LETTERS, MAX_ACS_VALUE, \
    MOB_TYPER_REPORT_HEADER, MOB_CLUSTER_INFO_HEADER, default_database_dir


def init_console_logger(lvl):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    report_lvl = logging_levels[lvl]

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)
    return logging


def parse_args():
    "Parse the input arguments, use '-h' for help"
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass
    parser = ArgumentParser(
        description="MOB-Cluster: Generate and update existing plasmid clusters' version: {}".format(__version__), formatter_class=CustomFormatter)
    parser.add_argument('-m', '--mode', type=str, required=True,
                        help='Build: Create a new database from scratch, Update: Update an existing database with one or more sequences')
    parser.add_argument('-f', '--infile', type=str, required=True,
                        help='Fasta file of sequences to cluster')
    parser.add_argument('-p', '--mob_typer_file', type=str, required=True,
                        help='MOB-typer report file for new sequences')
    parser.add_argument('-t', '--taxonomy', type=str, required=True,
                        help='TSV file for new sequences with the fields "id, organism"')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output Directory to put results')
    parser.add_argument('-c', '--ref_cluster_file', type=str, required=False,
                        help='Existing MOB-cluster file to add the new sequences to')
    parser.add_argument('-r', '--ref_fasta_file', type=str, required=False,
                        help='Existing MOB-cluster fasta file of sequences contained in the MOB-cluster file')
    parser.add_argument('-d', '--database_directory',
                        default=default_database_dir,
                        required=False,
                        help='Directory you want to use for your databases. If the databases are not already '
                             'downloaded, they will be downloaded automatically. Defaults to {}'.format(
                            default_database_dir))
    parser.add_argument('--num_threads', type=int, required=False, help='Number of threads to be used', default=1)
    parser.add_argument('--primary_cluster_dist', type=float, required=False,
                        help='Mash distance for assigning primary cluster id 0 - 1', default=0.06)
    parser.add_argument('--secondary_cluster_dist', type=float, required=False,
                        help='Mash distance for assigning primary cluster id 0 - 1', default=0.025)
    parser.add_argument('--debug', required=False, help='Show debug information', action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    return parser.parse_args()


'''
    Input: Accession string
    Return: True if string is in valid format : [A-Z][A-Z][0-9][0-9][0-9] else False
    
'''


def validate_acs_format(accession):
    if not isinstance(accession, str):
        logging.warning('Error provided accession number is not a string: {}'.format(accession))
        return False
    if not re.search("[A-Z][A-Z]\d\d\d", accession):
        logging.debug(
            'Error provided accession number is not in the correct format ->[A-Z][A-Z][0-9][0-9][0-9]<-: {}'.format(
                accession))
        return False
    else:
        return True


'''
    Input: Accession string ie. AA000
    Return: returns positive numerical integer value of Accession string, else returns -1
    

'''


def acs_to_int(accession):
    is_valid = validate_acs_format(accession)
    if not is_valid:
        logging.error('Error cannot continue due to invalid accession number {}'.format(accession))
        return -1

    return ACS_LETTER_VALUES[accession[0]] * ACS_FORMAT_VALUES[0] + \
           ACS_LETTER_VALUES[accession[1]] * ACS_FORMAT_VALUES[1] + \
           int(accession[2:5]) * ACS_FORMAT_VALUES[2]


'''
    Input: Integer
    Return: True if positive integer from 0 - MAX_ACS_VALUE

'''


def validate_int(number):
    if not isinstance(number, int):
        logging.error('Error provided numerical id is not a valid integer: {}'.format(number))
        return False
    if number < 0:
        logging.error('Error provided a negative number which is not a valid id: {}'.format(number))
        return False
    if number > MAX_ACS_VALUE:
        logging.error(
            'Error provided a number greater than what the existing ACS format can accommodate: {}'.format(number))
        return False
    else:
        return True


'''
    Input: integer of desired accession less than MAX Integer value the accession can accomodate
    Return: returns accession formated string or empty string on error


'''


def int_to_acs(numerical_id):
    is_valid = validate_int(numerical_id)
    if not is_valid:
        logging.error('Error cannot continue due to invalid id number {}'.format(numerical_id))
        return ""

    acs = [0, 0, 0]
    remainder = numerical_id
    for i in range(0, len(ACS_FORMAT_VALUES)):
        x = int(remainder / ACS_FORMAT_VALUES[i])
        remainder = remainder - (x * ACS_FORMAT_VALUES[i])
        acs[i] = x

    return "{}{}{}".format(ACS_VALUES_TO_LETTERS[acs[0]], ACS_VALUES_TO_LETTERS[acs[1]], str(acs[2]).zfill(3))


'''
    Input: Path to TSV file with MOB_CLUSTER_MINIMAL_INFO fields as the header lines
    Output: Dictionary of sequence indexed by sequence identifier
'''


def read_user_new_sequence_info(file):
    if os.path.getsize(file) == 0:
        return dict()
    data = pd.read_csv(file, sep='\t', header=0, names=MOB_TYPER_REPORT_HEADER, index_col=0)
    sequences = dict()
    header = list(data.head())

    for index, row in data.iterrows():
        if not index in sequences:
            sequences[index] = {}

        for i in range(0, len(header)):
            if row[header[i]] == 'nan':
                sequences[index][header[i]] = ''
            else:
                sequences[index][header[i]] = row[header[i]]

    return sequences


'''
    Input: fast file, mash reference sketch file, output file for mash distances
    Output: Mash distance results

'''


def calcDistances(input_fasta, ref_sketch_db, output_file):
    m = mash()
    mash_results = dict()
    with open(output_file, 'w', encoding="utf-8") as oh:
        m.run_mash(ref_sketch_db, input_fasta, oh)
        mash_lines = m.read_mash(output_file)
        for line in mash_lines:
            row = line.strip("\n").split("\t")
            ref_id = row[0]
            query_id = row[1]
            distance = float(row[2])
            if not query_id in mash_results:
                mash_results[query_id] = dict()
            mash_results[query_id][ref_id] = distance
    return mash_results


def write_clusters(file, cluster_assignments, header):
    with open(file, 'w', encoding="utf-8") as outfile:
        for id in cluster_assignments:
            outfile.write(str(id) + "\t" + "\t".join(cluster_assignments[id]))
        outfile.close()


def build_cluster_db(distance_matrix_file, distances):

    data = pd.read_csv(distance_matrix_file, sep='\t', header=0,
                       index_col=0)

    distance_matrix = data.to_numpy()
    condensed_matrix = scipy.spatial.distance.squareform(distance_matrix)
    Z = scipy.cluster.hierarchy.linkage(condensed_matrix, method='complete')

    clust_assignments = dict()

    for dist in distances:
        index = 0
        clusters = fcluster(Z, dist, criterion='distance')
        for id in data.columns.values:
            if not id in clust_assignments:
                clust_assignments[id] = list()
            clust_assignments[id].append(str(clusters[index]))
            index += 1
    return clust_assignments


def writeClusterAssignments(output_file, header, cluster_assignmnets):
    with open(output_file, 'w', encoding="utf-8") as out:
        out.write("\t".join(map(str, header)) + "\n")
        for id in cluster_assignmnets:
            line = []
            for h in header:
                if h == 'id':
                    continue
                if h in cluster_assignmnets[id]:
                    if cluster_assignmnets[id][h] == 'nan':
                        cluster_assignmnets[id][h] = '-'

                    line.append(str(cluster_assignmnets[id][h]))
                else:
                    line.append('')
            out.write("{}\n".format("\t".join(line)))
        out.close()


def appendFasta(new_seq_fasta, out_fasta):
    fh = open(out_fasta, 'a')
    for record in SeqIO.parse(new_seq_fasta, "fasta"):
        fh.write(">{}\n{}\n".format(record.id, record.seq))
    fh.close()


def updateFastaFile(in_fasta_file, out_fasta_file, cluster_assignments):
    out = open(out_fasta_file, 'w', encoding="utf-8")
    with open(in_fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            row = str(record.id).split('|')
            id = row[0]

            if not id in cluster_assignments:
                continue

            out.write(">{}\n{}\n".format(id, record.seq))
        handle.close()
    out.close()


def selectCluster(clust_assignments, column):
    out = dict()
    for id in clust_assignments:
        row = str(id).split('|')
        out[row[0]] = clust_assignments[id][column]
    return out


def getMashDistances(mash_results_file, max_thresh=1):
    mashfile_handle = open(mash_results_file, 'r', encoding="utf-8")
    query_mash_distances = dict()

    for line in mashfile_handle:
        row = line.split('\t')
        query_id = row[1]
        ref_id = row[0].split('|')[0]
        distance = float(row[2])

        if distance > max_thresh:
            continue
        if not query_id in query_mash_distances:
            query_mash_distances[query_id] = dict()

        query_mash_distances[query_id][ref_id] = distance
    return query_mash_distances


def get_field_values(dictionary, field_key):
    values = []
    for key in dictionary:
        if field_key in dictionary[key]:
            values.append(dictionary[key][field_key])
    return values


def update_existing_db(new_seq_info, mob_cluster_seq_info, clust_assignments, primary_distance, secondary_distance,
                       num_threads=1):
    primary_acs = get_field_values(mob_cluster_seq_info, 'primary_cluster_id')
    primary_acs.sort(reverse=True)
    secondary_acs = get_field_values(mob_cluster_seq_info, 'secondary_cluster_id')
    secondary_acs.sort(reverse=True)

    acs = max(acs_to_int(primary_acs[0]), acs_to_int(secondary_acs[0])) + 1

    pri_acs_mapping = {}
    sec_acs_mapping = {}
    for seq_id in clust_assignments:
        p_acs = clust_assignments[seq_id][0]
        s_acs = clust_assignments[seq_id][1]
        if p_acs not in pri_acs_mapping:
            pri_acs_mapping[p_acs] = int_to_acs(acs)
            acs += 1
        if s_acs not in sec_acs_mapping:
            sec_acs_mapping[s_acs] = int_to_acs(acs)
            acs += 1

    for seq_id in new_seq_info:
        if seq_id in mob_cluster_seq_info:
            continue
        data = new_seq_info[seq_id]
        data['primary_dist'] = primary_distance
        data['secondary_dist'] = secondary_distance
        neighbor_id = data['mash_nearest_neighbor']
        dist = data['mash_neighbor_distance']
        if dist > primary_distance or neighbor_id not in mob_cluster_seq_info:
            data['primary_cluster_id'] = pri_acs_mapping[clust_assignments[seq_id][0]]
        else:
            data['primary_cluster_id'] = mob_cluster_seq_info[neighbor_id]['primary_cluster_id']

        if dist > secondary_distance or neighbor_id not in mob_cluster_seq_info:
            data['secondary_cluster_id'] = sec_acs_mapping[clust_assignments[seq_id][1]]
        else:
            data['secondary_cluster_id'] = mob_cluster_seq_info[neighbor_id]['secondary_cluster_id']

        mob_cluster_seq_info[seq_id] = data

    return mob_cluster_seq_info


'''
    Input: Two dictionaries
    Output: the key:value pairs which are not present in both dictionaries
'''


def find_different_dict_keys(first_dict, second_dict):
    set1 = first_dict.keys()
    set2 = second_dict.keys()

    return (set1 ^ set2)


'''
    Input: Dictionary indexed by md5 hash for sequences for new sequences and existing sequences.
    Output: Returns one sequence per md5 key which exists in the new set
'''


def getUniqSeqKeys(md5_new_seq_lookup, md5_mobcluster_seq_lookup):
    unique_keys = []

    for md5 in md5_new_seq_lookup:
        if md5 in md5_mobcluster_seq_lookup:
            continue
        for id in md5_new_seq_lookup[md5]:
            unique_keys.append(id)
            break

    return unique_keys


'''
    Input: Two dictionaries with identical sets of keys with seqStats produced by calcFastaStatsIndividual 
    and userMeta produced by read_user_new_sequence_info
    Output: Combined dictionary with the elements from both dictionaries
'''


def combine_seqStats_with_userMeta(seqStats, userMeta):
    if len(seqStats.keys()) != len(userMeta):
        logging.error(
            'Error cannot continue due to difference number of sequences in fasta and metadata files: Num Fasta = {}, Num Meta = {}'.format(
                len(seqStats), len(userMeta)))
        return {}

    missing_keys = find_different_dict_keys(seqStats, userMeta)

    if len(missing_keys) > 0:
        logging.error(
            'Error cannot continue due to different ids in fasta and metadata files, keys are: '.format(missing_keys))
        return {}

    for id in userMeta:
        values = userMeta[id]
        for e in values:
            seqStats[id][e] = values[e]

    return seqStats


def convert_num_to_acs(clust_assignments):
    primary_keys = []
    secondary_keys = []

    for id in clust_assignments:

        p = clust_assignments[id][0]
        is_valid = validate_acs_format(p)
        if is_valid:
            p = int(acs_to_int(p))
        primary_keys.append(p)
        if len(clust_assignments[id]) > 1:
            s = clust_assignments[id][1]
            is_valid = validate_acs_format(s)
            if is_valid:
                s = int(acs_to_int(s))
            secondary_keys.append(s)
        clust_assignments[id][0] = int(p)
        clust_assignments[id][1] = int(s)

    primary_keys = list(set(primary_keys))
    unique_ints = []

    for p in primary_keys:
        is_valid = validate_acs_format(p)
        if is_valid:
            p = int(acs_to_int(p))
        unique_ints.append(int(p))

    highest_id = max(unique_ints)

    for id in clust_assignments:
        primary_acs = int_to_acs(clust_assignments[id][0])
        if len(clust_assignments[id]) > 1:
            secondary_acs = int_to_acs(clust_assignments[id][1] + highest_id)
        clust_assignments[id][0] = primary_acs
        clust_assignments[id][1] = secondary_acs

    return clust_assignments


def main():
    args = parse_args()
    if args.debug:
        logging = init_console_logger(3)
    else:
        logging = init_console_logger(2)
    logging.info('Running Mob-Suite Clustering toolkit v. {}'.format(__version__))
    logging.info('Processing fasta file {}'.format(args.infile))
    logging.info('Analysis directory {}'.format(args.outdir))

    check_dependencies(logging)

    input_fasta = args.infile
    if not os.path.isfile(input_fasta):
        logging.error('Error, input fasta specified does not exist: {}'.format(input_fasta))
        sys.exit()

    mob_typer_report_file = args.mob_typer_file
    if not os.path.isfile(mob_typer_report_file):
        logging.error('Error, input metadata file specified does not exist: {}'.format(mob_typer_report_file))
        sys.exit()

    mode = str(args.mode).lower()
    if mode not in ('update', 'build'):
        logging.error('Error you have not entered a valid mode of build or update, you entered: {}'.format(mode))
        sys.exit()

    out_dir = args.outdir
    num_threads = args.num_threads

    if not (args.primary_cluster_dist >= 0 and args.primary_cluster_dist <= 1):
        logging.error('Error distance thresholds must be between 0 - 1: {}'.format(args.primary_cluster_dist))
        sys.exit()
    else:
        primary_distance = args.primary_cluster_dist

    if not (args.secondary_cluster_dist >= 0 and args.secondary_cluster_dist <= 1):
        logging.error('Error distance thresholds must be between 0 - 1: {}'.format(args.secondary_cluster_dist))
        sys.exit()
    else:
        secondary_distance = args.secondary_cluster_dist

    if not os.path.isdir(out_dir):
        logging.info('Creating directory {}'.format(args.outdir))
        os.makedirs(out_dir, 0o755)
    tmp_dir = os.path.join(out_dir, '__tmp')
    if not os.path.isdir(tmp_dir):
        logging.info('Creating directory {}'.format(args.outdir))
        os.makedirs(tmp_dir, 0o755)

    taxonomy_file = args.taxonomy
    database_dir = os.path.abspath(args.database_directory)

    records = read_file_to_dict(mob_typer_report_file, MOB_TYPER_REPORT_HEADER, separater="\t")

    seq_ids = []
    new_seq_info = {}
    duplicate_keys = []
    for record in records:
        seq_ids.append(record['sample_id'])
        if not record['sample_id'] in new_seq_info:
            new_seq_info[record['sample_id']] = record
        else:
            duplicate_keys.append(record['sample_id'])

    if len(duplicate_keys) > 0:
        logging.error(
            "Duplicate sequence identifiers in fasta file. Please make every sequence id unique in the input file before using this tool")
        logging.error(
            "Duplicate sequence ids: {}".format(",".join(duplicate_keys)))
        sys.exit()

    record_identifications = read_file_to_dict(taxonomy_file, ['sample_id', 'organism'], separater="\t")
    organisms = []
    for record in record_identifications:
        organism = record['organism']
        if organism == 'unknown' or organism == '' or organism == 'Unknown':
            organism = 'Bacteria'
        organisms.append(organism)
        seq_id = record['sample_id']
        if seq_id in new_seq_info:
            new_seq_info[seq_id]['organism'] = organism

    ETE3DBTAXAFILE = os.path.abspath(database_dir + "/taxa.sqlite")
    taxids = NamesToTaxIDs(organisms, ETE3DBTAXAFILE, database_dir)
    del(organisms)

    for seq_id in new_seq_info:
        organism = new_seq_info[seq_id]['organism']
        if organism in taxids:
            new_seq_info[seq_id]['taxid'] = taxids[organism][0]
        else:
            new_seq_info[seq_id]['taxid'] = 2

    if len(new_seq_info) == 0:
        logging.error('Error no MOB-typer results for sequences. Sequences must be typed with MOB-typer first')
        sys.exit()

    fasta_dict = read_fasta_dict(input_fasta)

    if len(fasta_dict) == 0:
        logging.error('Error no sequences found in input fasta: {}..cannot continue'.format(input_fasta))
        sys.exit()

    key_set_1 = set(seq_ids)
    key_set_2 = set(list(fasta_dict.keys()))

    if len(list(key_set_1 ^ key_set_2)) > 0:
        logging.error(
            'Error MOB-typer results: {} and input fasta: {} do not have the same set of identifiers, these must match in order to proceed'.format(
                mob_typer_report_file, input_fasta))
        logging.error(
            'Keys present in  MOB-typer results: {} and not in input fasta: {} are: {}'.format(mob_typer_report_file,
                                                                                               input_fasta, list(
                    key_set_1 - key_set_2)))
        logging.error(
            'Keys present in  MOB-typer results: {} and not in input fasta: {} are: {}'.format(mob_typer_report_file,
                                                                                               input_fasta, list(
                    key_set_2 - key_set_1)))
        sys.exit()

    tmp_cluster_file = os.path.join(out_dir, 'clusters.txt')
    tmp_ref_fasta_file = os.path.join(tmp_dir, 'references_tmp.fasta')
    update_fasta = os.path.join(out_dir, 'references_updated.fasta')

    # Sketch and calculate distances within update sequences
    if len(fasta_dict) > 1:
        mashObj = mash()
        mashObj.mashsketch(input_fasta, input_fasta + ".msh", num_threads=num_threads)
        distance_matrix_file = os.path.join(tmp_dir, 'mash_dist_matrix.txt')
        mashfile_handle = open(distance_matrix_file, 'w', encoding="utf-8")
        mashfile_handle.write(
            mashObj.run_mash(input_fasta + '.msh', input_fasta + '.msh', table=True, num_threads=num_threads).decode())
        mashfile_handle.close()
        clust_assignments = build_cluster_db(distance_matrix_file, (primary_distance, secondary_distance))
    else:
        seq_id = next(iter(fasta_dict))
        clust_assignments = {seq_id: [0, 1]}

    logging.info('Running MOB-cluster in {} mode'.format(mode))
    if mode == 'update':

        if args.ref_cluster_file is None:
            logging.error('Reference fasta file must be specified, please check help for parameter reference')
            sys.exit()

        ref_fasta = args.ref_fasta_file
        if not os.path.isfile(ref_fasta):
            logging.error('Reference fasta file specified does not exist: {}'.format(ref_fasta))
            sys.exit()

        if args.ref_cluster_file is None:
            logging.error('Reference cluster file must be specified, please check help for parameter reference')
            sys.exit()

        ref_cluster_file = args.ref_cluster_file
        if not os.path.isfile(ref_cluster_file):
            logging.error('Reference cluster file specified does not exist: {}'.format(ref_cluster_file))
            sys.exit()

        mob_cluster_seq_info = read_sequence_info(ref_cluster_file, MOB_CLUSTER_INFO_HEADER)

        logging.info('Running mob-cluster in update mode with input file: {}'.format(input_fasta))
        logging.info('Running mob-cluster in update mode with output directory: {}'.format(out_dir))
        logging.info('Running mob-cluster in update mode on reference fasta file: {}'.format(ref_fasta))
        logging.info('Reading previous cluster reference assignments from : {}'.format(ref_cluster_file))

        shutil.copy(ref_cluster_file, tmp_cluster_file)
        shutil.copy(ref_fasta, tmp_ref_fasta_file)
        logging.info('Creating new cluster assignments')
        new_seq_info = update_existing_db(new_seq_info,
                                          mob_cluster_seq_info,
                                          clust_assignments,
                                          primary_distance,
                                          secondary_distance,
                                          num_threads)

        cluster_assignments = {**mob_cluster_seq_info, **new_seq_info}
        logging.info('Writting cluster assignments to : {}'.format(tmp_cluster_file))
        writeClusterAssignments(tmp_cluster_file, MOB_CLUSTER_INFO_HEADER, cluster_assignments)
        shutil.copy(tmp_ref_fasta_file, update_fasta)


    else:
        cluster_acs = convert_num_to_acs(clust_assignments)
        for id in cluster_acs:
            primary_key = cluster_acs[id][0]
            secondary_key = cluster_acs[id][1]
            new_seq_info[id]['primary_cluster_id'] = primary_key
            new_seq_info[id]['primary_dist'] = primary_distance
            new_seq_info[id]['secondary_cluster_id'] = secondary_key
            new_seq_info[id]['secondary_dist'] = secondary_distance

        writeClusterAssignments(tmp_cluster_file, MOB_CLUSTER_INFO_HEADER, new_seq_info)
        shutil.copy(input_fasta, update_fasta)

    logging.info("Sketching new fasta {}".format(update_fasta))
    mash_db_file = "{}.msh".format(update_fasta)
    mObj = mash()
    mObj.mashsketch(update_fasta, mash_db_file, num_threads=num_threads)
    logging.info("Building blastdb {}".format(update_fasta))
    blast_runner = BlastRunner(update_fasta, '')
    blast_runner.makeblastdb(update_fasta, 'nucl', logging=logging)
    logging.info("Removing temporary directory")
    shutil.rmtree(tmp_dir)
    logging.info("MOB-cluster completed, analysis results written to {}".format(out_dir))

# call main function
if __name__ == '__main__':
    main()
