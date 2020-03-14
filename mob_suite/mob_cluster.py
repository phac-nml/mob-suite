import logging
import os
import shutil
import sys
import re
from argparse import (ArgumentParser)
from mob_suite.version import __version__
import pandas as pd
import scipy
import scipy.cluster.hierarchy as sch
from Bio import SeqIO
import glob
from scipy.cluster.hierarchy import fcluster
from mob_suite.mob_typer import getRepliconContigs
from mob_suite.blast import BlastRunner
from mob_suite.utils import \
    read_fasta_dict, \
    calcFastaStatsIndividual, \
    fix_fasta_header, \
    filterFastaByIDs, \
    run_mob_typer, \
    write_individual_fasta, \
    replicon_blast, \
    check_dependencies,\
    mob_blast, \
    read_sequence_info, \
    dict_from_alt_key

from mob_suite.wrappers import mash

ACS_LETTER_VALUES = {
                    'A':0, 'B':1, 'C':2, 'D':3, 'E':4, 'F':5, 'G':6, 'H':7, 'I':8, 'J':9,
                    'K':10, 'L':11, 'M':12, 'N':13, 'O':14, 'P':15, 'Q':16, 'R':17, 'S':18,
                    'T':19, 'U':20, 'V':21, 'W':22, 'X':23, 'Y':24, 'Z':25
                     }

ACS_VALUES_TO_LETTERS = {
                        0:'A', 1:'B', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I', 9:'J',
                        10:'K', 11:'L', 12:'M', 13:'N', 14:'O', 15:'P', 16:'Q', 17:'R', 18:'S',
                        19:'T', 20:'U', 21:'V', 22:'W', 23:'X', 24:'Y', 25:'Z'
                        }


MAX_ACS_VALUE = 650999
ACS_FORMAT_VALUES = [26000,1000,1]

MOB_CLUSTER_INFO_HEADER = ['id','size','gc_content','md5','organism','primary_cluster_id','primary_dist','secondary_cluster_id','secondary_dist',"rep_type(s)","rep_type_accession(s)","relaxase_type(s)","relaxase_type_accession(s)"]
MOB_CLUSTER_MINIMAL_INFO = ['id','organism']


LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'

def init_console_logger(lvl):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    report_lvl = logging_levels[lvl]

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)
    return logging

def parse_args():
    "Parse the input arguments, use '-h' for help"
    default_database_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'databases')
    parser = ArgumentParser(description="MOB-Cluster: Generate and update existing plasmid clusters' version: {}".format(__version__))
    parser.add_argument('-m','--mode', type=str, required=True, help='Build: Create a new database from scratch, Update: Update an existing database with one or more sequences')
    parser.add_argument('-o','--outdir', type=str, required=True, help='Output Directory to put results')
    parser.add_argument('-i','--infile', type=str, required=True, help='Input fasta file of one or more closed plasmids to process')
    parser.add_argument('-t', '--sequence_meta', type=str, required=True,
                        help='TSV file for new sequences with the fields "id, header"')
    parser.add_argument('-c','--ref_cluster_file', type=str, required=False, help='Reference mob-cluster file')
    parser.add_argument('-r','--ref_fasta_file', type=str, required=False, help='Reference mob-cluster fasta file')
    parser.add_argument('--ref_mash_db', type=str, required=False, help='Reference mob-cluster mash sketch file')
    parser.add_argument('--num_threads', type=int, required=False, help='Number of threads to be used', default=1)
    parser.add_argument('--primary_cluster_dist', type=int, required=False, help='Mash distance for assigning primary cluster id 0 - 1', default=0.06)
    parser.add_argument('--secondary_cluster_dist', type=int, required=False, help='Mash distance for assigning primary cluster id 0 - 1',
                        default=0.025)
    parser.add_argument('--min_rep_evalue', type=str, required=False,
                        help='Minimum evalue threshold for replicon blastn',
                        default=0.00001)

    parser.add_argument('--min_mob_evalue', type=str, required=False,
                        help='Minimum evalue threshold for relaxase tblastn',
                        default=0.00001)

    parser.add_argument('--min_rep_ident', type=int, required=False, help='Minimum sequence identity for replicons',
                        default=80)

    parser.add_argument('--min_mob_ident', type=int, required=False, help='Minimum sequence identity for relaxases',
                        default=80)
    parser.add_argument('--min_rep_cov', type=int, required=False,
                        help='Minimum percentage coverage of replicon query by input assembly',
                        default=80)

    parser.add_argument('--min_mob_cov', type=int, required=False,
                        help='Minimum percentage coverage of relaxase query by input assembly',
                        default=80)
    parser.add_argument('--min_overlap', type=int, required=False,
                        help='Minimum overlap of fragments',
                        default=10)
    parser.add_argument('--plasmid_replicons', type=str, required=False, help='Fasta of plasmid replicons',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/rep.dna.fas'))
    parser.add_argument('--plasmid_mob', type=str, required=False, help='Fasta of plasmid relaxases',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/mob.proteins.faa'))
    parser.add_argument('-w','--overwrite',  required=False, help='Overwrite the MOB-suite databases with results', action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-d','--database_directory',default=default_database_dir,
                        required=False,
                        help='Directory you want to use for your databases. If the databases are not already '
                             'downloaded, they will be downloaded automatically. Defaults to {}. '
                             'If you change this from the default, will override --plasmid_mash_db, '
                             '--plasmid_replicons, --plasmid_mob, --plasmid_mpf, and '
                             '--plasmid_orit'.format(default_database_dir))
    return parser.parse_args()



'''
    Input: Accession string
    Return: True if string is in valid format : [A-Z][A-Z][0-9][0-9][0-9] else False
    
'''
def validate_acs_format(accession):
    if not isinstance(accession,str):
        logging.warning('Error provided accession number is not a string: {}'.format(accession))
        return False
    if not re.search("[A-Z][A-Z]\d\d\d",accession):
        logging.warning('Error provided accession number is not in the correct format ->[A-Z][A-Z][0-9][0-9][0-9]<-: {}'.format(accession))
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
        logging.error('Error provided a number greater than what the existing ACS format can accommodate: {}'.format(number))
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
    data = pd.read_csv(file, sep='\t', header=0,names=MOB_CLUSTER_MINIMAL_INFO,index_col=0)
    sequences = dict()
    header = list(data.head())

    for index, row in data.iterrows():
        if not index in sequences:
            sequences[index] = {}

        for i in range(0,len(header)):
            if row[header[i]] == 'nan' :
                sequences[index][header[i]] = ''
            else:
                sequences[index][header[i]] = row[header[i]]

    return sequences


'''
    Input: fast file, mash reference sketch file, output file for mash distances
    Output: Mash distance results

'''

def calcDistances(input_fasta,ref_sketch_db,output_file):
    m = mash()
    mash_results = dict()
    with open(output_file, 'w',encoding="utf-8") as oh:
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


def write_clusters(file,cluster_assignments,header):
    with open(file, 'w', encoding="utf-8") as outfile:
        for id in cluster_assignments:
            outfile.write(str(id) + "\t" + "\t".join(cluster_assignments[id]))
        outfile.close()


def build_cluster_db(distance_matrix_file,distances):
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



def writeClusterAssignments(output_file,header,cluster_assignmnets):
    with open(output_file, 'w', encoding="utf-8") as out:
        out.write("\t".join(map(str,header)) + "\n")
        for id in cluster_assignmnets:
            line = [id]
            for h in header:
                if h == 'id':
                    continue
                if h in cluster_assignmnets[id]:
                    if cluster_assignmnets[id][h] == 'nan':
                        cluster_assignmnets[id][h] = ''

                    line.append(str(cluster_assignmnets[id][h]))
                else:
                    line.append('')
            out.write("{}\n".format("\t".join(line)))
        out.close()

def appendFasta(new_seq_fasta,out_fasta):
    fh = open(out_fasta,'a')
    for record in SeqIO.parse(new_seq_fasta, "fasta"):
        fh.write(">{}\n{}\n".format(record.id,record.seq))
    fh.close()




def updateFastaFile(in_fasta_file,out_fasta_file,cluster_assignments):
    out = open(out_fasta_file,'w', encoding="utf-8")
    with open(in_fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            row = str(record.id).split('|')
            id = row[0]

            if not id in cluster_assignments:
                continue

            out.write(">{}\n{}\n".format(id,record.seq))
        handle.close()
    out.close()

def selectCluster(clust_assignments,column):
    out = dict()
    for id in clust_assignments:
        row = str(id).split('|')
        out[row[0]] = clust_assignments[id][column]
    return out

def getMashDistances(mash_results_file,max_thresh=1):
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


def update_existing_db(new_seq_info,existing_seq_info,unique_new_seq,ref_mashdb,tmp_dir,primary_dist,secondary_dist,num_threads=1):

    existing_acs = []
    for id in existing_seq_info:
        print(existing_seq_info[id])
        p = existing_seq_info[id]['primary_cluster_id']
        is_valid = validate_acs_format(p)

        if is_valid:
            p = int(acs_to_int(p))
        existing_acs.append(p)

        p = existing_seq_info[id]['secondary_cluster_id']
        is_valid = validate_acs_format(p)

        if is_valid:
            p = int(acs_to_int(p))
        existing_acs.append(p)

    highest_id = max(existing_acs) + 1

    md5_cluster_assignments = {}


    for id in unique_new_seq:

        mash_results_file  = os.path.join(tmp_dir, "{}.mash.txt".format(id))
        fasta_file = os.path.join(tmp_dir,"{}.fasta".format(id))

        # Run mash to get distances to the reference
        out_dir = os.path.dirname(fasta_file)
        mashObj = mash()
        mashObj.mashsketch(fasta_file, output_path=fasta_file, sketch_ind=True, num_threads=1, kmer_size=21,
                           sketch_size=1000)
        mashfile_handle = open(mash_results_file, 'w', encoding="utf-8")
        mashObj.run_mash(ref_mashdb, fasta_file + '.msh', mashfile_handle, table=False, num_threads=num_threads)
        mashfile_handle.close()
        query_mash_distances = getMashDistances(mash_results_file,primary_dist)

        if len(query_mash_distances) == 0:
            primary_query_acs = highest_id
            highest_id += 1
            secondary_query_acs = highest_id
            highest_id += 1

        for query_id in query_mash_distances:
            min_distance = min(query_mash_distances[query_id].values())
            min_dist_key = min(query_mash_distances[query_id], key=query_mash_distances[query_id].get)

            if min_distance > primary_dist:
                primary_query_acs = highest_id
                highest_id += 1
                secondary_query_acs = highest_id
                highest_id+=1
            elif min_distance < secondary_dist:
                primary_query_acs = existing_seq_info[min_dist_key]['primary_cluster_id']
                secondary_query_acs = existing_seq_info[min_dist_key]['secondary_cluster_id']
            elif min_distance > secondary_dist:
                primary_query_acs = existing_seq_info[min_dist_key]['primary_cluster_id']
                secondary_query_acs = highest_id
                highest_id+=1
        md5 = new_seq_info[id]['md5']
        if not md5 in md5_cluster_assignments:
            md5_cluster_assignments[md5]  = {}

        new_seq_info[id]['primary_cluster_id'] = int_to_acs(primary_query_acs)
        new_seq_info[id]['secondary_cluster_id'] = int_to_acs(secondary_query_acs)

        md5_cluster_assignments[md5]['primary_cluster_id'] = new_seq_info[id]['primary_cluster_id']
        md5_cluster_assignments[md5]['secondary_cluster_id'] = new_seq_info[id]['secondary_cluster_id']

    for id in new_seq_info:
        md5 = new_seq_info[id]['md5']
        if md5 in md5_cluster_assignments:
            p = md5_cluster_assignments[md5]['primary_cluster_id']
            s = md5_cluster_assignments[md5]['secondary_cluster_id']
            new_seq_info[id]['primary_cluster_id'] = p
            new_seq_info[id]['secondary_cluster_id'] = s
    return new_seq_info


'''
    Input: Two dictionaries
    Output: the key:value pairs which are not present in both dictionaries
'''
def find_different_dict_keys(first_dict,second_dict):

    set1 = first_dict.keys()
    set2 = second_dict.keys()

    return(set1 ^ set2)

'''
    Input: Dictionary indexed by md5 hash for sequences for new sequences and existing sequences.
    Output: Returns one sequence per md5 key which exists in the new set
'''
def getUniqSeqKeys(md5_new_seq_lookup,md5_mobcluster_seq_lookup):

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


def combine_seqStats_with_userMeta(seqStats,userMeta):
    if len(seqStats.keys()) != len(userMeta):
        logging.error('Error cannot continue due to difference number of sequences in fasta and metadata files: Num Fasta = {}, Num Meta = {}'.format(len(seqStats),len(userMeta)))
        return {}

    missing_keys = find_different_dict_keys(seqStats,userMeta)

    if len(missing_keys) > 0:
        logging.error('Error cannot continue due to different ids in fasta and metadata files, keys are: '.format(missing_keys))
        return {}

    for id  in userMeta :
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

    primary_keys = list(set(primary_keys ))
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

    return  clust_assignments





def main():
    default_database_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'databases')
    args = parse_args()
    logging = init_console_logger(3)
    logging.info('Running Mob-Suite Clustering toolkit v. {}'.format(__version__))
    logging.info('Processing fasta file {}'.format(args.infile))
    logging.info('Analysis directory {}'.format(args.outdir))

    check_dependencies(logging)

    database_dir = os.path.abspath(args.database_directory)
    if database_dir == default_database_dir:
        mob_ref = args.plasmid_mob
        replicon_ref = args.plasmid_replicons
    else:
        mob_ref = os.path.join(database_dir, 'mob.proteins.faa')
        replicon_ref = os.path.join(database_dir, 'rep.dna.fas')

    min_rep_ident = float(args.min_rep_ident)
    min_mob_ident = float(args.min_mob_ident)
    min_rep_evalue = float(args.min_rep_evalue)
    min_mob_evalue = float(args.min_mob_evalue)
    min_rep_cov = float(args.min_rep_cov)
    min_mob_cov = float(args.min_mob_cov)

    input_fasta = args.infile
    if not os.path.isfile(input_fasta):
        logging.error('Error, input fasta specified does not exist: {}'.format(input_fasta ))
        sys.exit()

    user_sequence_metadata_file = args.sequence_meta
    if not os.path.isfile(input_fasta):
        logging.error('Error, input metadata file specified does not exist: {}'.format(user_sequence_metadata_file ))
        sys.exit()

    mode = str(args.mode).lower()
    if mode not in ('update','build'):
        logging.error('Error you have not entered a valid mode of build or update, you entered: {}'.format(mode))
        sys.exit()

    out_dir = args.outdir
    num_threads = args.num_threads
    if not os.path.isdir(out_dir):
        logging.info('Creating directory {}'.format(args.outdir))
        os.mkdir(out_dir, 0o755)
    tmp_dir = os.path.join(out_dir, '__tmp')
    if not os.path.isdir(tmp_dir):
        logging.info('Creating directory {}'.format(args.outdir))
        os.mkdir(tmp_dir, 0o755)

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

    fixed_fasta = os.path.join(tmp_dir,"fixed.fasta")

    fix_fasta_header(input_fasta, fixed_fasta)
    fastaSeqStats = calcFastaStatsIndividual(fixed_fasta)
    new_sequence_metadata = read_user_new_sequence_info(user_sequence_metadata_file)
    new_seq_info = combine_seqStats_with_userMeta(fastaSeqStats,new_sequence_metadata)
    for id in new_seq_info:
        new_seq_info[id]['primary_dist'] = primary_distance
        new_seq_info[id]['secondary_dist'] = secondary_distance


    if len(new_seq_info) == 0:
        logging.error('Error something has gone wrong with update sequence or metadata files, check log for info')
        sys.exit()


    md5_new_seq_lookup = dict_from_alt_key(new_seq_info,'md5')

    unique_new_sequence_file = os.path.join(tmp_dir, "unique_seqs.fasta")
    unique_seq_keys = getUniqSeqKeys(md5_new_seq_lookup, {})
    filterFastaByIDs(fixed_fasta, unique_new_sequence_file, unique_seq_keys)

    write_individual_fasta(unique_new_sequence_file,tmp_dir)

    #Get Replicon and relaxase typing for all new sequences
    for seq_id in unique_seq_keys:
        fasta_file = os.path.join(tmp_dir,"{}.fasta".format(seq_id))
        blast_results = os.path.join(tmp_dir,"{}.blast".format(seq_id))
        replicon_contigs = (getRepliconContigs(
            replicon_blast(

                        replicon_ref,
                fasta_file,
                           min_rep_ident,
                           min_rep_cov,
                           min_rep_evalue,
                           tmp_dir,
                            blast_results,
                           num_threads=num_threads)))
        found_replicons = dict()

        for contig_id in replicon_contigs:
            for hit in replicon_contigs[contig_id]:
                acs, type = hit.split('|')
                found_replicons[acs] = type

        #Apply result to all identical sequences
        md5 = new_seq_info[seq_id]['md5']
        for sid in md5_new_seq_lookup[md5]:

            new_seq_info[sid]['rep_type(s)'] = ",".join(list(found_replicons.values()))
            new_seq_info[sid]['rep_type_accession(s)'] = ",".join(list(found_replicons.keys()))

        mob_contigs = getRepliconContigs(
            mob_blast(mob_ref, fasta_file, min_mob_ident, min_mob_cov, min_mob_evalue, tmp_dir, blast_results,
                      num_threads=num_threads))
        found_mob = dict()
        for contig_id in mob_contigs:
            for hit in mob_contigs[contig_id]:
                acs, type = hit.split('|')
                found_mob[acs] = type
        for sid in md5_new_seq_lookup[md5]:
            new_seq_info[sid]['relaxase_type(s)'] = ",".join(list(found_mob.values()))
            new_seq_info[sid]['relaxase_type_accession(s)'] = ",".join(list(found_mob.keys()))



    tmp_cluster_file = os.path.join(out_dir, 'clusters.txt')
    tmp_ref_fasta_file = os.path.join(tmp_dir, 'references_tmp.fasta')
    update_fasta = os.path.join(out_dir, 'references_updated.fasta')

    if mode == 'update':
        if args.ref_cluster_file is None:
            logging.error('Reference fasta file must be specified, please check help for parameter reference')
            sys.exit()

        ref_fasta = args.ref_fasta_file
        if not os.path.isfile(ref_fasta ):
            logging.error('Reference fasta file specified does not exist: {}'.format(ref_fasta))
            sys.exit()

        if args.ref_cluster_file is None:
            logging.error('Reference cluster file must be specified, please check help for parameter reference')
            sys.exit()

        ref_cluster_file = args.ref_cluster_file
        if not os.path.isfile(ref_cluster_file):
            logging.error('Reference cluster file specified does not exist: {}'.format(ref_cluster_file))
            sys.exit()

        if args.ref_mash_db is None:
            logging.error('Reference mash sketch file must be specified, please check help for parameter reference')
            sys.exit()

        ref_mash_db = args.ref_mash_db
        if not os.path.isfile(ref_mash_db):
            logging.error('Reference mash file specified does not exist: {}'.format(ref_mash_db))
            sys.exit()

        mob_cluster_seq_info = read_sequence_info(ref_cluster_file)
        md5_mobcluster_seq_lookup = dict_from_alt_key(mob_cluster_seq_info , 'md5')
        unique_new_sequence_file = os.path.join(tmp_dir,"unique_seqs.fasta")
        unique_seq_keys = getUniqSeqKeys(md5_new_seq_lookup,md5_mobcluster_seq_lookup)
        logging.info('Writting unique sequences to file: {}'.format(unique_new_sequence_file))
        filterFastaByIDs(fixed_fasta, unique_new_sequence_file, unique_seq_keys)

        logging.info('Running mob-cluster in update mode with input file: {}'.format(input_fasta))
        logging.info('Running mob-cluster in update mode with output directory: {}'.format(out_dir))
        logging.info('Running mob-cluster in update mode on reference fasta file: {}'.format(ref_fasta))
        logging.info('Reading previous cluster reference assignments from : {}'.format(ref_cluster_file))

        shutil.copy(ref_cluster_file, tmp_cluster_file)
        shutil.copy(ref_fasta, tmp_ref_fasta_file)
        logging.info('Creating new cluster assignments')
        new_seq_info = update_existing_db(new_seq_info,
                                          mob_cluster_seq_info,
                                          unique_seq_keys,
                                          ref_mash_db,
                                          tmp_dir,
                                          primary_distance,
                                          secondary_distance,
                                          num_threads)

        cluster_assignments = {**mob_cluster_seq_info, **new_seq_info}
        logging.info('Writting cluster assignments to : {}'.format(tmp_cluster_file))
        writeClusterAssignments(tmp_cluster_file, MOB_CLUSTER_INFO_HEADER, cluster_assignments)

        appendFasta(fixed_fasta, tmp_ref_fasta_file)

        shutil.copy(tmp_ref_fasta_file, os.path.join(out_dir,update_fasta))
        mash_db_file = "{}.msh".format(update_fasta)
        mObj = mash()
        mObj.mashsketch(update_fasta, mash_db_file, num_threads=num_threads)
        blast_runner = BlastRunner(update_fasta, '')
        blast_runner.makeblastdb(update_fasta, 'nucl')

        if args.overwrite:
            shutil.copy(update_fasta,ref_fasta)
            shutil.copy(tmp_cluster_file, ref_cluster_file)
            shutil.copy(mash_db_file,ref_mash_db)
            blast_runner = BlastRunner(ref_fasta, '')
            blast_runner.makeblastdb(ref_fasta, 'nucl')




    else:
        mashObj = mash()
        mashObj.mashsketch(input_fasta,input_fasta+".msh",num_threads=num_threads)
        distance_matrix_file = os.path.join(tmp_dir,'mash_dist_matrix.txt')
        mashfile_handle = open(distance_matrix_file,'w',encoding="utf-8")

        mashObj.run_mash(input_fasta+'.msh', input_fasta+'.msh', mashfile_handle,table=True,num_threads=num_threads)
        clust_assignments = build_cluster_db(distance_matrix_file, (primary_distance, secondary_distance))
        cluster_acs = convert_num_to_acs(clust_assignments)

        for id in cluster_acs:
            primary_key = cluster_acs[id][0]
            secondary_key = cluster_acs[id][1]
            new_seq_info[id]['primary_cluster_id'] = primary_key
            new_seq_info[id]['primary_dist'] = primary_distance
            new_seq_info[id]['secondary_cluster_id'] = secondary_key
            new_seq_info[id]['secondary_dist'] = secondary_distance

        writeClusterAssignments(tmp_cluster_file, MOB_CLUSTER_INFO_HEADER, new_seq_info)
        shutil.copy(fixed_fasta, update_fasta)
    shutil.rmtree(tmp_dir)










'''

def main():
    args = parse_args()
    logging = init_console_logger(3)
    logging.info('Running Mob-Suite Clustering toolkit v. {}'.format(__version__))
    logging.info('Processing fasta file {}'.format(args.infile))
    logging.info('Analysis directory {}'.format(args.outdir))

    input_fasta = args.infile
    if not os.path.isfile(input_fasta):
        logging.error('Error, input fasta specified does not exist: {}'.format(input_fasta ))
        sys.exit()

    out_dir = args.outdir
    num_threads = args.num_threads
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir, 0o755)
    tmp_dir = os.path.join(out_dir, '__tmp')
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir, 0o755)

    mode = str(args.mode).lower()

    fixed_fasta = os.path.join(tmp_dir,"fixed.fasta")

    fix_fasta_header(input_fasta, fixed_fasta)
    fastaSeqStats = calcFastaStatsIndividual(fixed_fasta)

    if mode not in ('update','build'):
        logging.error('Error you have not entered a valid mode of build or update, you entered: {}'.format(mode))
        print(('Error you have not entered a valid mode of build or update, you entered: {}'.format(mode)))
        sys.exit()

    header = ('id', 0.06, 0.025)
    tmp_cluster_file = os.path.join(out_dir, 'clusters.txt')
    tmp_ref_fasta_file = os.path.join(tmp_dir, 'references_tmp.fasta')
    update_fasta = os.path.join(out_dir, 'references_updated.fasta')

    if mode == 'update':
        if args.ref_cluster_file is None:
            logging.error('Reference fasta file must be specified, please check help for parameter reference')
            sys.exit()

        ref_fasta = args.ref_fasta_file

        if not os.path.isfile(ref_fasta ):
            logging.error('Reference fasta file specified does not exist: {}'.format(ref_fasta))
            sys.exit()

        if args.ref_cluster_file is None:
            logging.error('Reference cluster file must be specified, please check help for parameter reference')
            sys.exit()

        ref_cluster_file = args.ref_cluster_file

        if not os.path.isfile(ref_cluster_file):
            logging.error('Reference cluster file specified does not exist: {}'.format(ref_cluster_file))
            sys.exit()

        if args.ref_mash_db is None:
            logging.error('Reference mash sketch file must be specified, please check help for parameter reference')
            sys.exit()

        ref_mash_db = args.ref_mash_db
        if not os.path.isfile(ref_mash_db):
            logging.error('Reference mash file specified does not exist: {}'.format(ref_mash_db))
            sys.exit()

        logging.info('Running mob-cluster in update mode with input file: {}'.format(input_fasta))
        logging.info('Running mob-cluster in update mode with output directory: {}'.format(out_dir))
        logging.info('Running mob-cluster in update mode on reference fasta file: {}'.format(ref_fasta))
        logging.info('Reading previous cluster reference assignments from : {}'.format(ref_cluster_file))

        shutil.copy(ref_cluster_file, tmp_cluster_file)
        shutil.copy(ref_fasta, tmp_ref_fasta_file)
        update_existing(input_fasta, tmp_dir, ref_mash_db, tmp_cluster_file, header, tmp_ref_fasta_file, update_fasta)

        if args.overwrite:
            shutil.move(update_fasta,ref_fasta)
            shutil.move(tmp_cluster_file,ref_cluster_file)
            mash_db_file = "{}.msh".format(input_fasta)
            mObj = mash()
            mObj.mashsketch(input_fasta, mash_db_file, num_threads=num_threads)
            blast_runner = BlastRunner(ref_fasta, '')
            blast_runner.makeblastdb(ref_fasta, 'nucl')
    else:
        mashObj = mash()
        mashObj.mashsketch(input_fasta,input_fasta+".msh",num_threads=num_threads)
        distance_matrix_file = os.path.join(tmp_dir,'mash_dist_matrix.txt')
        mashfile_handle = open(distance_matrix_file,'w',encoding="utf-8")

        mashObj.run_mash(input_fasta+'.msh', input_fasta+'.msh', mashfile_handle,table=True,num_threads=num_threads)
        clust_assignments = build_cluster_db(distance_matrix_file, (0.06, 0.025))
        writeClusterAssignments(tmp_cluster_file, header, clust_assignments)
        clust_dict = selectCluster(clust_assignments, 0)
        shutil.copy(input_fasta, tmp_ref_fasta_file)
        updateFastaFile(tmp_ref_fasta_file ,update_fasta, clust_dict)




'''


# call main function
if __name__ == '__main__':
    main()


