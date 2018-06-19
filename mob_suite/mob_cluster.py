import logging
import os
import shutil
import sys
from argparse import (ArgumentParser)
from mob_suite.version import __version__
import pandas as pd
import scipy
import scipy.cluster.hierarchy as sch
from Bio import SeqIO
from scipy.cluster.hierarchy import fcluster
from mob_suite.blast import BlastRunner
from mob_suite.utils import \
    read_fasta_dict
from mob_suite.wrappers import mash

LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'

def init_console_logger(lvl):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    report_lvl = logging_levels[lvl]

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)
    return logging

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Mob-Suite: Generate and update existing plasmid clusters')
    parser.add_argument('-m','--mode', type=str, required=True, help='Build: Create a new database from scratch, Update: Update an existing database with one or more sequences')
    parser.add_argument('-o','--outdir', type=str, required=True, help='Output Directory to put results')
    parser.add_argument('-i','--infile', type=str, required=True, help='Input fasta file of one or more closed plasmids to process')
    parser.add_argument('--ref_cluster_file', type=str, required=False, help='Reference mob-cluster file')
    parser.add_argument('--ref_fasta_file', type=str, required=False, help='Reference mob-cluster fasta file')
    parser.add_argument('--ref_mash_db', type=str, required=False, help='Reference mob-cluster mash sketch file')
    parser.add_argument('--num_threads', type=int, required=False, help='Number of threads to be used', default=1)
    parser.add_argument('-w','--overwrite',  required=False, help='Overwrite the MOB-suite databases with results', action='store_true')
    return parser.parse_args()

def read_cluster_assignments(file):
    if os.path.getsize(file) == 0:
        return dict()
    data = pd.read_csv(file, sep='\t', header=0,index_col=0)
    clusters = dict()
    for index, row in data.iterrows():
        record = list()
        for r in row:
            record.append(r)
        clusters[index] = record

    return clusters

def calcDistances(input_fasta,ref_sketch_db,output_file):
    m = mash()
    mash_results = dict()
    with open(output_file, 'w') as oh:
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
    with open(file, 'w') as outfile:
        for id in cluster_assignments:
            outfile.write(str(id) + "\t" + "\t".join(cluster_assignments[id]))
        outfile.close()


def build_cluster_db(distance_matrix_file,distances):
    data = pd.read_csv(distance_matrix_file, sep='\t', header=0,
                       index_col=0)
    distance_matrix = data.as_matrix()
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



def add_new_record(fasta_file,ref_mashdb,mash_results_file,ref_cluster_file,distances_thresholds,num_threads):

    ref_clusters = read_cluster_assignments(ref_cluster_file)

    cluster_pointers = list()
    for id in ref_clusters:
        for i in range(0,len(ref_clusters[id])):
            if len(cluster_pointers) == 0:
                cluster_pointers = list(ref_clusters[id])
                continue
            if ref_clusters[id][i] > cluster_pointers[i]:
                cluster_pointers[i] = ref_clusters[id][i]

    for i in range(0,len(cluster_pointers)):
        cluster_pointers[i]+=1

    max_thresh = max(distances_thresholds)
    min_thresh = min(distances_thresholds)

    # Run mash to get distances to the reference
    out_dir = os.path.dirname(fasta_file)
    mashObj = mash()
    mashObj.mashsketch(fasta_file, output_path=fasta_file, sketch_ind=True, num_threads=1, kmer_size=21, sketch_size=1000)
    mashfile_handle = open(mash_results_file,'w')
    mashObj.run_mash(ref_mashdb, fasta_file + '.msh', mashfile_handle,table=False,num_threads=num_threads)
    mashfile_handle.close()

    mashfile_handle = open(mash_results_file, 'r')
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


    for query_id in query_mash_distances:

        min_distance = min(query_mash_distances[query_id].values())
        min_dist_key = min(query_mash_distances[query_id], key=query_mash_distances[query_id].get)
        clusters = list()
        if float(min_distance) > float(max_thresh):
            for i in range(0,len(distances_thresholds)):
                clusters.append(cluster_pointers[i])
                cluster_pointers[i]+=1

        elif float(min_distance) <= float(min_thresh):

            clusters = list(ref_clusters[min_dist_key])

        else:

            clusters = list(ref_clusters[min_dist_key])
            for i in range(0, len(distances_thresholds)):
                dist = distances_thresholds[i]
                if float(min_distance) <= float(dist):
                    continue
                clusters[i] = cluster_pointers[i]
                cluster_pointers[i] += 1
        ref_clusters[query_id] = clusters

    return ref_clusters


def writeClusterAssignments(output_file,header,cluster_assignmnets):
    with open(output_file, 'w') as out:
        out.write("\t".join(map(str,header)) + "\n")
        for id in cluster_assignmnets:
            out.write("{}\t{}".format(id,"\t".join(map(str,cluster_assignmnets[id]))) + "\n")
        out.close()

def updateFastaFile(in_fasta_file,out_fasta_file,cluster_assignments):
    out = open(out_fasta_file,'w')
    with open(in_fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            row = str(record.id).split('|')
            id = row[0]

            if not id in cluster_assignments:
                continue

            out.write(">{}|{}\n{}\n".format(id,cluster_assignments[id],record.seq))
        handle.close()
    out.close()

def selectCluster(clust_assignments,column):
    out = dict()
    for id in clust_assignments:
        row = str(id).split('|')
        out[row[0]] = clust_assignments[id][column]
    return out


def update_existing(input_fasta,tmp_dir,ref_mash_db,tmp_cluster_file,header,tmp_ref_fasta_file,update_fasta,num_threads=1):
    sequences = read_fasta_dict(input_fasta)

    for id in sequences:
        seq = sequences[id]
        tmp_fasta = os.path.join(tmp_dir,id + '_tmp.fasta')
        with open(tmp_fasta , "w") as fh:
            fh.write("\n>{}\n{}\n".format(id,seq))
            fh.close()
        tmp_mash = os.path.join(tmp_dir,id + '_tmp.txt')
        clust_assignments = add_new_record(tmp_fasta, ref_mash_db ,tmp_mash, tmp_cluster_file,(0.05, 0.0001),num_threads)

        writeClusterAssignments(tmp_cluster_file , header, clust_assignments)
        with open(tmp_ref_fasta_file , "a") as fh:
            fh.write("\n>{}\n{}\n".format(id,seq))
            fh.close()

    clust_dict = selectCluster(clust_assignments, 1)

    updateFastaFile(tmp_ref_fasta_file ,update_fasta, clust_dict)

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

    if mode not in ('update','build'):
        logging.error('Error you have not entered a valid mode of build or update, you entered: {}'.format(mode))
        print(('Error you have not entered a valid mode of build or update, you entered: {}'.format(mode)))
        sys.exit()

    header = ('id', 0.05, 0.0001)
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
        mashfile_handle = open(distance_matrix_file,'w')

        mashObj.run_mash(input_fasta+'.msh', input_fasta+'.msh', mashfile_handle,table=True,num_threads=num_threads)
        clust_assignments = build_cluster_db(distance_matrix_file, (0.05, 0.0001))
        writeClusterAssignments(tmp_cluster_file, header, clust_assignments)
        clust_dict = selectCluster(clust_assignments, 1)
        shutil.copy(input_fasta, tmp_ref_fasta_file)
        updateFastaFile(tmp_ref_fasta_file ,update_fasta, clust_dict)







# call main function
if __name__ == '__main__':
    main()


