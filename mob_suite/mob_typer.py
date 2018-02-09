#!/usr/bin/env python
from __future__ import print_function

from collections import OrderedDict
import logging
import os
import shutil
import sys
from argparse import (ArgumentParser, FileType)

from blast import BlastRunner
from blast import BlastReader
from wrappers import circlator
from wrappers import mash
from classes.mcl import mcl
from utils import \
    fixStart, \
    read_fasta_dict, \
    write_fasta_dict, \
    filter_overlaping_records, \
    replicon_blast, \
    mob_blast, \
    getRepliconContigs, \
    fix_fasta_header, \
    getMashBestHit,\
    calcFastaStats

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Mob Typing: Typing and prediction of plasmid mibility')
    parser.add_argument('--outdir', type=str, required=True, help='Output Directory to put results')
    parser.add_argument('--infile', type=str, required=True, help='Input assembly fasta file to process')
    parser.add_argument('--num_threads', type=int, required=False, help='Number of threads to be used', default=1)
    parser.add_argument('--evalue', type=str, required=False, help='Minimum evalue threshold for blast', default=0.00001)
    parser.add_argument('--min_ident', type=str, required=False, help='Minimum sequence identity', default=80)
    parser.add_argument('--min_cov', type=str, required=False, help='Minimum percentage coverage of assembly contig by the plasmid reference database to be considered', default=65)
    parser.add_argument('--keep_tmp', type=str, required=False,
                        help='Do not delete temporary file directory', default=False)
    parser.add_argument('--plasmid_mash_db', type=str, required=False,
                        help='Companion Mash database of reference database',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/ncbi_plasmid_full_seqs.fas.msh'))
    parser.add_argument('--plasmid_replicons', type=str, required=False, help='Fasta of plasmid replicons',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/rep.dna.fas'))
    parser.add_argument('--plasmid_mob', type=str, required=False, help='Fasta of plasmid relaxases',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/mob.proteins.faa'))
    parser.add_argument('--plasmid_mpf', type=str, required=False, help='Fasta of known plasmid mate-pair proteins',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/mpf.proteins.faa'))
    parser.add_argument('--plasmid_orit', type=str, required=False, help='Fasta of known plasmid oriT dna sequences',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/orit.fas'))
    return parser.parse_args()


def determine_mpf_type(hits):
    types = dict()
    for hit in hits:
        type = hits[hit]
        if not type in types:
            types[type] = 0
        types[type]+=1

    return max(types, key=lambda i: types[i])




def main():
    logging.info('Running Mob-typer v. {}'.format('0.1'))
    args = parse_args()

    if not args.outdir:
        logging.info('Error, no output directory specified, please specify one')
        sys.exit()

    if not args.infile:
        logging.info('Error, no fasta specified, please specify one')
        sys.exit()

    if not os.path.isfile(args.infile):
        logging.info('Error, fasta file does not exist')
        sys.exit()

    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir, 0755)

    if not isinstance(args.num_threads,int):
        logging.info('Error number of threads must be an integer, you specified "{}"'.format(args.num_threads))

    #Script arguments
    input_fasta = args.infile
    input_fasta = args.infile
    out_dir = args.outdir
    num_threads = args.num_threads
    evalue = args.evalue
    keep_tmp = args.keep_tmp
    mob_ref = args.plasmid_mob
    mpf_ref = args.plasmid_mpf
    orit_ref = args.plasmid_orit
    mash_db = args.plasmid_mash_db


    tmp_dir = os.path.join(out_dir, '__tmp')
    file_id = os.path.basename(input_fasta)
    fixed_fasta = os.path.join(tmp_dir, 'fixed.input.fasta')
    replicon_ref = args.plasmid_replicons
    replicon_blast_results = os.path.join(tmp_dir, 'replicon_blast_results.txt')
    mob_blast_results = os.path.join(tmp_dir, 'mob_blast_results.txt')
    mpf_blast_results = os.path.join(tmp_dir, 'mpf_blast_results.txt')
    orit_blast_results = os.path.join(tmp_dir, 'orit_blast_results.txt')
    report_file = os.path.join(out_dir, 'mobtyper_' + file_id + '_report.txt')
    mash_file = os.path.join(tmp_dir, 'mash_' + file_id + '.txt')



    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir, 0755)

    fix_fasta_header(input_fasta, fixed_fasta)


    #run individual marker blasts
    logging.info('Running replicon blast on {}'.format(replicon_ref))
    replicon_contigs = getRepliconContigs(replicon_blast(replicon_ref, fixed_fasta, 85, 85, args.evalue, tmp_dir,replicon_blast_results,num_threads=num_threads))
    found_replicons = dict()
    for contig_id in replicon_contigs:
        for hit in replicon_contigs[contig_id]:
            acs, type = hit.split('|')
            found_replicons[acs] = type

    #print(found_replicons)

    logging.info('Running relaxase blast on {}'.format(mob_ref))
    mob_contigs = getRepliconContigs(mob_blast(mob_ref, fixed_fasta, 85, 85, args.evalue, tmp_dir,mob_blast_results,num_threads=num_threads ))
    found_mob = dict()
    for contig_id in mob_contigs:
        for hit in mob_contigs[contig_id]:
            acs,type = hit.split('|')
            found_mob[acs] = type

    #print (found_mob)

    logging.info('Running mpf blast on {}'.format(mob_ref))
    mpf_contigs = getRepliconContigs(mob_blast(mpf_ref, fixed_fasta, 85, 85, args.evalue, tmp_dir, mpf_blast_results,num_threads=num_threads))
    found_mpf = dict()
    for contig_id in mpf_contigs:
        for hit in mpf_contigs[contig_id]:
            acs, type = hit.split('|')
            found_mpf[acs] = type

    #print(found_mpf)

    logging.info('Running orit blast on {}'.format(replicon_ref))
    orit_contigs = getRepliconContigs(replicon_blast(orit_ref, fixed_fasta, 90, 90, args.evalue, tmp_dir,orit_blast_results,num_threads=num_threads))
    found_orit = dict()
    for contig_id in orit_contigs:
        for hit in orit_contigs[contig_id]:
            acs, type = hit.split('|')
            found_orit[acs] = type

    #print(found_orit)


    #Get closest neighbor by mash distance
    m = mash()
    mash_distances = dict()
    mashfile_handle = open(mash_file, 'w')
    m.run_mash(mash_db, fixed_fasta, mashfile_handle)
    mash_results = m.read_mash(mash_file)
    mash_top_hit = getMashBestHit(mash_results)


    results_fh = open(report_file, 'w')
    results_fh.write("file_id\tnum_contigs\ttotal_length\tgc\t"\
                    "rep_type(s)\trep_type_accession(s)\t"\
                    "relaxase_type(s)\trelaxase_type_accession(s)\t"\
                    "mpf_type\tmpf_type_accession(s)\t"\
                    "orit_type(s)\torit_accession(s)\tPredictedMobility"\
                    "mash_nearest_neighbor\tmash_neighbor_distance\tmash_neighbor_cluster\n")

    if len(found_replicons) > 0:
        rep_types = ",".join(found_replicons.values())
        rep_acs = ",".join(found_replicons.keys())
    else:
        rep_types = "-"
        rep_acs = "-"

    if len(found_mob) > 0:
        mob_types = ",".join(found_mob.values())
        mob_acs = ",".join(found_mob.keys())
    else:
        mob_types = "-"
        mob_acs = "-"

    if len(found_mpf) > 0:
        mpf_type = determine_mpf_type(found_mpf)
        mpf_acs = ",".join(found_mpf.keys())
    else:
        mpf_type = "-"
        mpf_acs = "-"

    if len(found_orit) > 0:
        orit_types = ",".join(found_orit.values())
        orit_acs = ",".join(found_orit.keys())
    else:
        orit_types = "-"
        orit_acs = "-"
    stats = calcFastaStats(fixed_fasta)
    predicted_mobility = 'Non-mobilizable'

    if mob_acs != '-' or orit_acs != '-':
        predicted_mobility = 'Mobilizable'

    if mob_acs != '-' and mpf_acs != '-':
        predicted_mobility = 'Conjugative'



    string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(file_id,stats['num_seq'],stats['size'],stats['gc_content'],
                                                                               rep_types,rep_acs,mob_types,
                                                                               mob_acs,mpf_type,mpf_acs,
                                                                               orit_types,orit_acs,predicted_mobility,
                                                                               mash_top_hit['top_hit'],mash_top_hit['mash_hit_score'],mash_top_hit['clustid'])
    results_fh.write(string)
    print(string)

    if not keep_tmp:
        shutil.rmtree(tmp_dir)




# call main function
if __name__ == '__main__':
    main()

