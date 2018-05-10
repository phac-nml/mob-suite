#!/usr/bin/env python
import logging, os, sys
from argparse import (ArgumentParser, FileType)
from mob_suite.blast import BlastReader


def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Filter overlapping queries')
    parser.add_argument('--infile', type=str, required=True, help='Input file to process')
    parser.add_argument('--outdir', type=str, required=True, help='Output directory')
    parser.add_argument('--min_overlap', type=str, required=False, help='Minimum bp overlap', default=5)
    return parser.parse_args()

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


def fixStart(blast_df):
    for index, row in blast_df.iterrows():
        sstart = blast_df.at[index, 'sstart']
        send = blast_df.at[index, 'send']
        # print "{}\t{}".format(sstart,send)
        if send < sstart:
            temp = sstart
            blast_df.at[index, 'sstart'] = send
            blast_df.at[index, 'send'] = temp
        # print "====>{}\t{}".format(self.blast_df.at[index, 'sstart'], self.blast_df.at[index, 'send'])
        qstart = blast_df.at[index, 'qstart']
        qend = blast_df.at[index, 'qend']
        if qend < qstart:
            temp = qstart
            blast_df.at[index, 'qstart'] = qend
            blast_df.at[index, 'qend'] = temp
    return blast_df

def filter_blast(blast_results_file, min_ident, min_cov, evalue, overlap):
    if os.path.getsize(blast_results_file) == 0:
        return dict()
    blast_df = BlastReader(blast_results_file).df
    blast_df = blast_df.loc[blast_df['pident'] >= min_ident]
    blast_df = blast_df.loc[blast_df['qcovhsp'] >= min_cov]
    blast_df = fixStart(blast_df)
    blast_df = blast_df.sort_values(['sseqid','sstart', 'send', 'bitscore'], ascending=[True, True, True, False])
    blast_df = blast_df.reset_index(drop=True)
    size = str(len(blast_df))
    prev_size = 0
    while size != prev_size:
        blast_df = filter_overlaping_records(blast_df, overlap, 'sseqid', 'sstart', 'send', 'bitscore')
        prev_size = size
        size = str(len(blast_df))

    return blast_df




def main():
    logging.info('Running plasmid detector v. {}'.format('0.1'))
    args = parse_args()
    if not args.infile:
        logging.info('Error, no blast file specified, please specify one')
        sys.exit()
    if not args.outdir:
        logging.info('Error, no output directory specified, please specify one')
        sys.exit()
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir, 0o755)

    blast_file = args.infile
    base_file_name = os.path.splitext(os.path.basename(blast_file))[0]
    out_dir = args.outdir
    blast_results_file = os.path.join(out_dir, base_file_name+'_blast_results.txt')
    processed_blast_results = filter_blast(blast_file, 95, 95, 0.00001, 5)
    if isinstance(processed_blast_results,dict):
        results_fh = open(blast_results_file, 'w')
        results_fh.write('')
        results_fh.close()
    else:
        processed_blast_results.to_csv(blast_results_file, sep='\t', header=True, line_terminator='\n', index=False)



main()