import pandas as pd
import numpy as np
from pandas.io.common import EmptyDataError
from blast import BlastReader
import sys
from collections import OrderedDict
from operator import itemgetter




class mge_predict:
    circular_contigs = list()
    valid_clusters = dict()
    cluster_scores = dict()
    blast_df =  pd.DataFrame()
    seq_sizes = dict()
    seq_scores = dict()
    contig_coverage = dict()
    database_seq_coverage = dict()
    mask_ref_contig_ids = dict()
    mask_ref_cluster_ids = dict()
    plasmid_contigs = dict()

    def __init__(self,blast_results):
        blast_reader = BlastReader(blast_results)
        self.blast_df = blast_reader.df
        self.filter_dataframe(1000,400000,60)

        self.fixStart()
        self.blast_df = self.blast_df.reset_index(drop=True)
        self.blast_df = self.blast_df.sort_values(['qseqid', 'qstart','qend'], ascending=[True, True, True])
        #print self.blast_df

        self.blast_df = self.blast_df.reset_index(drop=True)
        self.init_contig_info()
        contig_covered_bases = self.calc_covered_seq_bases('qseqid', 'qstart', 'qend')
        #print contig_covered_bases
        #print contig_covered_bases
        self.contig_coverage = self.calc_perc_coverage(contig_covered_bases)
        self.blast_df = self.blast_df.sort_values(['sseqid', 'sstart', 'send'], ascending=[True, True, True])
        self.blast_df = self.blast_df.reset_index(drop=True)
        database_covered_bases = self.calc_covered_seq_bases('sseqid', 'sstart', 'send')
        #print database_covered_bases
        self.database_seq_coverage = self.calc_perc_coverage(database_covered_bases)

        self.filter_low_cov_hits(60)
        self.blast_df = self.blast_df.reset_index(drop=True)

        for id in self.mask_ref_contig_ids:
            self.blast_df = self.blast_df[self.blast_df["sseqid"]  != id]
            if id in self.seq_sizes:
                del(self.seq_sizes[id])
                del(self.database_seq_coverage[id])
        self.blast_df = self.blast_df.reset_index(drop=True)

        for contig_id in self.seq_scores:
            scores = dict()
            for match_id in self.seq_scores[contig_id]:
                if match_id in self.mask_ref_contig_ids:
                    continue
                scores[match_id] = self.seq_scores[contig_id][match_id]

            self.seq_scores[contig_id] = OrderedDict(sorted(list(scores.items()), key=itemgetter(1), reverse=True))
            top_match = -1
            scores = dict()
            for match_id in self.seq_scores[contig_id]:
                if self.seq_scores[contig_id][match_id] > top_match:
                    top_match = self.seq_scores[contig_id][match_id]
                if self.seq_scores[contig_id][match_id] / top_match * 100 < 90:
                    continue
                scores[match_id] = self.seq_scores[contig_id][match_id]
            self.seq_scores[contig_id] = OrderedDict(sorted(list(scores.items()), key=itemgetter(1), reverse=True))


        temp = dict()
        for contig_id in self.seq_scores:
            if len(self.seq_scores[contig_id]) > 0:
                temp[contig_id] = self.seq_scores[contig_id]
                for match_id in self.seq_scores[contig_id]:
                    score = self.seq_scores[contig_id][match_id]
                    acs,cluster_id = match_id.split('|')
                    if not cluster_id in self.cluster_scores:
                        self.cluster_scores[cluster_id] = 0
                    self.cluster_scores[cluster_id]+= score
        self.seq_scores = temp

        for contig_id in self.seq_scores:
            if len(self.seq_scores[contig_id]) == 1:
                continue
            prev_cluster_id = -1
            temp = dict()
            for match_id in self.seq_scores[contig_id]:
                acs, cluster_id = match_id.split('|')
                score = self.cluster_scores[cluster_id]
                if prev_cluster_id == -1:
                    prev_cluster_id = cluster_id
                    prev_cluster_score = score
                    temp[match_id] = self.seq_scores[contig_id][match_id]

                if float(score)/prev_cluster_score*100 < 90:
                    continue

                temp[match_id] = self.seq_scores[contig_id][match_id]
            self.seq_scores[contig_id] = temp

        for contig_id in self.seq_scores:
            placement = 'single'
            matches = self.seq_scores[contig_id]


            clusters = dict()
            match_ids = dict()
            sizes = list()
            for match_id in matches:
                acs, cluster_id = match_id.split('|')
                clusters[cluster_id] = acs
                sizes.append(self.seq_sizes[match_id])

            mininimum = min(sizes)
            maximum = max(sizes)
            if mininimum != maximum:
                size = str(mininimum) +'-' + str(maximum)
            else:
                size = mininimum
            cluster_ids = ','.join(list(clusters.keys()))
            if len(clusters) > 1:
                placement = 'multiple'
            self.plasmid_contigs[contig_id] = {
                'contig_length': self.seq_sizes[contig_id],
                'plasmid_membership': placement,
                'plasmid_cluster_ids': cluster_ids,
                'contig_cov': self.contig_coverage[contig_id],
                'cluster_size_range':size
            }





    def filter_dataframe(self, query_len_min, query_len_max,query_cov ):

        self.blast_df = self.blast_df.loc[self.blast_df['qlen'] <= query_len_max]
        self.blast_df = self.blast_df.loc[self.blast_df['qlen'] >= query_len_min]
        self.blast_df = self.blast_df.loc[self.blast_df['qcovs'] >= query_cov]
        self.blast_df = self.blast_df.reset_index(drop=True)


    def fixStart(self):
        for index, row in self.blast_df.iterrows():
            sstart = self.blast_df.at[index, 'sstart']
            send = self.blast_df.at[index, 'send']
            #print "{}\t{}".format(sstart,send)
            if send < sstart:
                temp = sstart
                self.blast_df.at[index, 'sstart'] = send
                self.blast_df.at[index, 'send'] = temp
            #print "====>{}\t{}".format(self.blast_df.at[index, 'sstart'], self.blast_df.at[index, 'send'])
            qstart = self.blast_df.at[index, 'qstart']
            qend = self.blast_df.at[index, 'qend']
            if qend < qstart:
                temp = qstart
                self.blast_df.at[index, 'qstart'] = qend
                self.blast_df.at[index, 'qend'] = temp

    def get_seq_cov_ranges(self,id_col_name,start_col_name,end_col_name):
        ranges = dict()

        for index, row in self.blast_df.iterrows():
            seq_id = row[id_col_name]
            start = row[start_col_name]
            end = row[end_col_name]
            #print "{}\t{}\t{}".format(seq_id,start,end)

            if not seq_id in ranges:
                ranges[seq_id] = list()
            ranges[seq_id].append((start,end))
        #print ranges

        return ranges

    def summarize_ranges(self,ranges):
        prev_start = -1
        prev_end = -1
        summary = list()

        for start,end in ranges:
            #print "{}\t{}\t{}\t{}".format(prev_start,prev_end,start,end)
            if prev_start == -1 :
                prev_start = start
                prev_end = end
                continue

            if start >= prev_start and start <= prev_end:
                if end > prev_end:
                    prev_end = end
                continue

            if start > prev_end:
                summary.append((prev_start,prev_end))
                prev_start = -1
                prev_end= -1
                
        if prev_start == -1:
            summary.append((start, end))
        else:
        	summary.append((prev_start, prev_end))

            

        #print summary
        #sys.exit()
        return summary


    def calc_covered_seq_bases(self,id_col_name,start_col_name,end_col_name):
        covered_bases = dict()
        ranges = self.get_seq_cov_ranges(id_col_name,start_col_name,end_col_name)
        #print ranges

        for seq_id in ranges:
            summary = self.summarize_ranges(ranges[seq_id])
            #print summary

            sum = 0
            for start,end in summary:
                sum += end - start
            covered_bases[seq_id] = sum+1
        #sys.exit()
        return covered_bases

    def calc_perc_coverage(self, covered_bases):
        for seq_id in covered_bases:
            if seq_id in self.seq_sizes:

                covered_bases[seq_id] = float(covered_bases[seq_id]) / self.seq_sizes[seq_id]*100
            else:
                covered_bases[seq_id] = -1


        return covered_bases

    def filter_low_cov_hits(self,cov):
        for id in self.database_seq_coverage:
            if self.database_seq_coverage[id] < cov:
                self.mask_ref_contig_ids[id] = 'Coverage Too Low'




    def init_contig_info(self):
        seq_scores = dict()
        seq_sizes = dict()
        for index, row in self.blast_df.iterrows():
            if not row['qseqid'] in seq_scores:
                seq_scores[row['qseqid']] = dict()


            if not row['sseqid'] in seq_scores[row['qseqid']] :
                seq_scores[row['qseqid']][row['sseqid']] = 0

            if seq_scores[row['qseqid']][row['sseqid']] < row['bitscore']:
                seq_scores[row['qseqid']][row['sseqid']] = row['bitscore']

            seq_sizes[row['sseqid']] = row['slen']
            seq_sizes[row['qseqid']] = row['qlen']
        self.seq_sizes = seq_sizes
        self.seq_scores = seq_scores

    def get_plasmid_contigs(self):
        return self.plasmid_contigs



