from subprocess import Popen, PIPE, check_call
import os, sys, logging
from mob_suite.blast import BlastRunner
from mob_suite.blast import BlastReader

class mash:
    def __init__(self):
        return

    def run_mash(self, reference_db, input_fasta, table=False, num_threads=1):

        if table:
            p = Popen(['mash', "dist", "-t", "-p", str(num_threads), reference_db, input_fasta],
                      stdout=PIPE,
                      stderr=PIPE)
            (stdout, stderr) = p.communicate()
        else:
            p = Popen(['mash', "dist", "-p", str(num_threads), reference_db, input_fasta],
                      stdout=PIPE,
                      stderr=PIPE)
            (stdout, stderr) = p.communicate()
        if len(str(stderr)) > 0 and str(stderr) != "b''":
            logging.info('{}'.format(stderr))
        return stdout

    def run_mash_screen(self, reference_db, input_fasta, winner_take_all=True, num_threads=1):

        if winner_take_all:
            p = Popen(['mash', "screen", "-p", str(num_threads), '-w', '-i','0', reference_db, input_fasta],
                      stdout=PIPE,
                      stderr=PIPE)
            (stdout, stderr) = p.communicate()
        else:
            p = Popen(['mash', "screen", "-p", str(num_threads), '-i','0', reference_db, input_fasta],
                      stdout=PIPE,
                      stderr=PIPE)
            (stdout, stderr) = p.communicate()

        logging.info(
            '{}'.format(stderr))
        return stdout


    def mashsketch(self, input_fasta, output_path, sketch_ind=True, num_threads=1, kmer_size=21, sketch_size=1000):
        if output_path == '':
            os.path.dirname(input_fasta)

        if sketch_ind:
            p = Popen(['mash', "sketch",
                       "-p", str(num_threads),
                       "-i",
                       "-o", output_path,
                       "-k", str(kmer_size),
                       "-s", str(sketch_size), input_fasta],
                      stdout=PIPE,
                      stderr=PIPE)
        else:
            p = Popen(['mash', "sketch",
                       "-p", str(num_threads),
                       "-o", output_path,
                       "-k", str(kmer_size),
                       "-s", str(sketch_size), input_fasta],
                      stdout=PIPE,
                      stderr=PIPE)
        p.wait()
        stdout = p.stdout.read()
        stderr = p.stderr.read()


class detectCircularity:
    ### Method adapted from Berokka https://github.com/tseemann/berokka by Torsten Seemann
    def __init__(self):
        return

    def run(self,input_fasta, output_path, logging, min_cov=1, min_ident=1, evalue=1, num_threads=1, min_length=25):
        blast_results_file = os.path.join(output_path, 'circularize.blast.txt')
        self.run_blast(input_fasta=input_fasta,
                       output_path=output_path,
                       blast_results_file=blast_results_file,
                       logging=logging,
                       min_cov=min_cov,
                       min_ident=min_ident,
                       evalue=evalue,
                       num_threads=num_threads,
                       min_length=min_length)
        return (self.overhangDetection(blast_results_file,logging))

    def run_blast(self,input_fasta,output_path,blast_results_file,logging,min_cov=1,min_ident=1,evalue=1,num_threads=1,min_length=25):
        blast_runner = BlastRunner(input_fasta, output_path)
        blast_runner.makeblastdb(input_fasta, 'nucl',logging)
        blast_runner.run_blast(query_fasta_path=input_fasta, blast_task='megablast', db_path=input_fasta,
                               db_type='nucl', min_cov=min_cov, min_ident=min_ident, evalue=evalue,
                               blast_outfile=blast_results_file, num_threads=num_threads, word_size=11,logging=logging)

        if os.path.getsize(blast_results_file) == 0:
            fh = open(blast_results_file, 'w', encoding="utf-8")
            fh.write('')
            fh.close()
            return dict()

        blast_df = BlastReader(blast_results_file,logging).df
        blast_df = blast_df.loc[blast_df['length'] >= min_length]
        blast_df = blast_df.reset_index(drop=True)
        blast_df.to_csv(blast_results_file, sep='\t', header=False, line_terminator='\n', index=False)

    def overhangDetection(self,blast_results_file,logging,min_length=25):
        if os.path.getsize(blast_results_file) == 0:
            return dict()

        blast_df = BlastReader(blast_results_file,logging).df.sort_values(['qseqid', 'qstart', 'qend', 'bitscore'], ascending=[True, True, True, False])

        circular_contigs = {}


        for index, row in blast_df.iterrows():
            contig_id_query = row['qseqid']
            contig_id_subject = row['sseqid']
            contig_start_subject = int(row['sstart'])
            contig_end_subject = int(row['send'])
            contig_start_query = int(row['qstart'])
            contig_end_query = int(row['qend'])
            contig_length = int(row['qlen'])
            length = int(row['length'])


            if contig_id_query != contig_id_subject and contig_id_subject != "ref|{}|".format(contig_id_query):
                continue

            if contig_start_query != 1 or length < min_length:

                continue

            if contig_start_query  == contig_start_subject and contig_end_query == contig_end_subject:

                continue



            if contig_start_query == 1  and contig_end_subject == contig_length:
                circular_contigs[contig_id_query] = 'Circular: Overlap {} bp'.format(length)
        return circular_contigs





