from subprocess import Popen, PIPE, check_call
import os, sys, logging


class circlator:
    def __init__(self):
        return

    def run_minimus(self, input_fasta, output_prefix):
        p = Popen(['circlator', "minimus2", input_fasta, output_prefix],
                  stdout=PIPE,
                  stderr=PIPE)
        (stdout, stderr) = p.communicate()
        logging.info(
            '{}'.format(stderr))

    def parse_minimus(self, minimuslog_file):
        if not os.path.isfile(minimuslog_file):
            logging.info('Error minimus2 file  {} not found, confirm that circulator minimus2 is functional'.format(
                minimuslog_file))
            return list()
        circularized_ids = list()
        fh = open(minimuslog_file, 'r')
        content = fh.readlines()
        content = [x.strip() for x in content]
        for line in content:
            if line[0:21] == 'Circularised contigs:':
                circularized_ids = line.split("\t")

        return circularized_ids[1:]


class mash:
    def __init__(self):
        return

    def run_mash(self, reference_db, input_fasta, output_filehandle, table=False, num_threads=1):

        if table:
            p = Popen(['mash', "dist", "-t", "-p", str(num_threads), reference_db, input_fasta],
                      stdout=output_filehandle,
                      stderr=PIPE)
            (stdout, stderr) = p.communicate()
        else:
            p = Popen(['mash', "dist", "-p", str(num_threads), reference_db, input_fasta],
                      stdout=output_filehandle,
                      stderr=PIPE)
            (stdout, stderr) = p.communicate()

        logging.info(
            '{}'.format(stderr))
        output_filehandle.close()

    def read_mash(selfs, mashfile):
        fh = open(mashfile, 'r')
        return fh.readlines()

    def mashsketch(self, input_fasta, output_path, sketch_ind=True, num_threads=1, kmer_size=21, sketch_size=1000):
        if output_path == '':
            os.path.dirname(input_fasta)
        p = Popen(['mash', "sketch",
                   "-p", str(num_threads),
                   "-i",
                   "-o", output_path,
                   "-k", str(kmer_size),
                   "-s", str(sketch_size), input_fasta],
                  stdout=PIPE,
                  stderr=PIPE)
        p.wait()
        stdout = p.stdout.read()
        stderr = p.stderr.read()





