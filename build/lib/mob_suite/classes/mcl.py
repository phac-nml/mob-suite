#!/usr/bin/env python

from subprocess import Popen, PIPE
import os

class mcl:


    def __init__(self, blast_results_file, working_dir,inflation=1.5):
        blast_abc = os.path.join(working_dir, 'seq.abc')
        cols = [0,1,14]
        self.prep_blast(blast_results_file, blast_abc, cols)
        mci_outfile = os.path.join(working_dir, 'seq.mci')
        tab_outfile = os.path.join(working_dir, 'seq.tab')
        clust_outfile = os.path.join(working_dir, 'seq.clust')
        self.mcxload(blast_abc,mci_outfile, tab_outfile)
        self.run_mcl(mci_outfile,tab_outfile,clust_outfile,inflation,1)
        self.clusters = self.parse_mcl(clust_outfile)

    def getclusters(self):
        return self.clusters


    def prep_blast(self,blast_results_file,outfile,col_numbers):
        with open(blast_results_file) as f:
            content = f.readlines()
            content = [x.strip() for x in content]
            f.close()

        with open(outfile,'w') as f:
            for line in content:
                row = line.split("\t")
                selection = list()
                for col in col_numbers:
                    selection.append(row[col])
                f.write("\t".join(selection) + "\n")
            f.close()



    def  mcxload(self, blast_results_file,mci_outfile, tab_outfile):
        p = Popen(['mcxload',
                   '-abc', blast_results_file,
                   '--stream-mirror',
                   '--stream-neg-log10',
                   '-o', mci_outfile,
                   '-write-tab', tab_outfile,],
                  stdout=PIPE,
                  stderr=PIPE)


        p.wait()
        stdout = p.stdout.read()
        print(stdout)
        stderr = p.stderr.read()
        print(stderr)

    def run_mcl(self,mci_outfile,tab_outfile,cluster_outfile,inflation,num_threads):
        p = Popen(['mcl',
                   mci_outfile,
                   '-I',str(inflation),
                   '-use-tab', tab_outfile,
                   '-o',cluster_outfile],
                  stdout=PIPE,
                  stderr=PIPE)

        p.wait()
        stdout = p.stdout.read()
        stderr = p.stderr.read()

    def parse_mcl(self,cluster_outfile):
        with open(cluster_outfile) as f:
            content = f.readlines()
            content = [x.strip() for x in content]
        f.close()
        count = 0
        clusters = dict()
        for line in content:
            members = line.split("\t")
            for member in members:
                clusters[member] = count
            count+=1
        return clusters

