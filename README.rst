## WARNING: This is still in active testing, please report any bugs to james.robertson@canada.ca. This should only be used for testing purposes.

## Introduction

Plasmids are mobile genetic elements (MGEs), which allow for rapid evolution and adaption of bacteria to new niches through horizontal transmission of novel traits to different genetic backgrounds. Plasmids and other MGE’s are of great concern to public health since they can encode for traits such as antimicrobial resistance (AMR), toxin production, and virulence. The mobile nature of plasmids can allow for rapid dissemination of AMR through a population, and many plasmids are readily transmissible. Understanding the transmission dynamics of AMR plasmids is critical for public health professionals in order to accurately respond to the rising threat of multi-drug resistant bacteria. Mobility is key to the long-term survival of plasmids and understanding the epidemiology of plasmid-borne traits requires insight into how these vectors are transmitted between hosts. Classification schemes for plasmids based on their mobility have been the subject of several papers which have broadly grouped plasmids based on their relaxase gene and type IV coupling protein (T4CP) into the following categories i) conjugative ii) mobilizable iii) non-mobilizable. Six major relaxase families (MOB) have been previously been identified by relaxase typing based on specific sequence motifs. Currently, there are no tools available to type plasmid assemblies automatically according to MOB groups and predict their mobility. This includes tools such as Plasmid Finder https://cge.cbs.dtu.dk/services/PlasmidFinder/, which focus on the typing of plasmids within draft bacterial genomes, based on the detection of the replication gene(s). Plasmid family information derived from replicon typing can provide information on host range, frequency of AMR elements, among others but it cannot identify whether a given member of a plasmid family possesses the necessary components to transfer itself into another host through conjugation. 

## Release 0.1

MOB-Suite is a collection of three tools, which are designed to assist with mobilome research.

1)	MOB-cluster: This tool creates plasmid similarity groups using fast genomic distance estimation using MASH.  Plasmids are grouped into clusters using single-linkage clustering and the cluster codes provided by the tool provide an approximation of operational taxonomic units OTU’s 
2)	MOB-recon: This tool reconstructs individual plasmid sequences from draft genome assemblies using the clustered plasmid reference databases provided by MOB-cluster.
3)	MOB-typer: provides in silico predictions of the replicon family, relaxase type, mate-pair formation type and predicted transferability of the plasmid

## Dependancies
Blast+ https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
Mash
http://mash.readthedocs.io/en/latest/index.html

Circlator:
http://sanger-pathogens.github.io/circlator/

## Installation

1) Clone the mob-suite repository
2) Install dependencies programs, blast+ and circlator
3) download the needed databases from: https://figshare.com/s/a4c92dd84f17b2cefea6
4) Uncompress the files and copy them into the mob-suite repository mob-suite/mob_suite/databases
5) Navigate to the mob-suite repository
6) Install the program using python setup.py install

##MOB-cluster Use Instructions
Use this tool only to update the plasmid databases or build a new one and should only be completed with closed high quality plasmids

Build a new database:
mob_cluster --mode build --infile plasmid.fasta --outdir output_directory

Add a sequence to an existing database
mob_cluster --mode update --infile plasmid.fasta --outdir output_directory

If you want your database to replace the existing mob-suite database you will need to either replace the existing db and rebuild the blast and mash databases
makeblastdb -in ncbi_plasmid_full_seqs.fas -dbtype nucl
mash sketch -i ncbi_plasmid_full_seqs.fas ncbi_plasmid_full_seqs.fas

If you just want to test a new plasmid database you can specify the parameters below and give the locations of your new databases
--plasmid_mash_db
--plasmid_db

##MOB-typer Use Instructions
Use this tool on fasta files which contain only one plasmid. The plasmid may be complete or represented in fragments

Basic Use:
mob_typer --infile plasmid.fasta --outdir output_directory

##MOB-recon Use Instructions
Use this tool to reconstruct plasmids from mixed genome assemblies of plasmids and chromosomes

Basic Use:
mob_recon --infile assembly.fasta --outdir output_directory




## Contact

James Robertson - james.robertson@canada.ca

## License

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
