# MOB-suite: Software tools for clustering, reconstruction and typing of plasmids from draft assemblies

## Introduction ## 
Plasmids are mobile genetic elements (MGEs), which allow for rapid evolution and adaption of
bacteria to new niches through horizontal transmission of novel traits to different genetic
backgrounds. The MOB-suite is designed to be a modular set of tools for the typing and
reconstruction of plasmid sequences from WGS assemblies.

### MOB-cluster
This tool creates plasmid similarity groups using fast genomic distance estimation using MASH.  Plasmids are grouped into clusters using single-linkage clustering and the cluster codes provided by the tool provide an approximation of operational taxonomic units OTU’s 

### MOB-recon
This tool reconstructs individual plasmid sequences from draft genome assemblies using the clustered plasmid reference databases provided by MOB-cluster.

### MOB-typer
Provides in silico predictions of the replicon family, relaxase type, mate-pair formation type and predicted transferability of the plasmid

## Installation ## 

## Dependancies
```
% conda config --add channels defaults
% conda config --add channels conda-forge
% conda config --add channels bioconda
% conda install blast mash circlator
```


### Pip
Currently, the only installation method available is by installing through pip but there are plans to package it for bioconda
```
% git clone https://github.com/jrober84/mob-suite
% cd mob-suite
% python install setup.py

```
## Using MOB-typer to perform replicon and relaxase typing of complete plasmids and predict mobility

You can perform plasmid typing using a fasta formated file containing a single plasmid represented by one or more contigs. Do not include multiple unrelated plasmids in the file as they will be treated as a single plasmid.

```
# Basic Mode
% mob_typer.py --infile assembly.fasta --outdir my_out_dir

# Look for a file called mobtyper_(input_file)_report.txt
% cat my_out_dir/mobtyper_(input_file)_report.txt
```

## Using MOB-recon to reconstruct plasmids from draft assemblies
This procedure works with draft or complete genomes and is agnostic of assembler choice but if
unicycler is used, then the circularity information can be parsed directly from the header of the unmodified assembly.

```
### Basic Mode
% mob_recon.py --infile assembly.fasta --outdir my_out_dir
```

```
### Full Mode
# In this mode, MOB-typer will be run on each identified plasmid grouping and will produce a summary report
% mob_recon.py --infile assembly.fasta --outdir my_out_dir --run_typer
```

## Using MOB-cluster
Use this tool only to update the plasmid databases or build a new one and should only be completed with closed high quality plasmids. If you add in poor quality data it will severely impact MOB-recon

```
### Build a new database
% mob_cluster.py --mode build --infile plasmid.fasta --outdir output_directory
```

```
### Add a sequence to an existing database
% mob_cluster --infile update_sequences.fasta --ref_fasta_file original.fasta --ref_mash_db original.msh --ref_cluster_file original_clusters.txt 
```

```
### Test new plasmid database with MOB-recon
% makeblastdb -in path_to_plasmid_testing_db -dbtype nucl
% mash sketch -i path_to_plasmid_testing_db   <---- produces mash sketch file with format "path_to_plasmid_testing_db.msh"
% mob_recon.py --infile assembly.fasta --outdir my_out_dir --run_typer --plasmid_mash_db path_to_mash_testing_db --plasmid_db path_to_plasmid_testing_db
```


```
### Update MOB-suite plasmid databases
% mv new_mob_formated_db.fasta mob_db_path/ncbi_plasmid_full_seqs.fas
% makeblastdb -in mob_db_path/ncbi_plasmid_full_seqs.fas -dbtype nucl
% mash sketch -i mob_db_path/ncbi_plasmid_full_seqs.fas 
```




# Output files
| file | Description |
| ------------ | ------------ |
| contig_report.txt | This file describes the assignment of the contig to chromosome or a particular plasmid grouping |
| repetitive_blast_report | Summary information of contigs found to consist of nothing but a repetitive element |
| chromosome.fasta | Fasta file of all contigs found to belong to the chromosome |
| plasmid_(X).fasta | Each plasmid group is written to an individual fasta file which contains the assigned contigs |
| mobtyper_(input_file)_report.txt | Individual MOB-typer report files for each identified plasmid |
| mobtyper_aggregate_report.txt | Aggregate MOB-typer report files for all identified plasmid |

# MOB-recon contig report format
| field id | description |
| -------- | ------------|
| file_id | Name of the input file  |
| cluster_id | MOB-cluster type of reference match |
| contig_id | Unique identifier of the contig |
| contig_length | Length of the contig |
| circularity_status | Circular if Circlator or Unicycler find it to be circular, and incomplete if not |
| rep_type | Replicon types idenfied |
| rep_type_accession | Accessions of replicons identified |
| relaxase_type | Relaxase types identified |
| relaxase_type_accession | Accessions of relaxases identified |
| mash_nearest_neighbor | Mate-pair formation types identified |
| mash_neighbor_distance | Mate-pair formation type accessioons |
| repetitive_dna_id | Repetitive DNA match id |
| match_type | Repetitive element class |
| score | Blast bitscore of match |
| contig_match_start | Start of match on contig |
| contig_match_end | End of match on contig |



# MOB-typer report file format
| field name | description|
| -----------| -----------|
| file_id | Name of the input file |
| num_contigs | Number of sequences identified in the file |
| total_length | Total number of bases in all sequences |
| gc | GC% of all sequences |
| rep_type(s) | Replicon types idenfied |
| rep_type_accession(s) | Accessions of replicons identified |
| relaxase_type(s) | Relaxase types identified |
| relaxase_type_accession(s) | Accessions of relaxases identified |
| mpf_type | Mate-pair formation types identified |
| mpf_type_accession(s) | Mate-pair formation type accessioons |
| orit_type(s) | Relaxase type of oriT sequence |
| orit_accession(s) | Accession for oriT |
| PredictedMobility | Mobility prediction for the plasmid (Conjugative, Mobilizable, Non-mobilizable) |
| mash_nearest_neighbor | Accession of closest database match |
| mash_neighbor_distance | Mash distance from query to match |
| mash_neighbor_cluster | MOB-cluster type of reference match |



## Contact

James Robertson - james.robertson@canada.ca

## License

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in
compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is
distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
implied. See the License for the specific language governing permissions and limitations under the
License.
