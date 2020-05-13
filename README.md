# MOB-suite: Software tools for clustering, reconstruction and typing of plasmids from draft assemblies

## Introduction ## 
Plasmids are mobile genetic elements (MGEs), which allow for rapid evolution and adaption of
bacteria to new niches through horizontal transmission of novel traits to different genetic
backgrounds. The MOB-suite is designed to be a modular set of tools for the typing and
reconstruction of plasmid sequences from WGS assemblies.


The MOB-suite depends on a series of databases which are too large to be hosted in git-hub. They can be downloaded or updated by running mob_init or if running any of the tools for the first time, the databases will download and initialize automatically if you do not specify an alternate database location. However, they are quite large so the first run will take a long time depending on your connection and speed of your computer.
Databases can be manually downloaded from https://share.corefacility.ca/index.php/s/rYaAH7oxrSVtilN/download or https://zenodo.org/record/3786915/files/data.tar.gz?download=1

### MOB-init
On first run of MOB-typer or MOB-recon, MOB-init should run to download the databases from figshare, sketch the databases and setup the blast databases. However, it can be run manually if the databases need to be re-initialized OR if you want to initialize the databases in an alternative directory.

```
% mob_init
```

### MOB-cluster
This tool creates plasmid similarity groups using fast genomic distance estimation using Mash.  Plasmids are grouped into clusters using complete-linkage clustering and the cluster code accessions provided by the tool provide an approximation of operational taxonomic units OTU’s. The plasmid nomenclature is designed to group highly similar plasmids together which are unlikely to have multiple representatives within a single cell and have a strong concordance with replicon and relaxase typing but is universally applicable since it uses the complete sequence of the plasmid itself rather than specific biomarkers.

### MOB-recon
This tool reconstructs individual plasmid sequences from draft genome assemblies using the clustered plasmid reference databases provided by MOB-cluster. It will also automatically provide the full typing information provided by MOB-typer. It optionally can use a chromosome depletion strategy based on closed genomes or user supplied filter of sequences to ignore.

### MOB-typer
Provides in silico predictions of the replicon family, relaxase type, mate-pair formation type and predicted transferability of the plasmid. Using a combination of biomarkers and MOB-cluster codes, it will also provide an observed host-range of your plasmid based on its replicon, relaxase and cluster assignment. This is combined with information mined from the literature to provide a prediction of the taxonomic rank at which the plasmid is likely to be stably maintained but it does not provide source attribution predictions.

## Installation ##

## Requires
+ Python v. 3.7 +
+ ete3 >= 3
+ biopython >= 1.70
+ pytables  >= 3.3
+ pycurl >= 7.43
+ pyqt  >= 5
+ numpy >= 1.18.1
+ scipy >= 1.4.1

## Dependencies

+ blast+ v. 2.3.0
+ mash


## Installation
We recommend MOB-Suite installation as a conda package due to large number of dependencies. The package is available through bioconda channel.

```
% conda config --add channels defaults
% conda config --add channels conda-forge
% conda config --add channels bioconda
% conda install -c bioconda mob_suite
```


### Pip

We recommend installing MOB-Suite via bioconda but you can install it via pip using the command below

```
% pip3 install mob_suite

```

### Docker image
A docker image is also available at [https://hub.docker.com/r/kbessonov/mob_suite](https://hub.docker.com/r/kbessonov/mob_suite)
```
% docker pull kbessonov/mob_suite:2.0.0 
% docker run --rm -v $(pwd):/mnt/ "kbessonov/mob_suite:2.0.0 " mob_recon -i /mnt/assembly.fasta -t -o /mnt/mob_recon_output

```

### Singularity image
A singularity image could be built via singularity recipe donated by Eric Deveaud. 
The recipe (`recipe.singularity`) is located in the singularity folder of this repository. 
The docker image section also has instructions on how to create singularity image from a docker image.

```bash
% singularity build mobsuite.simg recipe.singularity
```

## Using MOB-typer to perform replicon and relaxase typing of complete plasmids and to predict mobility and replicative plasmid host-range

### Setuptools
Clone this repository and install via setuptools. 

```
% git clone https://github.com/phac-nml/mob-suite.git
% cd mob-suite
% python setup.py install
```

## Using MOB-typer to perform replicon and relaxase typing of complete plasmids and predict mobility

You can perform plasmid typing using a fasta formated file containing a single plasmid represented by one or more contigs or it can treat all of the sequences in the fasta file as independant. The default behaviour is to treat all sequences in a file as from one plasmid, do not include multiple unrelated plasmids in the file without specifying --multi as they will be treated as a single plasmid.


```
# Single plasmid
% mob_typer --infile assembly.fasta --out_file sample_mobtyper_results.txt

# Multiple independant plasmids
% mob_typer --multi --infile assembly.fasta --out_file sample_mobtyper_results.txt

```

## Using MOB-recon to reconstruct plasmids from draft assemblies
This procedure works with draft or complete genomes and is agnostic of assembler choice but if
unicycler is used, then the circularity information can be parsed directly from the header of the unmodified assembly using -u . MOB-typing information is automatically generated for all plasmids reconstructed by MOB-recon.

```
### Basic Mode
% mob_recon --infile assembly.fasta --outdir my_out_dir
```

As of v. 3.0.0, we have added the ability of users to provide their own specific set of sequences to remove from plasmid reconstruction. This should be performed with caution and with the knowlede of your organism.  Sequences which are frequently of plasmid origin but are not in your organism is the primary use case we envision for this feature.

```
### User sequence mask
% mob_recon --infile assembly.fasta --outdir my_out_dir --
```

As of v. 3.0.0, we have provided the ability to use a collection of closed genomes which will be quickly checked using Mash for genomes which are genetically close and limit blast searches to those chromosomes. This more nuanced and automatic approach is recommended for users where there are sequences which should be filtered in one genomic context but not another. We provide as an optional download as set of closed Enterobacteriacea genomes from NCBI which can be used to provide added accuracy for some organisms such as E. coli and Klebsiella where there are sequences which switch between chromosome and plasmids.
<br><br>
If reconstructed plasmids exceed the Mash distance for primary cluster assignment, then they will get assigned a name in the format novel_{md5} where the md5 hash is calculated based on all of the sequences belonging to that reconstructed plasmid. This will provide a unique name for them but any change will result in a changed in the md5 hash. It is inadvised to use these groups for further analyses. Rather they should be highlighted as cases where targeted long read sequencing is required to obtain a closer database representitive of that plasmid.
```
### Autodetected close genome filter
% mob_recon --infile assembly.fasta --outdir my_out_dir -g 2019-11-NCBI-Enterobacteriacea-Chromosomes.fasta
```
## Using MOB-cluster
Use this tool only to update the plasmid databases or build a new one and should only be completed with closed high quality plasmids. If you add in poor quality data it can severely impact MOB-recon. As od v. 3.0.0, MOB-cluster has been re-written to utilize the output from MOB-typer to greatly speed up the process of updating and builing plasmid databases by using pre-computed results. Clusters generated from earlier versions of MOB-suite are not compatibile with the new clusters. We have povided a mapping file of previous cluster assignments and their new cluster accessions. Each cluster code is unique and will not be re-used.

```
### Build a new database
% mob_cluster --mode build -f new_plasmids.fasta -p new_plasmids_mobtyper_report.txt -t new_plasmids_host_taxonomy.txt --outdir output_directory
```

```
### Add a sequence to an existing database
% mob_cluster --mode update -f new_plasmids.fasta -p new_plasmids_mobtyper_report.txt -t new_plasmids_host_taxonomy.txt --outdir output_directory -c existing_clusters.txt -r existing_sequences.fasta
```

```
### Update MOB-suite plasmid databases
% cp output_directory/clusters.txt
% mv output_directory/updated.fasta mob_db_path/ncbi_plasmid_full_seqs.fas
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
| mobtyper_results | Aggregate MOB-typer report files for all identified plasmid |

# MOB-recon contig report format
| field  | Description |
| --------- |  --------- | 
| sample_id | Sample ID specified by user or deault to filename |
| molecule_type | Plasmid or Chromosome |
| primary_cluster_id | primary MOB-cluster id of neighbor |
| secondary_cluster_id | secondary MOB-cluster id of neighbor |
| size | Length in base pairs |
| gc | GC % |
| md5 | md5 hash |
| circularity_status | Molecule is either circular, incomplete or not tested based on parameters used |
| rep_type(s) | Replion type(s) |
| rep_type_accession(s) | Replicon sequence accession(s) |
| relaxase_type(s) | Relaxase type(s) |
| relaxase_type_accession(s) | Relaxase sequence accession(s) |
| mpf_type | Mate-Pair formation type |
| mpf_type_accession(s) | Mate-Pair formation sequence accession(s) |
| orit_type(s) | Origin of transfer type |
| orit_accession(s) | Origin of transfer sequence accession(s) |
| predicted_mobility | Mobility prediction for the plasmid (Conjugative, Mobilizable, Non-mobilizable) |
| mash_nearest_neighbor | Accession of closest plasmid database match |
| mash_neighbor_distance | Mash distance from query to match |
| mash_neighbor_identification | Host taxonomy of the plasmid database match |
| repetitive_dna_id | Repetitive DNA match id |
| repetitive_dna_type| Repetitive element class |



# MOB-typer report file format
| field  | Description |
| --------- |  --------- | 
| sample_id | Sample ID specified by user or deault to filename |
| num_contigs | Number of sequences belonging to plasmid |
| size | Length in base pairs |
| gc | GC % |
| md5 | md5 hash |
| rep_type(s) | Replion type(s) |
| rep_type_accession(s) | Replicon sequence accession(s) |
| relaxase_type(s) | Relaxase type(s) |
| relaxase_type_accession(s) | Relaxase sequence accession(s) |
| mpf_type | Mate-Pair formation type |
| mpf_type_accession(s) | Mate-Pair formation sequence accession(s) |
| orit_type(s) | Origin of transfer type |
| orit_accession(s) | Origin of transfer sequence accession(s) |
| predicted_mobility | Mobility prediction for the plasmid (Conjugative, Mobilizable, Non-mobilizable) |
| mash_nearest_neighbor | Accession of closest plasmid database match |
| mash_neighbor_distance | Mash distance from query to match |
| mash_neighbor_identification | Host taxonomy of the plasmid database match |
| primary_cluster_id | primary MOB-cluster id of neighbor |
| secondary_cluster_id | secondary MOB-cluster id of neighbor |
| predicted_host_range_overall_rank | Taxon rank of convergence between observed and reported host ranges |
| predicted_host_range_overall_name | Taxon name of convergence between observed and reported host ranges |
| observed_host_range_ncbi_rank | Taxon rank of convergence of plasmids in MOB-suite plasmid DB |
| observed_host_range_ncbi_name | Taxon name of convergence of plasmids in MOB-suite plasmid DB |
| reported_host_range_lit_rank | Taxon rank of convergence of literature reported host ranges |
| reported_host_range_lit_name | Taxon name of convergence of literature reported host ranges |
| associated_pmid(s) | PubMed ID(s) associated with records |

# MOB-cluster sequence cluster information file
| field  | Description |
| --------- |  --------- | 
| sample_id | Sample ID specified by user or deault to filename |
| size | Length in base pairs |
| gc | GC % |
| md5 | md5 hash |
| organism | Host taxon name |
| taxid | Host NCBI taxon id |
| rep_type(s) | Replion type(s) |
| rep_type_accession(s) | Replicon sequence accession(s) |
| relaxase_type(s) | Relaxase type(s) |
| relaxase_type_accession(s) | Relaxase sequence accession(s) |
| mpf_type | Mate-Pair formation type |
| mpf_type_accession(s) | Mate-Pair formation sequence accession(s) |
| orit_type(s) | Origin of transfer type |
| orit_accession(s) | Origin of transfer sequence accession(s) |
| predicted_mobility | Mobility prediction for the plasmid (Conjugative, Mobilizable, Non-mobilizable) |
| primary_cluster_id | primary MOB-cluster id of plasmid |
| primary_dist | primary MOB-cluster distance cutoff to generate cluster |
| secondary_cluster_id | secondary MOB-cluster id of plasmid |
| secondary_dist | secondary MOB-cluster distance cutoff to generate cluster |



# blast report file format
| field name | description|
| -----------| -----------|
| qseqid | query sequence id |
| sseqid | subject sequence id |
| qlen | query length |
| slen | subject length |
| qstart | match start query |
| qend | match end query |
| sstart | match subject start|
| send | match subject end|
| length | length of alignment|
| mismatch | number of mismatches|
| pident | identity|
| qcovhsp | query coverage by hsp|
| qcovs | query coverage by subject|
| sstrand | strad of hit in subject|
| evalue | evalue of match|
| bitscore | bitscore of match |



## Contact

James Robertson - james.robertson@canada.ca
Kyrylo Bessonov - kyrylo.bessonov@canada.ca

## License

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in
compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is
distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
implied. See the License for the specific language governing permissions and limitations under the
License.
