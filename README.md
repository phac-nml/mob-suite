![](https://img.shields.io/conda/dn/bioconda/mob_suite)
![](https://img.shields.io/docker/pulls/kbessonov/mob_suite)
![](https://img.shields.io/pypi/dm/mob-suite)
![](https://img.shields.io/github/v/release/phac-nml/mob-suite?include_prereleases)
![](https://img.shields.io/github/last-commit/phac-nml/mob-suite)
![](https://img.shields.io/github/issues/phac-nml/mob-suite)

# MOB-suite: Software tools for clustering, reconstruction and typing of plasmids from draft assemblies

## Introduction
Plasmids are mobile genetic elements (MGEs), which allow for rapid evolution and adaption of
bacteria to new niches through horizontal transmission of novel traits to different genetic
backgrounds. The MOB-suite is designed to be a modular set of tools for the typing and
reconstruction of plasmid sequences from WGS assemblies.


The MOB-suite depends on a series of databases which are too large to be hosted in git-hub. They can be downloaded or updated by running mob_init or if running any of the tools for the first time, the databases will download and initialize automatically if you do not specify an alternate database location. However, they are quite large so the first run will take a long time depending on your connection and speed of your computer.
Databases can be manually downloaded from [here](https://zenodo.org/records/10304948/files/data.tar.gz?download=1). <br>
Our new automatic chromosome depletion feature in MOB-recon can be based on any collection of closed chromosome sequences.

## Citations
Below are the manuscripts describing the algorithmic approaches used in the MOB-suite. 

1) Robertson, James, and John H E Nash. “MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies.” Microbial genomics vol. 4,8 (2018): e000206. doi:10.1099/mgen.0.000206

2) Robertson, James et al. “Universal whole-sequence-based plasmid typing and its utility to prediction of host range and epidemiological surveillance.” Microbial genomics vol. 6,10 (2020): mgen000435. doi:10.1099/mgen.0.000435

### MOB-init
On first run of MOB-typer or MOB-recon, MOB-init (invoked by `mob_init` command) should run to download the databases from figshare, sketch the databases and setup the blast databases. However, it can be run manually if the databases need to be re-initialized OR if you want to initialize the databases in an alternative directory.


### MOB-cluster
This tool creates plasmid similarity groups using fast genomic distance estimation using Mash.  Plasmids are grouped into clusters using complete-linkage clustering and the cluster code accessions provided by the tool provide an approximation of operational taxonomic units OTU’s. The plasmid nomenclature is designed to group highly similar plasmids together which are unlikely to have multiple representatives within a single cell and have a strong concordance with replicon and relaxase typing but is universally applicable since it uses the complete sequence of the plasmid itself rather than specific biomarkers.

### MOB-recon
This tool reconstructs individual plasmid sequences from draft genome assemblies using the clustered plasmid reference databases provided by MOB-cluster. It will also automatically provide the full typing information provided by MOB-typer. It optionally can use a chromosome depletion strategy based on closed genomes or user supplied filter of sequences to ignore.

### MOB-typer
Provides in silico predictions of the replicon family, relaxase type, mate-pair formation type and predicted transferability of the plasmid. Using a combination of biomarkers and MOB-cluster codes, it will also provide an observed host-range of your plasmid based on its replicon, relaxase and cluster assignment. This is combined with information mined from the literature to provide a prediction of the taxonomic rank at which the plasmid is likely to be stably maintained but it does not provide source attribution predictions.

## Installation ##

## Requires
+ Python >= 3.7
+ ete3 >= 3.1.3 (due to updated taxonomy database init)
+ pandas >= 0.22.0,<=1.05
+ biopython >= 1.80
+ pytables  >= 3.3
+ pycurl >= 7.43
+ numpy >= 1.11.1
+ scipy >= 1.1.0
+ six >= 1.10

## Dependencies

+ blast+ v. 2.3.0
+ mash v. 2.0


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

### Source
To build from source code directly on Ubuntu Linux distro, follow these commands that include Python
libraries and other dependencies install
```bash
apt update && apt install python3-pip #installs gcc compiler for pycurl
apt install libcurl4-openssl-dev libssl-dev #for pycurl
pip3 install Cython
apt install mash ncbi-blast+
python3 setup.py install && mob_init #to install and init databases
```

### Docker image
A docker images are also available at [https://hub.docker.com/r/kbessonov/mob_suite](https://hub.docker.com/r/kbessonov/mob_suite) and at [https://quay.io/repository/biocontainers/mob_suite](https://quay.io/repository/biocontainers/mob_suite?tab=tags)

```
% latest_tag=$(curl -H "Authorization: Bearer X" -X GET "https://quay.io/api/v1/repository/biocontainers/mob_suite/tag/" | jq .tags[].name | head -1 | sed -e 's|\"||g')
% docker pull quay.io/biocontainers/mob_suite:${latest_tag} 
% docker run --rm -v $(pwd):/mnt/ "kbessonov/mob_suite:${latest_tag}" mob_recon -i /mnt/assembly.fasta -t -o /mnt/mob_recon_output
```

### Singularity image
A singularity image could be built locally via Singularity recipe donated by Eric Deveaud or pulled from one of the repositories.

The recipe (`recipe.singularity`) is located in the `singularity` folder of this repository and installs MOB-Suite via `conda`. 

```bash
% singularity build mobsuite.simg recipe.singularity
```

As the simplest alternative, Singularity image can be pulled from [BioContainers repository](https://biocontainers.pro/tools/mob_suite) where `<version>` is
the desired version (e.g. `3.0.3--py_0`) or [Quay.io](https://quay.io/biocontainers/mob_suite) or [Docker Hub](https://hub.docker.com/r/kbessonov/mob_suite) repositories.

The MOB-Suite image (`mob_suite.sif`) will be generated using one of the following 3 methods.Next MOB-Suite tools could be run a on mounted directory via `--bind` like so `singularity run mob_suite.sif --bind $PWD:/mnt  mob_recon -i /mnt/<input_fasta>  -o /mnt/<output_directory>`

```bash
# Method 1
% singularity pull mob_suite.sif https://depot.galaxyproject.org/singularity/mob_suite:<version>
# Method 2
% singularity pull mob_suite.sif docker://kbessonov/mob_suite:3.0.3 #or for the latest version
# Method 3 - recommended
% latest_version=$(curl -H "Authorization: Bearer X" -X GET "https://quay.io/api/v1/repository/biocontainers/mob_suite/tag/" | jq .tags[].name | head -1 | sed -e 's|\"||g')
% singularity pull mob_suite.sif docker://quay.io/biocontainers/mob_suite:${latest_version}
```

## Using MOB-typer to perform replicon and relaxase typing of complete plasmids and to predict mobility and replicative plasmid host-range

### Setuptools
Clone this repository and install via setuptools. 

```
% git clone https://github.com/phac-nml/mob-suite.git
% cd mob-suite
% python setup.py install
```

## MGE detection
As of v. 3.1.0, MOB-recon and MOB-typer can report the blast HSP of the repetive mask file for IS/TN and other MGE elements in a new report file called mge.report.txt which will report blast hits from both chromosome and plasmid contigs. The MGE report is generated by default in MOB-recon and can be toggled on in MOB-typer by using the '--mge_report_file' parameter and specifying an output file.  This is a very naieve implementation of detecting MGE's and further work will improve the utility for users based on feedback.




## Using MOB-typer to perform replicon and relaxase typing of complete plasmids and predict mobility

You can perform plasmid typing using a fasta formated file containing a single plasmid represented by one or more contigs or it can treat all of the sequences in the fasta file as independent. The default behaviour is to treat all sequences in a file as from one plasmid, so do not include multiple unrelated plasmids in the file without specifying --multi as they will be treated as a single plasmid.


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

As of v. 3.0.0, we have added the ability of users to provide their own specific set of sequences to remove from plasmid reconstruction. This should be performed with caution and with the knowledge of your organism.  Filtering of sequences which are frequently of plasmid origin but are not in your organism is the primary use case we envision for this feature.

```
### User sequence mask
% mob_recon --infile assembly.fasta --outdir my_out_dir --filter_db filter.fasta
```

As of v. 3.0.0, we have provided the ability to use a collection of closed genomes which will be quickly checked using Mash for genomes which are genetically close and limit blast searches to those chromosomes. This more nuanced and automatic approach is recommended for users where there are sequences which should be filtered in one genomic context but not another. We provide as [an optional download](https://doi.org/10.5281/zenodo.3785351) a set of closed Enterobacteriacea genomes from NCBI which can be used to provide added accuracy for some organisms such as E. coli and Klebsiella where there are sequences which switch between chromosome and plasmids.
<br><br>
If reconstructed plasmids exceed the Mash distance for primary cluster assignment, then they will be assigned a name in the format novel_{md5} where the md5 hash is calculated based on all of the sequences belonging to that reconstructed plasmid. This will provide a unique name for the plasmids but any change will result in a corresponding change in the md5 hash. It is therefore not advised to use these assigned names for further analyses. Rather they should be highlighted as cases where targeted long read sequencing is required to obtain a closer database representative of that plasmid.

```
### Autodetected close genome filter
% mob_recon --infile assembly.fasta --outdir my_out_dir -g 2019-11-NCBI-Enterobacteriacea-Chromosomes.fasta
```
## Using MOB-cluster
Use this tool only to update the plasmid databases or build a new one, however MOB-cluster should only be run with closed high quality plasmids. If you add in poor quality data it can severely impact MOB-recon. As of v3.0.0, MOB-cluster has been re-written to utilize the output from MOB-typer to greatly speed up the process of updating and building plasmid databases by using pre-computed results. Clusters generated from earlier versions of MOB-suite are not compatible with the new clusters. We have provided a mapping file of previous cluster assignments and their new cluster accessions. Each cluster code is unique and will not be re-used.

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
| mge.report.txt | Blast HSP of detected MGE's/repetitive elements with contextual information |
| chromosome.fasta | Fasta file of all contigs found to belong to the chromosome |
| plasmid_(X).fasta | Each plasmid group is written to an individual fasta file which contains the assigned contigs |
| mobtyper_results | Aggregate MOB-typer report files for all identified plasmid |

# MOB-recon contig report format
| field  | Description |
| --------- |  --------- | 
| sample_id | Sample ID specified by user or default to filename |
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
| sample_id | Sample ID specified by user or default to filename |
| num_contigs | Number of sequences belonging to plasmid |
| size | Length in base pairs |
| gc | GC % |
| md5 | md5 hash |
| rep_type(s) | Replicon type(s) |
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
| sample_id | Sample ID specified by user or default to filename |
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

# MGE Report
| field  | Description |
| --------- |  --------- | 
| sample_id |  Sample ID specified by user or default to filename |
| molecule_type | Plasmid or Chromosome |
| primary_cluster_id | primary MOB-cluster id of neighbor |
| secondary_cluster_id | secondary MOB-cluster id of neighbor |
| contig_id | Sequence Identifier |
| size | Length in base pairs |
| gc | GC % |
| md5 | md5 hash |
| mge_id | Unique numeric id |
| mge_acs | GenBank Accession of MGE|
| mge_type | Primary type of the MGE |
| mge_subtype | Subtype of the MGE |
| mge_length | Length of the MGE query |
| mge_start | HSP start of MGE query |
| mge_end |  HSP end of MGE query |
| contig_start |  HSP start on contig  |
| contig_end | HSP end on contig |
| length | Length of HSP |
| sstrand | Stand of HSP |
| qcovhsp | Query coverage of HSP |
| pident | Sequence identity of HSP |
| evalue | Sequence evalue of HSP |
| bitscore | Sequence bitscore of HSP |

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

James Robertson - james.robertson@canada.ca <br>
Kyrylo Bessonov - kyrylo.bessonov@canada.ca

## License

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in
compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is
distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
implied. See the License for the specific language governing permissions and limitations under the
License.
