import os

MOB_CLUSTER_INFO_HEADER = [
    'id',
    'size',
    'gc',
    'md5',
    'organism',
    'primary_cluster_id',
    'primary_dist',
    'secondary_cluster_id',
    'secondary_dist',
    "rep_type(s)",
    "rep_type_accession(s)",
    "relaxase_type(s)",
    "relaxase_type_accession(s)"
]

MOB_TYPER_REPORT_HEADER = [
    'sample_id',
    'num_contigs',
    'size',
    'gc',
    'md5',
    'rep_type(s)',
    'rep_type_accession(s)',
    'relaxase_type(s)',
    'relaxase_type_accession(s)',
    'mpf_type',
    'mpf_type_accession(s)',
    'orit_type(s)',
    'orit_accession(s)',
    'predicted_mobility',
    'mash_nearest_neighbor',
    'mash_neighbor_distance',
    'mash_neighbor_identification',
    'primary_cluster_id',
    'secondary_cluster_id',
    'predicted_host_range_overall_rank',
    'predicted_host_range_overall_name',
    'observed_host_range_ncbi_rank',
    'observed_host_range_ncbi_name',
    'reported_host_range_lit_rank',
    'reported_host_range_lit_name',
    'associated_pmid(s)'
]

MOB_RECON_INFO_HEADER = [
    'sample_id',
    'molecule_type',
    'primary_cluster_id',
    'secondary_cluster_id',
    'contig_id',
    'size',
    'gc',
    'md5',
    'circularity_status',
    'rep_type(s)',
    'rep_type_accession(s)',
    'relaxase_type(s)',
    'relaxase_type_accession(s)',
    'mpf_type',
    'mpf_type_accession(s)',
    'orit_type(s)',
    'orit_accession(s)',
    'predicted_mobility',
    'mash_nearest_neighbor',
    'mash_neighbor_distance',
    'mash_neighbor_identification',
    'repetitive_dna_id',
    'repetitive_dna_type',
    'filtering_reason',

]

MGE_INFO_HEADER = [
    'sample_id',
    'molecule_type',
    'primary_cluster_id',
    'secondary_cluster_id',
    'contig_id',
    'size',
    'gc',
    'md5',
    'mge_id',
    'mge_acs',
    'mge_type',
    'mge_subtype',
    'mge_length',
    'mge_start',
    'mge_end',
    'contig_start',
    'contig_end',
    'length',
    'sstrand',
    'qcovhsp',
    'pident',
    'evalue',
    'bitscore'

]

NCBI_PLASMID_TAXONOMY_HEADER = [
    'sample_id',
    'num_contigs',
    'total_length',
    'gc',
    'md5',
    'rep_type(s)',
    'rep_type_accession(s)',
    'relaxase_type(s)',
    'relaxase_type_accession(s)',
    'mpf_type',
    'mpf_type_accession(s)',
    'orit_type(s)',
    'orit_accession(s)',
    'predicted_mobility',
    'mash_nearest_neighbor',
    'mash_neighbor_distance',
    'primary_cluster_id',
    'secondary_cluster_id',
    'organism',
    'taxid'
]

MOB_CLUSTER_INFO_HEADER = [
    'sample_id',
    'size',
    'gc',
    'md5',
    'organism',
    'taxid',
    'rep_type(s)',
    'rep_type_accession(s)',
    'relaxase_type(s)',
    'relaxase_type_accession(s)',
    'mpf_type',
    'mpf_type_accession(s)',
    'orit_type(s)',
    'orit_accession(s)',
    'predicted_mobility',
    'primary_cluster_id',
    'primary_dist',
    'secondary_cluster_id',
    'secondary_dist',
]

ACS_LETTER_VALUES = {
    'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'J': 9,
    'K': 10, 'L': 11, 'M': 12, 'N': 13, 'O': 14, 'P': 15, 'Q': 16, 'R': 17, 'S': 18,
    'T': 19, 'U': 20, 'V': 21, 'W': 22, 'X': 23, 'Y': 24, 'Z': 25
}

ACS_VALUES_TO_LETTERS = {
    0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G', 7: 'H', 8: 'I', 9: 'J',
    10: 'K', 11: 'L', 12: 'M', 13: 'N', 14: 'O', 15: 'P', 16: 'Q', 17: 'R', 18: 'S',
    19: 'T', 20: 'U', 21: 'V', 22: 'W', 23: 'X', 24: 'Y', 25: 'Z'
}

MAX_ACS_VALUE = 650999
ACS_FORMAT_VALUES = [26000, 1000, 1]

LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'

LIT_PLASMID_TAXONOMY_HEADER = [
    'sample_id', 'rep_type(s)', 'length', 'host_species', 'host_taxid', 'reported_host_range_taxid', 'pmid',
    'pmcid', 'doi', 'year', 'author', 'notes'
]



default_database_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'databases')

LOG_FORMAT = '%(asctime)s %(name)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
