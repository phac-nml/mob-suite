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

LIT_PLASMID_TAXONOMY_HEADER = [
                                'sample_id', 'rep_type(s)', 'length', 'host_species', 'host_taxid', 'reported_host_range_taxid', 'pmid',
                                'pmcid', 'doi', 'year', 'author', 'notes'
]

LIT_PLASMID_TAXONOMY_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),"databases/host_range_literature_plasmidDB.txt")




default_database_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'databases')

ETE3DBTAXAFILE = os.path.abspath(default_database_dir + "/taxa.sqlite")

NCBI_PLASMID_TAXONOMY_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),"databases/host_range_ncbirefseq_plasmidDB.txt")

ETE3_LOCK_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),"databases/ETE3_DB.lock")

LOG_FORMAT = '%(asctime)s %(name)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'