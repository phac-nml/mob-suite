#!/usr/bin/env python3

import logging
import os, re, pandas, collections, shutil, sys, time, tempfile
from argparse import (ArgumentParser, FileType)
from collections import OrderedDict,Counter
from operator import itemgetter
from ete3 import NCBITaxa
from mob_suite.version import __version__
import mob_suite.mob_init
from mob_suite.wrappers import mash
from mob_suite.utils import \
    replicon_blast, \
    mob_blast, \
    getRepliconContigs, \
    fix_fasta_header, \
    getMashBestHit, \
    calcFastaStats, \
    verify_init, \
    check_dependencies, \
    read_sequence_info, \
    read_file_to_dict,\
    default_database_dir, \
    dict_from_alt_key_list, \
    get_data_associated_with_key, \
    writeReport, \
    getMashBestHitMultiSeq, \
    calcFastaStatsIndividual


MOB_TYPER_REPORT_HEADER = ['sample_id', 'num_contigs', 'total_length', 'gc', 'md5','rep_type(s)',
                           'rep_type_accession(s)', 'relaxase_type(s)', 'relaxase_type_accession(s)',
                           'mpf_type', 'mpf_type_accession(s)', 'orit_type(s)', 'orit_accession(s)',
                           'predicted_mobility', 'mash_nearest_neighbor', 'mash_neighbor_distance',
                           'primary_cluster_id', 'secondary_cluster_id', 'predicted_host_range_overall_rank',
                           'predicted_host_range_overall_name', 'observed_host_range_ncbi_rank',
                           'observed_host_range_ncbi_name', 'reported_host_range_lit_rank',
                           'reported_host_range_lit_name', 'associated_pmid(s)' ]

database_directory = os.path.abspath(default_database_dir)
ETE3DBTAXAFILE = os.path.abspath(database_directory + "/taxa.sqlite")


NCBI_PLASMID_TAXONOMY_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),"databases/host_range_ncbirefseq_plasmidDB.txt")
NCBI_PLASMID_TAXONOMY_HEADER = [
                                'sample_id', 'num_contigs', 'total_length', 'gc', 'md5','rep_type(s)',
                                'rep_type_accession(s)', 'relaxase_type(s)', 'relaxase_type_accession(s)',
                                'mpf_type', 'mpf_type_accession(s)', 'orit_type(s)', 'orit_accession(s)',
                                'predicted_mobility', 'mash_nearest_neighbor', 'mash_neighbor_distance',
                                'primary_cluster_id', 'secondary_cluster_id', 'organism', 'taxid'
]

LIT_PLASMID_TAXONOMY_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),"databases/host_range_literature_plasmidDB.txt")
LIT_PLASMID_TAXONOMY_HEADER = [
                                'sample_id', 'rep_type(s)', 'length', 'host_species', 'host_taxid', 'reported_host_range_taxid', 'pmid',
                                'pmcid', 'doi', 'year', 'author', 'notes'
]




def init_console_logger(lvl=2):
    root = logging.getLogger()

    LOG_FORMAT = '%(asctime)s %(name)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]

    report_lvl = logging_levels[lvl]
    root.setLevel(report_lvl)  # set root logger level

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)
    return logging.getLogger(__name__)


def parse_args():
    "Parse the input arguments, use '-h' for help"
    default_database_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'databases')
    parser = ArgumentParser(
        description="MOB-Typer: Plasmid typing version: {}".format(
            __version__))

    parser.add_argument('-o', '--out_file', type=str, required=True, help='Output file to write results')
    parser.add_argument('-a', '--analysis_dir', type=str, required=False, help='Analysis directory')
    parser.add_argument('-i', '--infile', type=str, required=True, help='Input assembly fasta file to process')
    parser.add_argument('-s', '--sample_id', type=str, required=False, help='Sample name for report')
    parser.add_argument('-x', '--multi', required=False, help='Fasta file is multiple independant records',action='store_true',
                        default=False)

    parser.add_argument('-n', '--num_threads', type=int, required=False, help='Number of threads to be used', default=1)

    parser.add_argument('--min_rep_evalue', type=str, required=False,
                        help='Minimum evalue threshold for replicon blastn',
                        default=0.00001)

    parser.add_argument('--min_mob_evalue', type=str, required=False,
                        help='Minimum evalue threshold for relaxase tblastn',
                        default=0.00001)

    parser.add_argument('--min_con_evalue', type=str, required=False, help='Minimum evalue threshold for contig blastn',
                        default=0.00001)

    parser.add_argument('--min_ori_evalue', type=str, required=False,
                        help='Minimum evalue threshold for oriT elements blastn',
                        default=0.00001)
    parser.add_argument('--min_mpf_evalue', type=str, required=False,
                        help='Minimum evalue threshold for mpf elements blastn',
                        default=0.00001)

    parser.add_argument('--min_rep_ident', type=int, required=False, help='Minimum sequence identity for replicons',
                        default=80)

    parser.add_argument('--min_mob_ident', type=int, required=False, help='Minimum sequence identity for relaxases',
                        default=80)

    parser.add_argument('--min_ori_ident', type=int, required=False,
                        help='Minimum sequence identity for oriT elements', default=90)
    parser.add_argument('--min_mpf_ident', type=int, required=False,
                        help='Minimum sequence identity for mpf elements', default=80)

    parser.add_argument('--min_rep_cov', type=int, required=False,
                        help='Minimum percentage coverage of replicon query by input assembly',
                        default=80)

    parser.add_argument('--min_mob_cov', type=int, required=False,
                        help='Minimum percentage coverage of relaxase query by input assembly',
                        default=80)

    parser.add_argument('--min_ori_cov', type=int, required=False,
                        help='Minimum percentage coverage of oriT',
                        default=90)
    parser.add_argument('--min_mpf_cov', type=int, required=False,
                        help='Minimum percentage coverage of mpf',
                        default=80)

    parser.add_argument('--min_overlap', type=int, required=False,
                        help='Minimum overlap of fragments',
                        default=10)

    parser.add_argument('--keep_tmp', required=False,help='Do not delete temporary file directory', action='store_true')
    parser.add_argument('--debug', required=False, help='Show debug information', action='store_true')
    parser.add_argument('-m','--plasmid_meta', type=str, required=False,
                        help='MOB-cluster plasmid cluster formatted file matched to the reference plasmid db',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/clusters.txt'))
    parser.add_argument('--plasmid_mash_db', type=str, required=False,
                        help='Companion Mash database of reference database',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/ncbi_plasmid_full_seqs.fas.msh'))
    parser.add_argument('--plasmid_replicons', type=str, required=False, help='Fasta of plasmid replicons',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/rep.dna.fas'))
    parser.add_argument('--plasmid_mob', type=str, required=False, help='Fasta of plasmid relaxases',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/mob.proteins.faa'))
    parser.add_argument('--plasmid_mpf', type=str, required=False, help='Fasta of known plasmid mate-pair proteins',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/mpf.proteins.faa'))
    parser.add_argument('--plasmid_orit', type=str, required=False, help='Fasta of known plasmid oriT dna sequences',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'databases/orit.fas'))
    parser.add_argument('--primary_cluster_dist', type=int, required=False, help='Mash distance for assigning primary cluster id 0 - 1', default=0.06)
    parser.add_argument('--secondary_cluster_dist', type=int, required=False, help='Mash distance for assigning primary cluster id 0 - 1',
                        default=0.025)
    parser.add_argument('-d','--database_directory',default=default_database_dir,
                        required=False,
                        help='Directory you want to use for your databases. If the databases are not already '
                             'downloaded, they will be downloaded automatically. Defaults to {}. '
                             'If you change this from the default, will override --plasmid_mash_db, '
                             '--plasmid_replicons, --plasmid_mob, --plasmid_mpf, and '
                             '--plasmid_orit'.format(default_database_dir))

    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    return parser.parse_args()


def isETE3DBTAXAFILEexists():
    if not os.path.exists(ETE3DBTAXAFILE):
        return False
    else:
        return True
def initETE3Database():
    lockfilepath = os.path.join(database_directory, ".lock")

    if os.path.exists(lockfilepath) == False:
        open(file=lockfilepath, mode="w").close()
        logging.info("Placed lock file at {}".format(lockfilepath))
    else:
        while os.path.exists(lockfilepath):
            elapsed_time = time.time() - os.path.getmtime(lockfilepath)
            logging.info("Lock file found at {}. Waiting for other processes to finish ete3 database init ...".format(lockfilepath))
            logging.info("Elapsed time {} min. Will continue processing after 16 min mark.".format(int(elapsed_time/60)))
            if elapsed_time >= 1000:
                logging.info("Elapsed time {} min. Assuming previous process completed all init steps. Continue ...".format(int(elapsed_time/60)))
                try: #if previous process failed, no processes are running and > 16 min passed since the lock was created
                    os.remove(lockfilepath)
                except: #continue if file was removed by other process
                    pass
                break
            time.sleep(60) #recheck every 1 min if lock file was removed by other process
        logging.info("Lock file no longer exists. Assuming init process completed successfully")

    ncbi = NCBITaxa()
    ncbi.dbfile = ETE3DBTAXAFILE
    ncbi.update_taxonomy_database()

    try:
        os.remove(lockfilepath)
        logging.info("Lock file removed.")
    except:
        logging.warning("Lock file is already removed by some other process.")
        pass

    try:
        os.remove(os.path.join(os.getcwd(), "taxdump.tar.gz"))
        logging.info("Removed residual taxdump.tar.gz as ete3 is not doing proper cleaning job.")
    except:
        pass
    logging.info("ETE3 database init completed successfully.")

def determine_mpf_type(hits):
    types = dict()
    for hit in hits:
        if not hit in types:
            types[hit] = 0
        types[hit] += 1

    return max(types, key=lambda i: types[i])

def getAssocValues(query_list_values,look_up_key,value_key,data):
    values = []

    for q in query_list_values:
        values.extend(get_data_associated_with_key(look_up_key, q, value_key, data))
    return  values

def initMOBTyperReportTemplate(header):
    data = {}
    for i in header:
        data[i] = '-'
    return data

def getHeirarchy(taxid):
    if not isETE3DBTAXAFILEexists():
        logging.info("Did not find taxa.sqlite in {}. Initializaing ete3 taxonomy database".format(ETE3DBTAXAFILE))
        initETE3Database()

    ncbi = NCBITaxa(dbfile=ETE3DBTAXAFILE)
    if not isETE3DBTAXAFILEexists():
        logging.error("Tried ete3 init, but still was not able to find taxa.sqlite file for ete3 lib in {}. Aborting".format(ETE3DBTAXAFILE))
        return ['-','-']

    if not isinstance(taxid,int):
        return {'names': [], 'ranks': []}

    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    ranks = []
    for id in lineage:
        ranks.append(ncbi.get_rank(id))

    return {'names': names, 'ranks': names}




def filter_invalid_taxids(taxids):
    filtered = []
    for i in taxids:
        if i == 'NaN' or i == 'nan' or i == None or i == '-' :
            continue

        try:
            filtered.append(int(i))

        except ValueError:
            continue

    return filtered


def getHeirarchy(taxid):
    if not isETE3DBTAXAFILEexists():
        logging.info("Did not find taxa.sqlite in {}. Initializaing ete3 taxonomy database".format(ETE3DBTAXAFILE))
        initETE3Database()

    ncbi = NCBITaxa(dbfile=ETE3DBTAXAFILE)
    if not isETE3DBTAXAFILEexists():
        logging.error("Tried ete3 init, but still was not able to find taxa.sqlite file for ete3 lib in {}. Aborting".format(ETE3DBTAXAFILE))
        return ['-','-']

    if not isinstance(taxid,int):
        return {'names': [], 'ranks': []}

    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    ranks = []
    for id in lineage:
        ranks.append(ncbi.get_rank(id))

    return {'names': names, 'ranks': names}


def getTaxid(taxon):
    if not isETE3DBTAXAFILEexists():
        logging.info("Did not find taxa.sqlite in {}. Initializaing ete3 taxonomy database".format(ETE3DBTAXAFILE))
        initETE3Database()

    ncbi = NCBITaxa(dbfile=ETE3DBTAXAFILE)
    if not isETE3DBTAXAFILEexists():
        logging.error("Tried ete3 init, but still was not able to find taxa.sqlite file for ete3 lib in {}. Aborting".format(ETE3DBTAXAFILE))
        return ['-','-']

    taxid = ncbi.get_name_translator(taxon)

    return  taxid


def getTaxonConvergence(taxids):
    if not isETE3DBTAXAFILEexists():
        logging.info("Did not find taxa.sqlite in {}. Initializaing ete3 taxonomy database".format(ETE3DBTAXAFILE))
        initETE3Database()

    ncbi = NCBITaxa(dbfile=ETE3DBTAXAFILE)
    if not isETE3DBTAXAFILEexists():
        logging.error(
            "Tried ete3 init, but still was not able to find taxa.sqlite file for ete3 lib in {}. Aborting".format(
                ETE3DBTAXAFILE))
        return ['-', '-']

    tax_ranks = ['genus','family','order','class','phylum',  'species']

    phylla = []
    lineages = []
    if len(taxids) == 0:
        return ['-', '-']

    elif len(taxids) == 1:
        common_taxids = ncbi.get_lineage(taxids[0])

    else:
        for taxid in taxids:
            lineage = ncbi.get_lineage(taxid)
            ranks = ncbi.get_rank(lineage)
            names = ncbi.get_taxid_translator(lineage)
            for id in ranks:
                rank = ranks[id]
                name = names[id]
                if rank == 'phylum':
                    phylla.append(name)

            lineages.append(lineage)

        common_taxids = lineages[0]
        for i in range(1,len(lineages)):
            common_taxids = list(set(common_taxids) - set(set(common_taxids) - set(lineages[i])))



    if len(common_taxids) == 0:
        return ['-', '-']

    taxa = {}

    for i in range(0,len(common_taxids)):
        rank = ncbi.get_rank([common_taxids[i]])[common_taxids[i]]
        name = ncbi.get_taxid_translator([common_taxids[i]])[common_taxids[i]]

        if not rank in tax_ranks:
            continue

        taxa[rank] = name

    for rank in tax_ranks:
        if rank in taxa:
            return [rank,taxa[rank]]

    phylla = list(set(phylla))
    sorted(phylla)
    if len(phylla) > 1:
        return ('multi-phylla',",".join(phylla))
    elif len(phylla) == 1:
        return ('phylum', ",".join(phylla))

    return(['-','-'])




def hostrange(replion_types, relaxase_types, mob_cluster_id,ncbi,lit):

    host_range_predictions = {
        'observed_host_range_ncbi_name':'',
        'observed_host_range_ncbi_rank':'',
        'reported_host_range_lit_name':'',
        'reported_host_range_lit_rank':'',
        'predicted_host_range_overall_name':'',
        'predicted_host_range_overall_rank':'',
        'associated_pmid(s)':'-'

    }


    #Determine taxids associated with NCBI records
    if mob_cluster_id != '-' :
        ncbi_cluster_taxids = getAssocValues([mob_cluster_id], 'primary_cluster_id', 'taxid', ncbi)

    else:
        ncbi_cluster_taxids = []

    if len(replion_types) > 0 :
        ncbi_replicon_taxids = getAssocValues(replion_types,'rep_type(s)','taxid',ncbi)
        lit_replicon_taxids = getAssocValues(replion_types, 'rep_type(s)', 'reported_host_range_taxid', lit)
        lit_replicon_taxids = lit_replicon_taxids + getAssocValues(replion_types, 'rep_type(s)', 'host_taxid', lit)
        pmids = list(set(getAssocValues(replion_types, 'rep_type(s)', 'pmid', lit)))
        sorted(pmids)
        host_range_predictions['associated_pmid(s)'] = '; '.join(str(filter_invalid_taxids(x)) for x in pmids)

    else:
        ncbi_replicon_taxids = []
        lit_replicon_taxids = []

    if len(relaxase_types) > 0:
        ncbi_relaxase_taxids = getAssocValues(relaxase_types, 'relaxase_type_accession(s)', 'taxid', ncbi)
    else:
        ncbi_relaxase_taxids = []

    if  len(ncbi_cluster_taxids) == 0 and len(ncbi_replicon_taxids) == 0 and len(ncbi_relaxase_taxids) == 0:
        host_range_predictions['observed_host_range_ncbi_rank'] = '-'
        host_range_predictions['observed_host_range_ncbi_name'] = '-'
        ncbi_unique_taxids = []
    else:
        ncbi_unique_taxids = filter_invalid_taxids(list(set(ncbi_replicon_taxids + ncbi_cluster_taxids + ncbi_relaxase_taxids)))
        host_range_predictions['observed_host_range_ncbi_rank'], host_range_predictions['observed_host_range_ncbi_name'] = getTaxonConvergence(ncbi_unique_taxids)


    # Determine taxids associated with literature

    lit_unique_taxids = filter_invalid_taxids(list(set(lit_replicon_taxids )))

    host_range_predictions['reported_host_range_lit_rank'], host_range_predictions['reported_host_range_lit_name'] = getTaxonConvergence(lit_unique_taxids)

    #determine overall host range
    overall_taxids = filter_invalid_taxids(list(set(ncbi_unique_taxids + lit_unique_taxids)))
    host_range_predictions['predicted_host_range_overall_rank'], host_range_predictions['predicted_host_range_overall_name'] = getTaxonConvergence(overall_taxids)

    #move host-range prediction up to family when it is at genus or species level
    if host_range_predictions['predicted_host_range_overall_rank'] == 'genus' or host_range_predictions['predicted_host_range_overall_rank'] == 'species':
        taxid = getTaxid(host_range_predictions['predicted_host_range_overall_name'])
        heir = getHeirarchy(taxid)
        names = heir['names']
        ranks = heir['ranks']

        for i in range(0,len(ranks)):
            if ranks[i] == 'Family':
                host_range_predictions['predicted_host_range_overall_rank'] = 'family'
                host_range_predictions['predicted_host_range_overall_name'] = names[i]


    return host_range_predictions




def getBioMarkerContigs(contig_dict):

    biomarkers = {}

    for contig_id in contig_dict:
        if not contig_id in biomarkers:
            biomarkers[contig_id] = {'acs':[],'types':[]}
        hits = contig_dict[contig_id]

        if len(hits)== 0:
            continue

        for hit in hits:
            info = hit.split("|")
            biomarkers[contig_id]['acs'].append(info[0])
            biomarkers[contig_id]['types'].append(info[1])

    return  biomarkers

def sort_biomarkers(biomarker_dict):

    for id in biomarker_dict:
        acs = biomarker_dict[id]['acs']
        types = biomarker_dict[id]['types']

        if len(acs) == 0:
            continue

        tmp_dict = {}

        for i in range(0,len(acs)):
            tmp_dict[acs[i]] = types[i]

        tmp_dict = OrderedDict(sorted(tmp_dict.items(), key=itemgetter(1), reverse=False))

        biomarker_dict[id]['acs'] = list(tmp_dict.keys())
        biomarker_dict[id]['types'] = list(tmp_dict.values())

    return  biomarker_dict



def mergeBiomarkers(biomarkers,sample_id):
    acs = []
    types = []
    tmp_dict = {}
    for contig_id in biomarkers:
        acs = acs + biomarkers[contig_id]['acs']
        types = types + biomarkers[contig_id]['types']

    tmp_dict[sample_id] = {}
    tmp_dict[sample_id]['acs'] = acs
    tmp_dict[sample_id]['types'] = types

    return  tmp_dict

def predict_mobility(biomarker_dict):


    for sample_id in biomarker_dict:
        mobility = 'non-mobilizable'
        biomarkers = biomarker_dict[sample_id]
        if 'found_mob' in biomarkers and 'found_mpf' and 'found_orit':
            if len(biomarkers['found_mob']['acs']) > 0:
                if len(biomarkers['found_mpf']['acs']) > 0:
                    mobility = 'conjugative'
                else:
                    mobility = 'mobilizable'
            elif len(biomarkers['found_orit']['acs']) > 0:
                mobility = 'mobilizable'
        biomarker_dict[sample_id]['predicted_mobility'] = mobility

    return biomarker_dict







def main():
    default_database_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'databases')
    args = parse_args()

    if args.debug:
       logger = init_console_logger(3)
    else:
       logger = init_console_logger(2)

    logger.info('Running Mob-typer version {}'.format(__version__))


    if not os.path.isfile(args.infile):
        logger.info('Error, fasta file does not exist')
        sys.exit()

    if not args.analysis_dir:
        tmp_dir = tempfile.TemporaryDirectory( dir = tempfile.gettempdir()).name
    else:
        tmp_dir = args.analysis_dir

    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir, 0o755)

    if not isinstance(args.num_threads, int):
        logger.info('Error number of threads must be an integer, you specified "{}"'.format(args.num_threads))

    database_dir = os.path.abspath(args.database_directory)

    if args.sample_id is None:
        sample_id = re.sub(r"\.(fasta|fa|fas){1,1}", "", os.path.basename(args.infile))
    else:
        sample_id = args.sample_id

    verify_init(logger,database_dir)

    # Script arguments
    input_fasta = args.infile
    report_file = args.out_file
    num_threads = int(args.num_threads)
    keep_tmp = args.keep_tmp

    if args.multi:
        multi = True
    else:
        multi = False


    if not (args.primary_cluster_dist >= 0 and args.primary_cluster_dist <= 1):
        logging.error('Error distance thresholds must be between 0 - 1: {}'.format(args.primary_cluster_dist))
        sys.exit()
    else:
        primary_distance = float(args.primary_cluster_dist)

    if not (args.secondary_cluster_dist >= 0 and args.secondary_cluster_dist <= 1):
        logging.error('Error distance thresholds must be between 0 - 1: {}'.format(args.secondary_cluster_dist))
        sys.exit()
    else:
        secondary_distance = float(args.secondary_cluster_dist)

    if database_dir == default_database_dir:
        mob_ref = args.plasmid_mob
        mpf_ref = args.plasmid_mpf
        orit_ref = args.plasmid_orit
        mash_db = args.plasmid_mash_db
        replicon_ref = args.plasmid_replicons
        plasmid_meta = args.plasmid_meta
    else:
        mob_ref = os.path.join(database_dir, 'mob.proteins.faa')
        mpf_ref = os.path.join(database_dir, 'mpf.proteins.faa')
        orit_ref = os.path.join(database_dir, 'orit.fas')
        mash_db = os.path.join(database_dir, 'ncbi_plasmid_full_seqs.fas.msh')
        replicon_ref = os.path.join(database_dir, 'rep.dna.fas')
        plasmid_meta = os.path.join(database_dir, 'clusters.txt')


    fixed_fasta = os.path.join(tmp_dir, 'fixed.input.fasta')
    replicon_blast_results = os.path.join(tmp_dir, 'replicon_blast_results.txt')
    mob_blast_results = os.path.join(tmp_dir, 'mobtyper_blast_results.txt')
    mpf_blast_results = os.path.join(tmp_dir, 'mpf_blast_results.txt')
    orit_blast_results = os.path.join(tmp_dir, 'orit_blast_results.txt')

    if os.path.isfile(mob_blast_results):
        os.remove(mob_blast_results)
    if os.path.isfile(mpf_blast_results):
        os.remove(mpf_blast_results)
    if os.path.isfile(orit_blast_results):
        os.remove(orit_blast_results)
    if os.path.isfile(replicon_blast_results):
        os.remove(replicon_blast_results)


    mash_file = os.path.join(tmp_dir, "mash_{}.txt".format(sample_id))


    # Input numeric params

    min_rep_ident = float(args.min_rep_ident)
    min_mob_ident = float(args.min_mob_ident)
    min_ori_ident = float(args.min_ori_ident)
    min_mpf_ident = float(args.min_mpf_ident)

    idents = {'min_rep_ident': min_rep_ident, 'min_mob_ident': min_mob_ident, 'min_ori_ident': min_ori_ident}

    for param in idents:

        value = float(idents[param])

        if value < 60:
            logger.error("Error: {} is too low, please specify an integer between 70 - 100".format(param))
            sys.exit(-1)
        if value > 100:
            logger.error("Error: {} is too high, please specify an integer between 70 - 100".format(param))
            sys.exit(-1)


    min_rep_cov = float(args.min_rep_cov)
    min_mob_cov = float(args.min_mob_cov)
    min_ori_cov = float(args.min_ori_cov)
    min_mpf_cov = float(args.min_mpf_cov)



    covs = {'min_rep_cov': min_rep_cov, 'min_mob_cov': min_mob_cov, 'min_con_cov': min_ori_cov,
            'min_rpp_cov': min_ori_cov}

    for param in covs:

        value = float(covs[param])

        if value < 60:
            logger.error("Error: {} is too low, please specify an integer between 50 - 100".format(param))
            sys.exit(-1)
        if value > 100:
            logger.error("Error: {} is too high, please specify an integer between 50 - 100".format(param))
            sys.exit(-1)


    min_rep_evalue = float(args.min_rep_evalue)
    min_mob_evalue = float(args.min_mob_evalue)
    min_ori_evalue = float(args.min_ori_evalue)
    min_mpf_evalue = float(args.min_mpf_evalue)


    evalues = {'min_rep_evalue': min_rep_evalue, 'min_mob_evalue': min_mob_evalue, 'min_con_evalue': min_ori_evalue}

    for param in evalues:

        value = float(evalues[param])

        if value > 1:
            logger.error("Error: {} is too high, please specify an float evalue between 0 to 1".format(param))
            sys.exit(-1)


    check_dependencies(logger)

    needed_dbs = [replicon_ref, mob_ref, mash_db, mpf_ref]

    for db in needed_dbs:
        if (not os.path.isfile(db)):
            logger.info('Warning! Needed database missing "{}"'.format(db))
            mob_suite.mob_init.main()


    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir, 0o755)

    mobtyper_results = {}
    fix_fasta_header(input_fasta, fixed_fasta)
    if multi:
        fastaSeqStats = calcFastaStatsIndividual(fixed_fasta)
        for sid in fastaSeqStats:
            if not sid in mobtyper_results:
                mobtyper_results[sid] = initMOBTyperReportTemplate(MOB_TYPER_REPORT_HEADER)
                mobtyper_results[sid]['sample_id'] = sid
                mobtyper_results[sid]['num_contigs'] = 1
                for stat in fastaSeqStats[sid]:
                    mobtyper_results[sid][stat] = fastaSeqStats[sid][stat]
    else:
        mobtyper_results[sample_id] = initMOBTyperReportTemplate(MOB_TYPER_REPORT_HEADER)
        mobtyper_results[sample_id]['sample_id'] = sample_id
        fastaSeqStats = calcFastaStats(fixed_fasta)
        mobtyper_results[sample_id]['md5'] = fastaSeqStats['md5']
        mobtyper_results[sample_id]['total_length'] = fastaSeqStats['size']
        mobtyper_results[sample_id]['num_contigs'] = fastaSeqStats['num_seq']
        mobtyper_results[sample_id]['gc'] = fastaSeqStats['gc_content']


    #Get cluster information
    reference_sequence_meta = read_sequence_info(plasmid_meta)


    # run individual marker blasts
    logger.info('Running replicon blast on {}'.format(replicon_ref))
    replicon_contigs = getRepliconContigs(
        replicon_blast(replicon_ref, fixed_fasta, min_rep_ident, min_rep_cov, min_rep_evalue, tmp_dir, replicon_blast_results,
                       num_threads=num_threads))
    logger.info('Running relaxase blast on {}'.format(mob_ref))
    mob_contigs = getRepliconContigs(
        mob_blast(mob_ref, fixed_fasta, min_mob_ident, min_mob_cov, min_mob_evalue, tmp_dir, mob_blast_results, num_threads=num_threads))

    logger.info('Running mpf blast on {}'.format(mob_ref))
    mpf_contigs = getRepliconContigs(
        mob_blast(mpf_ref, fixed_fasta, min_mpf_ident, min_mpf_cov, min_mpf_evalue, tmp_dir, mpf_blast_results, num_threads=num_threads))

    logger.info('Running orit blast on {}'.format(replicon_ref))
    orit_contigs = getRepliconContigs(
        replicon_blast(orit_ref, fixed_fasta, min_ori_ident, min_ori_cov, min_ori_evalue, tmp_dir, orit_blast_results,
                       num_threads=num_threads))

    # Get closest neighbor by mash distance in the entire plasmid database

    m = mash()
    mashfile_handle = open(mash_file, 'w')
    mash_top_hit = {}

    if multi:
        tmp_mash_db_file = os.path.join(tmp_dir,"tmp_mash.msh")
        m.mashsketch(fixed_fasta, tmp_mash_db_file , num_threads=num_threads)
        m.run_mash(mash_db, tmp_mash_db_file , mashfile_handle)
        mash_results = m.read_mash(mash_file)
        mash_top_hit = getMashBestHitMultiSeq(mash_results)

    else:

        m.run_mash(mash_db, fixed_fasta, mashfile_handle)
        mash_results = m.read_mash(mash_file)
        mash_top_hit[sample_id] = getMashBestHit(mash_results)


    #Get MOB-cluster information for reference database matches
    mob_cluster_memberships = {}

    for s in mash_top_hit:

        if not s in mob_cluster_memberships:
            mob_cluster_memberships[s] = {'primary_cluster_id':'-','secondary_cluster_id':'-'}

        plasmid_primary_acs = '-'
        plasmid_secondary_acs = '-'

        if mash_top_hit[s]['top_hit'] in reference_sequence_meta:
            if mash_top_hit[s]['mash_hit_score'] <= primary_distance:
                plasmid_primary_acs = reference_sequence_meta[mash_top_hit[s]['top_hit']]['primary_cluster_id']
            if mash_top_hit[s]['mash_hit_score'] <= secondary_distance:
                plasmid_secondary_acs = reference_sequence_meta[mash_top_hit[s]['top_hit']]['secondary_cluster_id']

            mob_cluster_memberships[s]['primary_cluster_id'] = plasmid_primary_acs
            mob_cluster_memberships[s]['secondary_cluster_id'] = plasmid_secondary_acs
            mob_cluster_memberships[s]['mash_nearest_neighbor'] = mash_top_hit[s]['top_hit']
            mob_cluster_memberships[s]['mash_neighbor_distance'] = mash_top_hit[s]['mash_hit_score']

        else:
            logging.warning("Found sequence: {} without a cluster record, confirm cluster file matches mash file".format(mash_top_hit[s]['top_hit']))



    found_replicons = sort_biomarkers(getBioMarkerContigs(replicon_contigs))
    found_mob = sort_biomarkers(getBioMarkerContigs(mob_contigs))
    found_mpf = sort_biomarkers(getBioMarkerContigs(mpf_contigs))
    found_orit = sort_biomarkers(getBioMarkerContigs(orit_contigs))

    #merge all of the biomarkers if all contig sequences should be treated beloning to one plasmid
    biomarkers = {}

    if not multi:
        biomarkers[sample_id] = {'found_replicons':{'acs':[],'types':[]},
                                 'found_mob':{'acs':[],'types':[]},
                                 'found_mpf':{'acs':[],'types':[]},
                                 'found_orit':{'acs':[],'types':[]}}

        tmp = sort_biomarkers(mergeBiomarkers(found_replicons, sample_id))
        biomarkers[sample_id]['found_replicons']['acs'] = tmp[sample_id]['acs']
        biomarkers[sample_id]['found_replicons']['types'] = tmp[sample_id]['types']

        tmp = sort_biomarkers(mergeBiomarkers(found_mob, sample_id))
        biomarkers[sample_id]['found_mob']['acs'] = tmp[sample_id]['acs']
        biomarkers[sample_id]['found_mob']['types'] = tmp[sample_id]['types']

        tmp = sort_biomarkers(mergeBiomarkers(found_mpf, sample_id))
        biomarkers[sample_id]['found_mpf']['acs'] =  tmp[sample_id]['acs']
        biomarkers[sample_id]['found_mpf']['types'] = tmp[sample_id]['types']

        tmp = sort_biomarkers(mergeBiomarkers(found_orit, sample_id))
        biomarkers[sample_id]['found_orit']['acs'] = tmp[sample_id]['acs']
        biomarkers[sample_id]['found_orit']['types'] = tmp[sample_id]['types']
    else:
        for s in mob_cluster_memberships:
            if s not in biomarkers:
                biomarkers[s] = {'found_replicons':{'acs':[],'types':[]},
                                 'found_mob':{'acs':[],'types':[]},
                                 'found_mpf':{'acs':[],'types':[]},
                                 'found_orit':{'acs':[],'types':[]}}
            if s in found_replicons:
                biomarkers[s]['found_replicons'] = found_replicons[s]
            if s in found_mob:
                biomarkers[s]['found_mob'] = found_mob[s]
            if s in found_mpf:
                biomarkers[s]['found_mpf'] = found_mpf[s]
            if s in found_orit:
                biomarkers[s]['found_orit'] = found_orit[s]

    biomarkers = predict_mobility(biomarkers)


    ncbi = dict_from_alt_key_list(
        read_file_to_dict(NCBI_PLASMID_TAXONOMY_FILE, NCBI_PLASMID_TAXONOMY_HEADER, separater="\t"), "sample_id")
    lit = dict_from_alt_key_list(read_file_to_dict(LIT_PLASMID_TAXONOMY_FILE, LIT_PLASMID_TAXONOMY_HEADER, separater="\t"),"sample_id")

    for sid in biomarkers:
        rep_types = biomarkers[sid]['found_replicons']['types']
        relaxase_types = biomarkers[sid]['found_mob']['acs']
        mob_cluster_id = mob_cluster_memberships[sid]['primary_cluster_id']
        host_range = hostrange(rep_types, relaxase_types, mob_cluster_id,ncbi,lit)
        for key in host_range:
            mobtyper_results[sid][key] = host_range[key]

        for key in biomarkers[sid]:
            if key in mobtyper_results[sid]:
                mobtyper_results[sid][key] = biomarkers[sid][key]

        if sid in mob_cluster_memberships:
            for key in mob_cluster_memberships[sid]:
                if key in mobtyper_results[sid]:
                    mobtyper_results[sid][key] = mob_cluster_memberships[sid][key]



        if len(biomarkers[sid]['found_replicons']['types']) > 0:
            mobtyper_results[sid]['rep_type(s)'] = ",".join(biomarkers[sid]['found_replicons']['types'])
            mobtyper_results[sid]['rep_type_accession(s)'] = ",".join(biomarkers[sid]['found_replicons']['acs'])

        if len(biomarkers[sid]['found_mob']['types']) > 0:
            mobtyper_results[sid]['relaxase_type(s)'] = ",".join(biomarkers[sid]['found_mob']['types'])
            mobtyper_results[sid]['relaxase_type_accession(s)'] = ",".join(biomarkers[sid]['found_mob']['acs'])

        if len(biomarkers[sid]['found_mpf']['types']) > 0:
            mpf_type = determine_mpf_type(biomarkers[sid]['found_mpf']['types'])
            mobtyper_results[sid]['mpf_type'] = mpf_type

            mobtyper_results[sid]['mpf_type_accession(s)'] = ",".join(biomarkers[sid]['found_mob']['acs'])

        if len(biomarkers[sid]['found_orit']['types']) > 0:
            mobtyper_results[sid]['orit_type(s)'] = ",".join(biomarkers[sid]['found_orit']['types'])
            mobtyper_results[sid]['orit_accession(s)'] = ",".join(biomarkers[sid]['found_orit']['acs'])




    results = []
    for s in mobtyper_results:
        results.append(mobtyper_results[s])

    writeReport(results, MOB_TYPER_REPORT_HEADER, report_file)



    if not keep_tmp:
        shutil.rmtree(tmp_dir)
    logger.info("Run completed")


# call main function
if __name__ == '__main__':
    main()

