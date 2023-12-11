from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from mob_suite.blast import BlastRunner
from mob_suite.blast import BlastReader
import os, re, time
from subprocess import Popen, PIPE
import shutil, sys, logging, hashlib, random
import pandas as pd
from collections import OrderedDict
from operator import itemgetter
from ete3 import NCBITaxa
from mob_suite.constants import \
    MOB_TYPER_REPORT_HEADER, \
    MGE_INFO_HEADER


def getAssocValues(query_list_values, look_up_key, value_key, data):
    values = []

    for q in query_list_values:
        values.extend(get_data_associated_with_key(look_up_key, q, value_key, data))
    return values


def parseMashScreen(mash_results):
    mash_results = mash_results.decode('utf-8').split("\n")

    hits = {}

    for line in mash_results:
        if len(line) < 4:
            continue
        row = line.strip("\n").split("\t")

        ref_id = str(row[4])
        score = float(row[0])

        if ref_id not in hits:
            hits[ref_id] = {}

        hits[ref_id] = score

    return hits


def parseMash(mash_results):
    mash_results = mash_results.decode('utf-8').split("\n")

    hits = {}

    for line in mash_results:
        row = line.strip("\n").split("\t")
        if len(row) < 4:
            continue
        ref_id = str(row[0])
        query_id = str(row[1])
        score = float(row[2])

        if query_id not in hits:
            hits[query_id] = {}

        hits[query_id][ref_id] = score

    return hits


def filter_invalid_taxids(taxids):
    filtered = []
    for i in taxids:
        if i == 'NaN' or i == 'nan' or i == None or i == '-':
            continue

        try:
            filtered.append(int(i))

        except ValueError:
            continue

    return filtered


def getHeirarchy(taxid, ETE3DBTAXAFILE, database_directory):
    if not isETE3DBTAXAFILEexists(ETE3DBTAXAFILE):
        logging.info("Did not find taxa.sqlite in {}. Initializaing ete3 taxonomy database".format(ETE3DBTAXAFILE))
        initETE3Database(database_directory, ETE3DBTAXAFILE)

    ncbi = NCBITaxa(dbfile=ETE3DBTAXAFILE)
    if not isETE3DBTAXAFILEexists(ETE3DBTAXAFILE):
        logging.error(
            "Tried ete3 init, but still was not able to find taxa.sqlite file for ete3 lib in {}. Aborting".format(
                ETE3DBTAXAFILE))
        return ['-', '-']

    if not isinstance(taxid, int):
        return {'names': [], 'ranks': []}

    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    ranks = []
    for id in lineage:
        ranks.append(ncbi.get_rank(id))

    return {'names': names, 'ranks': names}


def getTaxid(taxon, ETE3DBTAXAFILE, database_directory):
    if not isETE3DBTAXAFILEexists(ETE3DBTAXAFILE):
        logging.info("Did not find taxa.sqlite in {}. Initializaing ete3 taxonomy database".format(ETE3DBTAXAFILE))
        initETE3Database(database_directory, ETE3DBTAXAFILE)

    ncbi = NCBITaxa(dbfile=ETE3DBTAXAFILE)
    if not isETE3DBTAXAFILEexists(ETE3DBTAXAFILE):
        logging.error(
            "Tried ete3 init, but still was not able to find taxa.sqlite file for ete3 lib in {}. Aborting".format(
                ETE3DBTAXAFILE))
        return ['-', '-']

    taxid = ncbi.get_name_translator(taxon)

    return taxid



def NamesToTaxIDs(names, ETE3DBTAXAFILE, database_directory):
    if not isETE3DBTAXAFILEexists(ETE3DBTAXAFILE):
        logging.info("Did not find taxa.sqlite in {}. Initializaing ete3 taxonomy database".format(ETE3DBTAXAFILE))
        initETE3Database(database_directory, ETE3DBTAXAFILE)

    ncbi = NCBITaxa(dbfile=ETE3DBTAXAFILE)

    if not isETE3DBTAXAFILEexists(ETE3DBTAXAFILE):
        logging.error(
            "Tried ete3 init, but still was not able to find taxa.sqlite file for ete3 lib in {}. Aborting".format(
                ETE3DBTAXAFILE))
        return {}

    return ncbi.get_name_translator(names)



def getTaxonConvergence(taxids, ETE3DBTAXAFILE, database_directory):
    if not isETE3DBTAXAFILEexists(ETE3DBTAXAFILE):
        logging.info("Did not find taxa.sqlite in {}. Initializaing ete3 taxonomy database".format(ETE3DBTAXAFILE))
        initETE3Database(database_directory, ETE3DBTAXAFILE)

    ncbi = NCBITaxa(dbfile=ETE3DBTAXAFILE)

    if not isETE3DBTAXAFILEexists(ETE3DBTAXAFILE):
        logging.error(
            "Tried ete3 init, but still was not able to find taxa.sqlite file for ete3 lib in {}. Aborting".format(
                ETE3DBTAXAFILE))
        return ['-', '-']

    tax_ranks = ['genus', 'family', 'order', 'class', 'phylum', 'species']

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
        for i in range(1, len(lineages)):
            common_taxids = list(set(common_taxids) - set(set(common_taxids) - set(lineages[i])))

    if len(common_taxids) == 0:
        return ['-', '-']

    taxa = {}

    for i in range(0, len(common_taxids)):
        rank = ncbi.get_rank([common_taxids[i]])[common_taxids[i]]
        name = ncbi.get_taxid_translator([common_taxids[i]])[common_taxids[i]]

        if not rank in tax_ranks:
            continue

        taxa[rank] = name

    for rank in tax_ranks:
        if rank in taxa:
            return [rank, taxa[rank]]

    phylla = list(set(phylla))
    sorted(phylla)
    if len(phylla) > 1:
        return ('multi-phylla', ",".join(phylla))
    elif len(phylla) == 1:
        return ('phylum', ",".join(phylla))

    return (['-', '-'])


def hostrange(replion_types, relaxase_types, mob_cluster_id, ncbi, lit, ETE3DBTAXAFILE, database_directory):
    host_range_predictions = {
        'observed_host_range_ncbi_name': '',
        'observed_host_range_ncbi_rank': '',
        'reported_host_range_lit_name': '',
        'reported_host_range_lit_rank': '',
        'predicted_host_range_overall_name': '',
        'predicted_host_range_overall_rank': '',
        'associated_pmid(s)': '-'

    }

    # Determine taxids associated with NCBI records
    if mob_cluster_id != '-':
        ncbi_cluster_taxids = getAssocValues([mob_cluster_id], 'primary_cluster_id', 'taxid', ncbi)

    else:
        ncbi_cluster_taxids = []

    if len(replion_types) > 0:
        ncbi_replicon_taxids = getAssocValues(replion_types, 'rep_type(s)', 'taxid', ncbi)
        lit_replicon_taxids = getAssocValues(replion_types, 'rep_type(s)', 'reported_host_range_taxid', lit)
        lit_replicon_taxids = lit_replicon_taxids + getAssocValues(replion_types, 'rep_type(s)', 'host_taxid', lit)
        pmids = list(set(getAssocValues(replion_types, 'rep_type(s)', 'pmid', lit)))
        sorted(pmids)
        pmids = filter_invalid_taxids(list(set(pmids)))
        host_range_predictions['associated_pmid(s)'] = '; '.join(str(x) for x in pmids)

    else:
        ncbi_replicon_taxids = []
        lit_replicon_taxids = []

    if len(relaxase_types) > 0:
        ncbi_relaxase_taxids = getAssocValues(relaxase_types, 'relaxase_type_accession(s)', 'taxid', ncbi)
    else:
        ncbi_relaxase_taxids = []

    if len(ncbi_cluster_taxids) == 0 and len(ncbi_replicon_taxids) == 0 and len(ncbi_relaxase_taxids) == 0:
        host_range_predictions['observed_host_range_ncbi_rank'] = '-'
        host_range_predictions['observed_host_range_ncbi_name'] = '-'
        ncbi_unique_taxids = []
    else:
        ncbi_unique_taxids = filter_invalid_taxids(
            list(set(ncbi_replicon_taxids + ncbi_cluster_taxids + ncbi_relaxase_taxids)))
        host_range_predictions['observed_host_range_ncbi_rank'], host_range_predictions[
            'observed_host_range_ncbi_name'] = getTaxonConvergence(ncbi_unique_taxids, ETE3DBTAXAFILE, database_directory)

    # Determine taxids associated with literature

    lit_unique_taxids = filter_invalid_taxids(list(set(lit_replicon_taxids)))

    host_range_predictions['reported_host_range_lit_rank'], host_range_predictions[
        'reported_host_range_lit_name'] = getTaxonConvergence(lit_unique_taxids, ETE3DBTAXAFILE, database_directory)

    # determine overall host range
    overall_taxids = filter_invalid_taxids(list(set(ncbi_unique_taxids + lit_unique_taxids)))
    host_range_predictions['predicted_host_range_overall_rank'], host_range_predictions[
        'predicted_host_range_overall_name'] = getTaxonConvergence(overall_taxids, ETE3DBTAXAFILE, database_directory)

    # move host-range prediction up to family when it is at genus or species level
    if host_range_predictions['predicted_host_range_overall_rank'] == 'genus' or host_range_predictions[
        'predicted_host_range_overall_rank'] == 'species':
        taxid = getTaxid(host_range_predictions['predicted_host_range_overall_name'], ETE3DBTAXAFILE, database_directory)
        heir = getHeirarchy(taxid, ETE3DBTAXAFILE, database_directory)
        names = heir['names']
        ranks = heir['ranks']

        for i in range(0, len(ranks)):
            if ranks[i] == 'Family':
                host_range_predictions['predicted_host_range_overall_rank'] = 'family'
                host_range_predictions['predicted_host_range_overall_name'] = names[i]

    unique_taxa =  sorted(list(set(host_range_predictions['observed_host_range_ncbi_name'].split(','))))
    taxa =  host_range_predictions['observed_host_range_ncbi_name'].split(",")
    ranks = host_range_predictions['observed_host_range_ncbi_rank'].split(",")
    unique_ranks = []
    for taxon in unique_taxa:
        for i in range(0,len(ranks)):
            if taxon == taxa[i]:
                unique_ranks.append(ranks[i])
                break

    host_range_predictions['observed_host_range_overall_rank'] = ",".join(unique_ranks)
    host_range_predictions['observed_host_range_ncbi_name'] = ",".join(unique_taxa)

    unique_taxa =  sorted(list(set(host_range_predictions['predicted_host_range_overall_name'].split(','))))
    taxa =  host_range_predictions['predicted_host_range_overall_name'].split(",")
    ranks = host_range_predictions['predicted_host_range_overall_rank'].split(",")
    unique_ranks = []
    for taxon in unique_taxa:
        for i in range(0,len(ranks)):
            if taxon == taxa[i]:
                unique_ranks.append(ranks[i])
                break
    host_range_predictions['predicted_host_range_overall_rank'] = ",".join(unique_ranks)
    host_range_predictions['predicted_host_range_overall_name'] = ",".join(unique_taxa)


    return host_range_predictions


def blastn(input_fasta, blastdb, min_ident, min_cov, evalue, min_length, out_dir, blast_results_file, logging,
           seq_filterfile=None, num_threads=1, max_length=400000, min_hsp_cov=1):
    blast_runner = BlastRunner(input_fasta, out_dir)
    blast_runner.run_blast(query_fasta_path=input_fasta, blast_task='megablast', db_path=blastdb,
                           db_type='nucl', min_cov=min_cov, min_ident=min_ident, evalue=evalue,
                           blast_outfile=blast_results_file, logging=logging, num_threads=num_threads, word_size=11,
                           seq_id_file=seq_filterfile)

    if os.path.getsize(blast_results_file) == 0:
        os.remove(blast_results_file)
        return False

    blast_df = BlastReader(blast_results_file, logging).df
    blast_df['length'] = blast_df['length'].astype('int64')
    blast_df['qlen'] = blast_df['qlen'].astype('int64')
    blast_df['qcovs'] = blast_df['qcovs'].astype('float64')
    blast_df['qcovhsp'] = blast_df['qcovhsp'].astype('float64')
    blast_df['evalue'] = blast_df['evalue'].astype('float64')
    blast_df['pident'] = blast_df['pident'].astype('float64')
    blast_df = blast_df.loc[blast_df['length'] >= min_length]
    blast_df = blast_df.loc[blast_df['qlen'] <= max_length]
    blast_df = blast_df.loc[blast_df['qcovs'] >= min_cov]
    blast_df = blast_df.loc[blast_df['qcovhsp'] >= min_hsp_cov]
    blast_df = blast_df.loc[blast_df['evalue'] <= evalue]
    blast_df = blast_df.loc[blast_df['pident'] >= min_ident]

    blast_df = blast_df.reset_index(drop=True)
    blast_df = fixStart(blast_df)
    blast_df.to_csv(blast_results_file, sep='\t', header=True, line_terminator='\n', index=False)

    return True


def tblastn(input_fasta, blastdb, min_ident, min_covs, evalue, out_dir, blast_results_file, logging, num_threads=1,
            min_covhsp=25, seq_id_file=None):
    blast_runner = BlastRunner(input_fasta, out_dir)

    blast_runner.run_tblastn(query_fasta_path=input_fasta, blast_task='megablast', db_path=blastdb,
                             db_type='protein', min_cov=min_covs, min_ident=min_ident, evalue=evalue,
                             blast_outfile=blast_results_file,
                             num_threads=num_threads, seq_id_file=seq_id_file, logging=logging)

    if os.path.getsize(blast_results_file) == 0:
        os.remove(blast_results_file)
        return False

    blast_df = BlastReader(blast_results_file, logging).df

    blast_df = blast_df.loc[blast_df['pident'] >= min_ident]
    blast_df = blast_df.loc[blast_df['qcovs'] >= min_covs]
    blast_df = blast_df.loc[blast_df['qcovhsp'] >= min_covhsp]
    blast_df = blast_df.loc[blast_df['evalue'] <= evalue]
    blast_df = fixStart(blast_df)
    blast_df = blast_df.sort_values(['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])
    blast_df = blast_df.reset_index(drop=True)
    blast_df.to_csv(blast_results_file, sep='\t', header=True, line_terminator='\n', index=False)

    return True


def isETE3DBTAXAFILEexists(ETE3DBTAXAFILE):
    if not os.path.exists(ETE3DBTAXAFILE):
        return False
    else:
        return True


def initETE3Database(database_directory, ETE3DBTAXAFILE):
    lockfilepath = os.path.join(database_directory, ".lock")

    if os.path.exists(lockfilepath) == False:
        open(file=lockfilepath, mode="w").close()
        logging.info("Placed lock file at {}".format(lockfilepath))
    else:
        while os.path.exists(lockfilepath):
            elapsed_time = time.time() - os.path.getmtime(lockfilepath)
            logging.info("Lock file found at {}. Waiting for other processes to finish ete3 database init ...".format(
                lockfilepath))
            logging.info(
                "Elapsed time {} min. Will continue processing after 16 min mark.".format(int(elapsed_time / 60)))
            if elapsed_time >= 1000:
                logging.info(
                    "Elapsed time {} min. Assuming previous process completed all init steps. Continue ...".format(
                        int(elapsed_time / 60)))
                try:  # if previous process failed, no processes are running and > 16 min passed since the lock was created
                    os.remove(lockfilepath)
                except:  # continue if file was removed by other process
                    pass
                break
            time.sleep(60)  # recheck every 1 min if lock file was removed by other process
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



def ETE3_db_status_check(taxid, lockfilepath, ETE3DBTAXAFILE, logging):
    """
    Place a lock file while using ETE3 taxonomy database (taxa.sqlite) to prevent accidental concurrent multiprocess update
    Parameters:
        taxid - the taxonomy id which is 1 by default for database health testing
        lockfilepath - path to the database lock file
        ETE3DBTAXAFILE - path to ETE3 taxa.sqlite file
        logging - logger object for logging messages
    Returns:
        Bool: True/False value with regards to database usage.
              If .lock file is not removed after 10 min, program exits
    """
    max_time = 600
    elapsed_time = 0

    while os.path.exists(lockfilepath) == True and elapsed_time < max_time:
        interval = int(random.randrange(10, 60))
        logging.info(
            "Taxonomy lock file {} exists..waiting {}s for other process to complete".format(lockfilepath, interval))
        time.sleep(interval)
        elapsed_time += interval

    if os.path.exists(lockfilepath) == True:
        # delete lock file older than 1 hour
        modTimesinceEpoc = os.path.getmtime(lockfilepath)
        elapsedTime = (time.time() - modTimesinceEpoc) / 60
        if elapsedTime >= 10:
            os.remove(lockfilepath)
        else:
            logging.error(
                "Taxonomy lock file {} still exists after several attempts. Please delete lock file to continue".format(
                    lockfilepath))
        return False

    else:
        logging.info("Creating Lock file {}".format(lockfilepath))

        #some file systems are read-only which will not support lock file writting
        try:
            open(file=lockfilepath, mode="w").close()
        except Exception as e:
            logging.info(e)
            pass

        logging.info("Testing ETE3 taxonomy db {}".format(ETE3DBTAXAFILE))
        ncbi = NCBITaxa(dbfile=ETE3DBTAXAFILE)

        lineage = ncbi.get_lineage(taxid)

        try:
            os.remove(lockfilepath)
            logging.info("Lock file removed.")
        except Exception as e:
            logging.warning("Lock file is already removed by some other process or read-only file system")
            logging.warning(e)

        if len(lineage) > 0:
            return True
        else:
            return False


'''
INPUT: Read file into pandas df
Returns: list of lines in file with a dict of key value pairs from the header
'''


def read_file_to_dict(file, header, separater="\t"):
    data = pd.read_csv(file, sep=separater, header=0, names=header, encoding="UTF-8")
    records = []
    header_len = len(header)
    for index, row in data.iterrows():
        report = {}
        for i in range(0, header_len):
            name = header[i]
            if name in row:
                report[name] = row[name]
            else:
                report[name] = '-'
        records.append(report)
    return records


'''
    Input: Dictionary with format key :{dict}
    Output: Dictionary indexed by alternate key specified by user
'''


def dict_from_alt_key(dictionary, new_key):
    new_dict = {}
    for id in dictionary:
        values = dictionary[id]
        if new_key in values:
            if not values[new_key] in new_dict:
                new_dict[values[new_key]] = {}
            new_dict[values[new_key]][id] = values
    return new_dict


'''
    Input: Dictionary with format key :{dict}
    Output: Dictionary indexed by alternate key specified by user
'''


def dict_from_alt_key_list(data, new_key):
    new_dict = {}

    for d in data:
        if new_key in d:
            if not d[new_key] in new_dict:
                new_dict[d[new_key]] = {}
            new_dict[d[new_key]] = d

    return new_dict


'''

'''


def get_data_associated_with_key(look_up_key, look_up_value, value_key, dictionary):
    asc_values = []

    for id in dictionary:
        values = dictionary[id]
        if look_up_key in values:
            terms = values[look_up_key].split(',')
            if not look_up_value in terms:
                continue
            if value_key in values and str(values[value_key]) != 'nan':
                asc_values.append(values[value_key])

    return asc_values


'''
    Input: Path to TSV file with MOB_CLUSTER_INFO_HEADER fields as the header lines
    Output: Dictionary of sequence indexed by sequence identifier
'''


def read_sequence_info(file, header):
    if os.path.getsize(file) == 0:
        return dict()
    data = pd.read_csv(file, sep='\t', header=0, names=header)
    sequences = dict()
    for index, row in data.iterrows():
        sample_id = row['sample_id']
        sequences[sample_id] = {}
        for i in range(0, len(header)):
            if header[i] == 'id':
                v = index
            else:
                if str(row[header[i]]) == 'nan':
                    v = '-'
                else:
                    v = row[header[i]]
            sequences[sample_id][header[i]] = v

    return sequences


def check_dependencies(logger):
    external_programs = ['blastn', 'makeblastdb', 'tblastn']
    missing = 0
    for program in external_programs:
        path = shutil.which(program)
        if path is None:
            missing += 1
            logger.error("ERROR: Missing program: {}".format(program, ))
        else:
            logger.info("SUCCESS: Found program {} at {}".format(program, path))
    if missing > 0:
        logger.error("Error, you are missing needed programs for mob-suite, please install them and retry")
        sys.exit(-1)


def calc_md5(seq):
    seq = str(seq).encode()
    md5 = hashlib.md5()
    md5.update(seq)
    return md5.hexdigest()


def fixStart(blast_df):
    for index, row in blast_df.iterrows():

        sstart = int(blast_df.at[index, 'sstart'])
        send = int(blast_df.at[index, 'send'])
        if send < sstart:
            temp = sstart
            blast_df.at[index, 'sstart'] = send
            blast_df.at[index, 'send'] = temp
        qstart = int(blast_df.at[index, 'qstart'])
        qend = int(blast_df.at[index, 'qend'])
        if qend < qstart:
            temp = qstart
            blast_df.at[index, 'qstart'] = qend
            blast_df.at[index, 'qend'] = temp

    blast_df = blast_df.astype({"qstart": 'int64', "qend": 'int64',
                                "sstart": 'int64', "send": 'int64',
                                "evalue": 'float', "bitscore": 'int64',
                                'length': 'int64', 'qlen': 'int64',
                                'slen': 'int64', 'mismatch': 'int64',
                                'pident': 'float', 'qcovs': 'float', 'qcovhsp': 'float'})

    return blast_df


def read_fasta_dict(fasta_file):
    seqs = dict()
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqs[str(record.id)] = str(record.seq)
    handle.close()
    return seqs


'''
    Input: Fasta file and out directory
    Output: writes individual files for each sequence record
'''


def write_individual_fasta(input_fasta, out_dir):
    with open(input_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            out_fasta = os.path.join(out_dir, "{}.fasta".format(record.id))
            with open(out_fasta, "w", encoding="utf-8") as fh:
                fh.write("\n>{}\n{}\n".format(record.id, record.seq))
                fh.close()


def write_fasta_dict(seqs, fasta_file):
    with open(fasta_file, "w") as handle:
        for id in seqs:
            handle.write(">{}\n{}\n".format(id, seqs[id]))
    handle.close()


def verify_init(logger, database_dir):
    mob_init_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'mob_init.py')
    status_file = os.path.join(database_dir, 'status.txt')
    if not os.path.isfile(status_file):
        logger.info('MOB-databases need to be initialized, this will take some time')
        p = Popen([sys.executable, mob_init_path, '-d', database_dir],
                  stdout=PIPE,
                  stderr=PIPE,
                  shell=False)

        stdout, stderr = p.communicate()
        return_code = p.returncode
        logger.info("Return code {}".format(return_code))

        # Verify if no errors were captured during the mob_init script run. Otherwise abort
        if len(re.findall("error", stderr.decode())) != 0 or return_code != 0:
            logger.error("Something went wrong with database download or unpacking"
                         "Check MOB-Suite databases directory and your Internet connection.")
            sys.exit(-1)


def remove_split_hits(blast_df, query_col, contig_id_col, contig_start_col, contig_end_col, bitscore_col):
    blast_df = blast_df.sort_values([query_col, contig_id_col, contig_start_col, contig_end_col, bitscore_col],
                                    ascending=[True, True, True, True, False])
    hits = {}
    for index, row in blast_df.iterrows():
        query_id = row[query_col]
        contig_id = row[contig_id_col]
        contig_start = int(row[contig_start_col])
        contig_end = int(row[contig_end_col])
        length = row['slen']
        score = row['bitscore']
        if contig_start == 1 or contig_end == length:

            if query_id not in hits:
                hits[query_id] = {}

            if contig_id not in hits[query_id]:
                hits[query_id][contig_id] = {'start': {}, 'end': {}}

            if contig_start == 1:
                hits[query_id][contig_id]['start'][index] = {'start': contig_start, 'end': contig_end, 'score': score}
            else:
                hits[query_id][contig_id]['end'][index] = {'start': contig_start, 'end': contig_end, 'score': score}

    filter = []

    for query_id in hits:
        for contig_id in hits[query_id]:
            positions = hits[query_id][contig_id]

            # skip records without both a begining and end hit
            if len(positions['start']) == 0 or len(positions['end']) == 0:
                continue

            start_keys = list(positions['start'].keys())
            end_keys = list(positions['end'].keys())
            max_bitscore = 0
            top_side = ''

            for index in start_keys:
                score = positions['start'][index]['score']
                if score > max_bitscore:
                    max_bitscore = score
                    top_side = 'start'

            for index in end_keys:
                score = positions['end'][index]['score']
                if score > max_bitscore:
                    max_bitscore = score
                    top_side = 'end'

            if top_side == 'start':
                filter += start_keys
            else:
                filter += end_keys

    blast_df = blast_df[~blast_df.index.isin(filter)]
    blast_df.reset_index(drop=True)
    return blast_df


def filter_overlaping_records(blast_df, overlap_threshold, contig_id_col, contig_start_col, contig_end_col,
                              bitscore_col):
    prev_contig_id = ''
    prev_index = -1
    prev_contig_start = -1
    prev_contig_end = -1
    prev_score = -1
    filter_indexes = list()
    exclude_filter = dict()

    for index, row in blast_df.iterrows():

        contig_id = row[contig_id_col]
        contig_start = row[contig_start_col]
        contig_end = row[contig_end_col]
        score = float(row[bitscore_col])


        if prev_contig_id == '':
            prev_index = index
            prev_contig_id = contig_id
            prev_contig_start = contig_start
            prev_contig_end = contig_end
            prev_score = score
            continue

        if contig_id != prev_contig_id:
            prev_index = index
            prev_contig_id = contig_id
            prev_contig_start = contig_start
            prev_contig_end = contig_end
            prev_score = score
            continue
        if prev_index in filter_indexes:
            continue
        if (contig_start >= prev_contig_start and contig_start <= prev_contig_end) or (
                contig_end >= prev_contig_start and contig_end <= prev_contig_end):
            overlap = abs(contig_start - prev_contig_end)


            if overlap > overlap_threshold:
                if prev_score < score:
                    filter_indexes.append(prev_index)
                else:
                    filter_indexes.append(index)

        prev_index = index
        prev_contig_id = contig_id
        prev_contig_start = contig_start
        prev_contig_end = contig_end
        prev_score = score

    for index in exclude_filter:
        filter_indexes.append(index)
    blast_df.drop(filter_indexes, inplace=True)

    return blast_df.reset_index(drop=True)


def recursive_filter_overlap_records(blast_df, overlap_threshold, contig_id_col, contig_start_col, contig_end_col,
                                     bitscore_col):
    size = len(blast_df)
    prev_size = 0

    while size != prev_size:
        blast_df = filter_overlaping_records(blast_df, overlap_threshold, contig_id_col, contig_start_col,
                                             contig_end_col, bitscore_col).sort_values(
            ['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])
        prev_size = size
        size = len(blast_df)
    return blast_df


def fix_fasta_header(in_fasta, out_fasta):
    ids = {}
    fh = open(out_fasta, 'w')
    incr=0
    with open(in_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id = str(record.description)
            status = 'false'
            if 'circular=true' in id or '_circ' in id:
                status = 'true'
            seq = str(record.seq.upper())
            md5 = calc_md5(seq)
            new_id = "{}_{}_circular={}".format(incr,md5,status)
            ids[new_id] = id
            fh.write(">{}\n{}\n".format(new_id,seq))
            incr+=1
    handle.close()
    fh.close()
    return ids


''''
    Accepts fasta file and returns size, number of sequence records and gc %
'''


def calcFastaStatsIndividual(fasta):
    stats = {}
    for record in SeqIO.parse(fasta, "fasta"):
        id = record.id
        seq = record.seq
        genome_size = len(seq)
        gc = gc_fraction(seq)
        md5 = calc_md5(seq)
        stats[id] = {
            'total_length': genome_size,
            'gc': gc,
            'md5': md5
        }

    return stats


''''
    Accepts fasta file and returns size, number of sequence records and gc %
'''


def calcFastaStats(fasta):
    num_seqs = 0;
    seq = ''
    for record in SeqIO.parse(fasta, "fasta"):
        num_seqs += 1
        seq = seq + record.seq
    genome_size = len(seq)
    gc = gc_fraction(seq)
    md5 = calc_md5(seq)

    return {
        'num_seq': num_seqs,
        'size': genome_size,
        'gc_content': gc,
        'md5': md5
    }


def filterFastaByIDs(in_fasta, out_fasta, include_list):
    fh = open(out_fasta, 'w')

    for record in SeqIO.parse(in_fasta, "fasta"):
        if record.id in include_list:
            fh.write(">{}\n{}\n".format(record.id, record.seq))

    fh.close()


def init_console_logger(lvl=2):
    LOG_FORMAT = '%(asctime)s %(name)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'

    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    report_lvl = logging_levels[min(lvl, 3)]

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)
    return logging.getLogger(__name__)


def writeReport(data_list, header, outfile):
    with open(outfile, 'w') as fh:
        fh.write("{}\n".format("\t".join(header)))
        for i in range(0, len(data_list)):
            data = data_list[i]
            row = []
            for h in header:
                if h in data:
                    if data[h] != '':
                        row.append(str(data[h]))
                    else:
                        row.append('-')
                else:
                    row.append("-")
            fh.write("{}\n".format("\t".join(row)))
        fh.close()


def sort_biomarkers(biomarker_dict):
    for id in biomarker_dict:
        acs = biomarker_dict[id]['acs']
        types = biomarker_dict[id]['types']

        if len(acs) == 0 or not isinstance(acs, list):
            continue

        tmp_dict = {}

        for i in range(0, len(acs)):
            tmp_dict[acs[i]] = types[i]

        tmp_dict = OrderedDict(sorted(tmp_dict.items(), key=itemgetter(1), reverse=False))

        biomarker_dict[id]['acs'] = list(tmp_dict.keys())
        biomarker_dict[id]['types'] = list(tmp_dict.values())

    return biomarker_dict


def determine_mpf_type(hits):
    if len(hits) == 0 or isinstance(hits, str):
        return hits
    types = dict()
    for hit in hits:
        if not hit in types:
            types[hit] = 0
        types[hit] += 1

    return max(types, key=lambda i: types[i])


def build_mobtyper_report(plasmid_contig_info, out_dir, outfile, seq_dict, ncbi, lit, ETE3DBTAXAFILE, database_directory):
    mob_typer_results = {}
    for clust_id in plasmid_contig_info:

        cluster_file = open(os.path.join(out_dir, "plasmid_{}.fasta".format(clust_id)), 'w')
        logging.info(
            "Writting plasmid sequences to {}".format(os.path.join(out_dir, "plasmid_{}.fasta".format(clust_id))))
        if clust_id not in mob_typer_results:
            mob_typer_results[clust_id] = {}
        for field in MOB_TYPER_REPORT_HEADER:
            if not field in mob_typer_results[clust_id]:
                mob_typer_results[clust_id][field] = []

        data = plasmid_contig_info[clust_id]

        # Put contig report data into MOB-typer report header
        # aggregating the data by cluster id
        cluster_seq = []
        for contig_id in data:
            if contig_id in seq_dict:
                cluster_file.write(">{}\n{}\n".format(contig_id, seq_dict[contig_id]))
                cluster_seq.append(seq_dict[contig_id])
            for field in MOB_TYPER_REPORT_HEADER:
                if field in data[contig_id]:
                    if isinstance(data[contig_id][field], list) and len(data[contig_id][field]) == 0:
                        continue
                    if data[contig_id][field] == '':
                        continue
                    mob_typer_results[clust_id][field].append(data[contig_id][field])

        # overwrite individual seq stat calculations with the overall
        cluster_seq = sorted(cluster_seq,key=len)
        seq = "".join(cluster_seq)
        mob_typer_results[clust_id]['md5'] = [calc_md5(seq)]
        mob_typer_results[clust_id]['gc'] = [gc_fraction(seq)]
        mob_typer_results[clust_id]['size'] = [len(seq)]
        mob_typer_results[clust_id]['num_contigs'] = len(cluster_seq)

        # Sort MOB-typer biomarker results
        replicon = sort_biomarkers({'rep': {'types': mob_typer_results[clust_id]['rep_type(s)'],
                                            'acs': mob_typer_results[clust_id]['rep_type_accession(s)']}})
        mob_typer_results[clust_id]['rep_type(s)'] = ",".join(replicon['rep']['types'])
        mob_typer_results[clust_id]['rep_type_accession(s)'] = ",".join(replicon['rep']['acs'])

        relaxase = sort_biomarkers({'mob': {'types': mob_typer_results[clust_id]['relaxase_type(s)'],
                                            'acs': mob_typer_results[clust_id]['relaxase_type_accession(s)']}})
        mob_typer_results[clust_id]['relaxase_type(s)'] = ",".join(relaxase['mob']['types'])
        mob_typer_results[clust_id]['relaxase_type_accession(s)'] = ",".join(relaxase['mob']['acs'])

        if len(mob_typer_results[clust_id]['mpf_type']) > 0:
            if isinstance(mob_typer_results[clust_id]['mpf_type'], str):
                mob_typer_results[clust_id]['mpf_type'] = determine_mpf_type(
                    mob_typer_results[clust_id]['mpf_type'].split(','))
            else:
                tmp = []
                for i in range(0, len(mob_typer_results[clust_id]['mpf_type'])):
                    tmp += mob_typer_results[clust_id]['mpf_type'][i].split(',')
                mob_typer_results[clust_id]['mpf_type'] = determine_mpf_type(tmp)
                del (tmp)
        else:
            mob_typer_results[clust_id]['mpf_type'] = ''

        mob_typer_results[clust_id]['mpf_type_accession(s)'] = ",".join(
            mob_typer_results[clust_id]['mpf_type_accession(s)'])
        mob_typer_results[clust_id]['orit_type(s)'] = ",".join(mob_typer_results[clust_id]['orit_type(s)'])
        mob_typer_results[clust_id]['orit_accession(s)'] = ",".join(mob_typer_results[clust_id]['orit_accession(s)'])

        # Assign mobility
        mob_typer_results[clust_id]['predicted_mobility'] = 'non-mobilizable'
        if len(mob_typer_results[clust_id]['relaxase_type(s)']) > 0 and len(
                mob_typer_results[clust_id]['mpf_type']) > 0:
            mob_typer_results[clust_id]['predicted_mobility'] = 'conjugative'
        elif (len(mob_typer_results[clust_id]['relaxase_type(s)']) > 0 or len(
                mob_typer_results[clust_id]['orit_type(s)']) > 0):
            mob_typer_results[clust_id]['predicted_mobility'] = 'mobilizable'

        if isinstance(mob_typer_results[clust_id]['primary_cluster_id'], list) and len(
                mob_typer_results[clust_id]['primary_cluster_id']) > 0:
            mob_typer_results[clust_id]['primary_cluster_id'] = mob_typer_results[clust_id]['primary_cluster_id'][0]
            mob_typer_results[clust_id]['mash_nearest_neighbor'] = mob_typer_results[clust_id]['mash_nearest_neighbor'][
                0]
            mob_typer_results[clust_id]['mash_neighbor_distance'] = \
            mob_typer_results[clust_id]['mash_neighbor_distance'][0]
            mob_typer_results[clust_id]['mash_neighbor_identification'] = \
            mob_typer_results[clust_id]['mash_neighbor_identification'][0]

        if isinstance(mob_typer_results[clust_id]['sample_id'], list) and len(
                mob_typer_results[clust_id]['sample_id']) > 0:
            mob_typer_results[clust_id]['sample_id'] = "{}:{}".format(mob_typer_results[clust_id]['sample_id'][0],
                                                                      mob_typer_results[clust_id]['primary_cluster_id'])

        if isinstance(mob_typer_results[clust_id]['secondary_cluster_id'], list) and len(
                mob_typer_results[clust_id]['secondary_cluster_id']) > 0:
            mob_typer_results[clust_id]['secondary_cluster_id'] = mob_typer_results[clust_id]['secondary_cluster_id'][0]

        if mob_typer_results[clust_id]['num_contigs'] > 1:
            mob_typer_results[clust_id]['circularity_status'] = 'incomplete'

        if len(replicon['rep']['types']) > 0:
            rep_types = replicon['rep']['types']
        else:
            rep_types = []

        if len(relaxase['mob']['acs']) > 0:
            relaxase_types = relaxase['mob']['acs']
        else:
            relaxase_types = []

        if mob_typer_results[clust_id]['primary_cluster_id'] != '' and not 'novel_' in mob_typer_results[clust_id][
            'primary_cluster_id']:
            mob_cluster_id = mob_typer_results[clust_id]['primary_cluster_id']
        else:
            mob_cluster_id = '-'

        host_range = hostrange(rep_types, relaxase_types, mob_cluster_id, ncbi, lit, ETE3DBTAXAFILE, database_directory)

        for field in host_range:
            mob_typer_results[clust_id][field] = host_range[field]

        for element in MOB_TYPER_REPORT_HEADER:
            if element in mob_typer_results[clust_id]:
                if isinstance(mob_typer_results[clust_id][element], list):
                    mob_typer_results[clust_id][element] = ', '.join(
                        str(x) for x in mob_typer_results[clust_id][element])

    results = []
    for clust_id in mob_typer_results:
        results.append(mob_typer_results[clust_id])

    writeReport(results, MOB_TYPER_REPORT_HEADER, outfile)
    return


def add_biomarker_results(biomarker_df, df_column_name_biomarker, df_column_name_seqid, contig_info,
                          contig_info_type_key, contig_info_acs_key, type_col_num=1, type_acs_num=0, delimeter='|'):
    results = {}
    for index, row in biomarker_df.iterrows():

        biomarker = row[df_column_name_biomarker].split(delimeter)

        if type_col_num + 1 > len(biomarker):
            logging.error(
                "Specified column for biomarker type: {} does not exist in biomarker id field: {}".format(type_col_num,
                                                                                                          row[
                                                                                                              df_column_name_biomarker]))
            continue

        if type_acs_num + 1 > len(biomarker):
            logging.error("Specified column for biomarker type: {} does not exist in biomarker id field: {}".format(
                type_acs_num, row[df_column_name_biomarker]))
            continue

        seqid = row[df_column_name_seqid]
        if seqid[len(seqid) - 1] == '|':
            tmp = seqid.split('|')
            #
            # tmp[0] in ['gb','emb','ddj','ref'] and
            if len(tmp) > 1:
                seqid = tmp[1]

        if not seqid in results:
            results[seqid] = {'types': [], 'acs': []}

        results[seqid]['types'].append(biomarker[type_col_num])
        results[seqid]['acs'].append(biomarker[type_acs_num])

    results = sort_biomarkers(results)

    for seqid in results:
        if not seqid in contig_info:
            continue

        if contig_info_type_key in contig_info[seqid]:
            contig_info[seqid][contig_info_type_key] = ','.join(results[seqid]['types'])
        else:
            logging.error(
                "Error: {} contig_info_type_key not found in contig_info, check name and fix".format(
                    contig_info_type_key))

        if contig_info_acs_key in contig_info[seqid]:
            contig_info[seqid][contig_info_acs_key] = ','.join(results[seqid]['acs'])
        else:
            logging.error(
                "Error: {} contig_info_acs_key not found in contig_info, check name and fix".format(
                    contig_info_type_key))

    return contig_info


def identify_biomarkers(contig_info, fixed_fasta, tmp_dir, min_length, logging,
                        replicon_ref, min_rep_ident, min_rep_cov, min_rep_evalue, replicon_blast_results,
                        mob_ref, min_mob_ident, min_mob_cov, min_mob_evalue, mob_blast_results,
                        mpf_ref, min_mpf_ident, min_mpf_cov, min_mpf_evalue, mpf_blast_results,
                        repetitive_mask_file, min_rpp_ident, min_rpp_cov, min_rpp_evalue,
                        plasmid_orit, orit_blast_results, repetitive_blast_results,
                        num_threads=1):
    # blast replicon database
    logging.info("Blasting replicon sequences {} against {}".format(replicon_ref, fixed_fasta))
    blastn(input_fasta=replicon_ref, blastdb=fixed_fasta, min_ident=min_rep_ident, min_cov=min_rep_cov,
           evalue=min_rep_evalue, min_length=80, out_dir=tmp_dir,
           blast_results_file=replicon_blast_results, num_threads=num_threads, logging=logging, min_hsp_cov=30)

    logging.info("Filtering replicon blast results {} ".format(replicon_blast_results))
    rep_blast_df = BlastReader(replicon_blast_results, logging=logging).df
    if len(rep_blast_df) > 0:
        rep_blast_df = rep_blast_df.drop(0)
        rep_blast_df = fixStart(rep_blast_df)
        rep_blast_df = remove_split_hits(rep_blast_df, 'qseqid', 'sseqid', 'sstart', 'send', 'bitscore').sort_values(
            ['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])

        rep_blast_df = recursive_filter_overlap_records(rep_blast_df, 5, 'sseqid', 'sstart', 'send', 'bitscore')

        contig_info = add_biomarker_results(biomarker_df=rep_blast_df, df_column_name_biomarker='qseqid',
                                            df_column_name_seqid='sseqid', contig_info=contig_info,
                                            contig_info_type_key='rep_type(s)',
                                            contig_info_acs_key='rep_type_accession(s)', delimeter='|')

    del (rep_blast_df)

    # blast relaxase database
    logging.info("Blasting relaxase sequences {} against {}".format(mob_ref, fixed_fasta))
    tblastn(input_fasta=mob_ref, blastdb=fixed_fasta, min_ident=min_mob_ident, min_covs=min_mob_cov,
            evalue=min_mob_evalue, out_dir=tmp_dir, logging=logging,
            blast_results_file=mob_blast_results, num_threads=num_threads)

    logging.info("Filtering relaxase blast results {} ".format(mob_blast_results))

    mob_blast_df = BlastReader(mob_blast_results, logging).df
    if len(mob_blast_df) > 0:
        mob_blast_df = fixStart(mob_blast_df.drop(0).sort_values(['sseqid', 'sstart', 'send', 'bitscore'],
                                                                 ascending=[True, True, True, False]))
        mob_blast_df = remove_split_hits(mob_blast_df, 'qseqid', 'sseqid', 'sstart', 'send', 'bitscore').sort_values(
            ['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])
        mob_blast_df = recursive_filter_overlap_records(mob_blast_df, 5, 'sseqid', 'sstart', 'send',
                                                        'bitscore')

        add_biomarker_results(biomarker_df=mob_blast_df, df_column_name_biomarker='qseqid',
                              df_column_name_seqid='sseqid', contig_info=contig_info,
                              contig_info_type_key='relaxase_type(s)', contig_info_acs_key='relaxase_type_accession(s)',
                              delimeter='|')

    del (mob_blast_df)

    # blast mpf database
    logging.info("Blasting MPF sequences {} against {}".format(mpf_ref, fixed_fasta))
    tblastn(input_fasta=mpf_ref, blastdb=fixed_fasta, min_ident=min_mpf_ident, min_covs=min_mpf_cov,
            evalue=min_mpf_evalue, out_dir=tmp_dir,
            blast_results_file=mpf_blast_results, num_threads=num_threads, logging=logging)

    mpf_blast_df = BlastReader(mpf_blast_results, logging).df

    if len(mpf_blast_df) > 0:
        mpf_blast_df = fixStart(mpf_blast_df.drop(0)).sort_values(['sseqid', 'sstart', 'send', 'bitscore'],
                                                                  ascending=[True, True, True, False])
        mpf_blast_df = remove_split_hits(mpf_blast_df, 'qseqid', 'sseqid', 'sstart', 'send', 'bitscore').sort_values(
            ['sseqid', 'sstart', 'send', 'bitscore'], ascending=[True, True, True, False])
        mpf_blast_df = recursive_filter_overlap_records(mpf_blast_df, 5, 'sseqid', 'sstart', 'send',
                                                        'bitscore')

        logging.info("Filtering MPF blast results {} ".format(mpf_blast_results))

        add_biomarker_results(biomarker_df=mpf_blast_df, df_column_name_biomarker='qseqid',
                              df_column_name_seqid='sseqid', contig_info=contig_info,
                              contig_info_type_key='mpf_type', contig_info_acs_key='mpf_type_accession(s)',
                              delimeter='|')
    del (mpf_blast_results)

    # blast orit database
    logging.info("Blasting orit sequences {} against {}".format(plasmid_orit, fixed_fasta))
    blastn(input_fasta=plasmid_orit, blastdb=fixed_fasta, min_ident=min_rep_ident, min_cov=min_rep_cov,
           evalue=min_rep_evalue, min_length=80, out_dir=tmp_dir,
           blast_results_file=orit_blast_results, num_threads=num_threads, logging=logging)

    logging.info("Filtering orit blast results {} ".format(orit_blast_results))

    orit_blast_df = BlastReader(orit_blast_results, logging).df
    if len(orit_blast_df) > 0:
        orit_blast_df = recursive_filter_overlap_records(fixStart(
            orit_blast_df.drop(0).sort_values(['sseqid', 'sstart', 'send', 'bitscore'],
                                              ascending=[True, True, True, False])), 5, 'sseqid', 'sstart', 'send',
                                                         'bitscore')
        add_biomarker_results(biomarker_df=orit_blast_df, df_column_name_biomarker='qseqid',
                              df_column_name_seqid='sseqid', contig_info=contig_info,
                              contig_info_type_key='orit_type(s)', contig_info_acs_key='orit_accession(s)',
                              delimeter='|')

    del (orit_blast_df)

    # blast repetitive database
    if repetitive_mask_file is not None:
        logging.info("Blasting contigs against repetitive sequences db: {}".format(repetitive_mask_file))
        blastn(input_fasta=fixed_fasta, blastdb=repetitive_mask_file, min_ident=min_rpp_ident, min_cov=min_rpp_cov,
               evalue=min_rpp_evalue, min_length=min_length, out_dir=tmp_dir,
               blast_results_file=repetitive_blast_results, num_threads=num_threads, logging=logging)
        logging.info("Filtering repetitive blast results {} ".format(repetitive_blast_results))

        repetitive_blast_df = BlastReader(repetitive_blast_results, logging).df
        if len(repetitive_blast_df) > 0:
            repetitive_blast_df = recursive_filter_overlap_records(fixStart(
                repetitive_blast_df.drop(0).sort_values(['sseqid', 'sstart', 'send', 'bitscore'],
                                                        ascending=[True, True, True, False])), 5, 'qseqid', 'qstart',
                                                                   'qend',
                                                                   'bitscore')

            repetitive_list = repetitive_blast_df['qseqid'].tolist()

            # add filtering flag to contigs which are primarially a repetitive element
            for contig_id in repetitive_list:
                if contig_id in contig_info:
                    logging.info('Filtering contig: {} due to repetitive sequence'.format(contig_id))
                    contig_info[contig_id]['filtering_reason'] = 'repetitve element'
                else:
                    logging.error('Contig: {} not found in contig_df this is likely an error'.format(contig_id))

            add_biomarker_results(biomarker_df=repetitive_blast_df, df_column_name_biomarker='sseqid',
                                  df_column_name_seqid='qseqid', contig_info=contig_info,
                                  contig_info_type_key='repetitive_dna_type', contig_info_acs_key='repetitive_dna_id',
                                  delimeter='|', type_col_num=2, type_acs_num=1)

        del (repetitive_blast_df)
    return contig_info

def blast_mge(contig_fasta, mge_fasta,tmp_dir, min_length, logging, min_rpp_ident, min_rpp_cov, min_rpp_evalue,num_threads=1):
    # blast repetitive database
    blast_results = os.path.join(tmp_dir,"mge.blast.results.txt")
    blastn(input_fasta=mge_fasta, blastdb=contig_fasta, min_ident=min_rpp_ident, min_cov=min_rpp_cov,
           evalue=min_rpp_evalue, min_length=min_length, out_dir=tmp_dir,
           blast_results_file=blast_results, num_threads=num_threads, logging=logging)

    #return empty dataframe if no blast results generated
    if not os.path.isfile(blast_results) or os.path.getsize(blast_results) == 0:
        return {}

    df = BlastReader(blast_results, logging).df
    reduced_df = recursive_filter_overlap_records(fixStart(
        df.drop(0).sort_values(['sseqid', 'sstart', 'send', 'bitscore'],
                                                ascending=[True, True, True, False])), 5, 'sseqid', 'sstart',
        'send',
        'bitscore')

    results = {}
    for index,row in reduced_df.iterrows():
        contig_id = row['sseqid']
        if not contig_id in results:
            results[contig_id] = []
        results[contig_id].append(row.to_dict())

    return results

def writeMGEresults(contig_membership,mge_results,outfile):
    out_string = ["\t".join(MGE_INFO_HEADER)]
    if len(mge_results) == 0:
        return
    for contig_id in contig_membership['chromosome']:
        if not contig_id in mge_results:
            continue

        for i in range(0, len(mge_results[contig_id])):
            row = {}
            for field in MGE_INFO_HEADER:
                row[field] = ''
                if field in mge_results[contig_id]:
                    row[field] = mge_results[contig_id][field]

            id = mge_results[contig_id][i]['qseqid'].split('|')
            row['mge_id'] = id[0]
            row['mge_acs'] = id[1]
            row['mge_type'] = id[2]
            row['mge_subtype'] = id[3]
            row['mge_length'] = mge_results[contig_id][i]['qlen']
            row['mge_start'] = mge_results[contig_id][i]['qstart']
            row['mge_end'] = mge_results[contig_id][i]['qend']
            row['contig_start'] = mge_results[contig_id][i]['sstart']
            row['contig_end'] = mge_results[contig_id][i]['send']

            for field in contig_membership['chromosome'][contig_id]:
                if field in row:
                    row[field] = contig_membership['chromosome'][contig_id][field]

            out_string.append("\t".join([str(x) for x in list(row.values())]))

    for mobcluster in contig_membership['plasmid']:
        for contig_id in contig_membership['plasmid'][mobcluster]:
            if not contig_id in mge_results:
                continue

            for i in range(0, len(mge_results[contig_id])):
                row = {}
                for field in MGE_INFO_HEADER:
                    row[field] = ''
                    if field in mge_results[contig_id][i]:
                        row[field] = mge_results[contig_id][i][field]
                id = mge_results[contig_id][i]['qseqid'].split('|')
                row['mge_id'] = id[0]
                row['mge_acs'] = id[1]
                row['mge_type'] = id[2]
                row['mge_subtype'] = id[3]
                row['mge_length'] =  mge_results[contig_id][i]['qlen']
                row['mge_start'] = mge_results[contig_id][i]['qstart']
                row['mge_end'] = mge_results[contig_id][i]['qend']
                row['contig_start'] = mge_results[contig_id][i]['sstart']
                row['contig_end'] = mge_results[contig_id][i]['send']

                for field in contig_membership['plasmid'][mobcluster][contig_id]:
                    if field in row:
                        row[field] = contig_membership['plasmid'][mobcluster][contig_id][field]

                out_string.append("\t".join([str(x) for x in list(row.values())]))


    fh = open(outfile,'w')
    fh.write("\n".join(out_string))
    fh.close()


def create_biomarker_dataframe(parameters,id_mapping,logging):
    logging.info("Creating plasmid biomarkers report ...")
    data_frames = []
    for label in parameters:
        file = parameters[label]['file']  
        if not os.path.isfile(file):
            continue
        blast_df = pd.read_csv(file, header=0, sep="\t")
     
        if len(blast_df) == 0:
            continue
        blast_df['sseqid'] = blast_df['sseqid'].replace(id_mapping)
        blast_df['length'].astype('int64')
        blast_df['qcovs'].astype('float64')
        blast_df['qcovhsp'].astype('float64')
        blast_df['evalue'].astype('float64')
        blast_df['pident'].astype('float64')
        blast_df = blast_df.loc[blast_df['length'] >= parameters[label]['min_length']]
        blast_df = blast_df.loc[blast_df['qcovs'] >= parameters[label]['min_cov']]
        blast_df = blast_df.loc[blast_df['qcovhsp'] >= parameters[label]['min_hsp_cov']]
        blast_df = blast_df.loc[blast_df['evalue'] <= parameters[label]['evalue']]
        blast_df = blast_df.loc[blast_df['pident'] >= parameters[label]['min_ident']]
        blast_df = blast_df.reset_index(drop=True)

        blast_df = fixStart(blast_df)
        if label == 'replicon' or label == 'relaxase' or label == 'mate-pair-formation':
            blast_df = remove_split_hits(blast_df,'qseqid', 'sseqid', 'sstart', 'send', 'bitscore')
        blast_df.sort_values(['sseqid', 'sstart', 'send', 'bitscore'],ascending=[True, True, True, False], inplace=True)
    
        blast_df = recursive_filter_overlap_records(blast_df, 5, 'sseqid', 'sstart', 'send', 'bitscore')
        blast_df = blast_df.reset_index(drop=True)
        blast_df['biomarker'] = label
        data_frames.append(blast_df)
    if len(data_frames) > 0:    
        return pd.concat(data_frames)
    else:
        return pd.DataFrame()
