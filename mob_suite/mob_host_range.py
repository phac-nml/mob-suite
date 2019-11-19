#!/usr/bin/env python
import pandas, re, logging, os,time
pandas.set_option('display.max_rows', 500)
pandas.set_option('display.max_columns', 500)
pandas.set_option('display.width', 1000)
import subprocess
from argparse import ArgumentParser
from collections import Counter,OrderedDict
from ete3 import NCBITaxa
from mob_suite.version import __version__
from mob_suite.utils import default_database_dir,init_console_logger

database_directory = os.path.abspath(default_database_dir)
ETE3DBTAXAFILE = os.path.abspath(database_directory + "/taxa.sqlite")

root_log_level = logging.getLogger().getEffectiveLevel()
logger = init_console_logger(root_log_level)

#pandas.options.display.float_format = '{:.1E}'.format #render scientific notation

#default init arguments
#args=ArgumentParser()
# args.multi_match = True
# args.exact_match = None
# args.loose_match = None
# args.render_tree = True
# args.write_newick = True
# args.debug = False
# args.relaxase_accession = None
# args.outdir = "./"
# args.inputseq = False

#OUTDIR = os.getcwd()+"/"

# def createLogger():
#     log = logging.getLogger('MOB-Suite')
#     formatter = logging.Formatter(
#         '%(asctime)s %(name)-12s %(message)s')
#     log.setLevel(logging.DEBUG)
#
#     console = logging.StreamHandler()
#     console.setFormatter(formatter)
#     console.setLevel(logging.DEBUG)
#     log.addHandler(console)
#     return log

#LOG = createLogger()

def isETE3DBTAXAFILEexists():
    if not os.path.exists(ETE3DBTAXAFILE):
        return False
    else:
        return True
def initETE3Database():
    lockfilepath = os.path.join(database_directory, ".lock")

    if os.path.exists(lockfilepath) == False:
        open(file=lockfilepath, mode="w").close()
        logger.info("Placed lock file at {}".format(lockfilepath))
    else:
        while os.path.exists(lockfilepath):
            elapsed_time = time.time() - os.path.getmtime(lockfilepath)
            logger.info("Lock file found at {}. Waiting for other processes to finish ete3 database init ...".format(lockfilepath))
            logger.info("Elapsed time {} min. Will continue processing after 16 min mark.".format(int(elapsed_time/60)))
            if elapsed_time >= 1000:
                logger.info("Elapsed time {} min. Assuming previous process completed all init steps. Continue ...".format(int(elapsed_time/60)))
                try: #if previous process failed, no processes are running and > 16 min passed since the lock was created
                    os.remove(lockfilepath)
                except: #continue if file was removed by other process
                    pass
                break
            time.sleep(60) #recheck every 1 min if lock file was removed by other process
        logger.info("Lock file no longer exists. Assuming init process completed successfully")

    ncbi = NCBITaxa()
    ncbi.dbfile = ETE3DBTAXAFILE
    ncbi.update_taxonomy_database()

    try:
        os.remove(lockfilepath)
        logger.info("Lock file removed.")
    except:
        logger.warning("Lock file is already removed by some other process.")
        pass

    try:
        os.remove(os.path.join(os.getcwd(), "taxdump.tar.gz"))
        logger.info("Removed residual taxdump.tar.gz as ete3 is not doing proper cleaning job.")
    except:
        pass
    logger.info("ETE3 database init completed successfully.")

def loadliteratureplasmidDB():
    literatureplasmidDB = pandas.read_csv(os.path.dirname(os.path.abspath(__file__))+"/databases/host_range_literature_plasmidDB.csv",
                            sep=",",encoding = "ISO-8859-1",dtype={"PMID":str, "TransferRate":float, "Year":str,"Size":float})


    return literatureplasmidDB

def findHitsInLiteratureDBbyReplicon(replicon_names,plasmid_lit_db):
    """
    Get me the indices/hits per each replicon name and adjust if the replicon name returns no hit, then use inc family search
    :param replicon_names: names of query replicons (e.g. IncFII)
    :param plasmid_lit_db: pandas dataframe with all literature information on plasmids
    :return: repliconsearchdict: dictionary replicon name / hit indices (i.e. rows) (e.g. {'IncF': [0, 1, 2, 3, 4, 5, 6, 7, 97, 123, 135]})
    """

    # Match by Inc families if exact match returns no hits (e.g. IncF instead of IncFII)
    repliconsearchdict = {}
    for replicon_name in replicon_names:
        db_hit_indices=[i for i in range(0, plasmid_lit_db.shape[0]) if plasmid_lit_db.iloc[i, :]["Replicon"] == replicon_name]
        if len(db_hit_indices) != 0:
            repliconsearchdict[replicon_name] = db_hit_indices
        elif len(re.findall("Inc[A-Z]{1}|ColRNA", replicon_name)) != 0:
            db_hit_indices = [i for i in range(0, plasmid_lit_db.shape[0]) if len(re.findall(replicon_name, plasmid_lit_db.iloc[i, :]["Replicon"])) != 0]
            repliconsearchdict[replicon_name] = db_hit_indices

        #re.findall("Inc[A-Z]{1}|Col", replicon_name)
        # search literature database by replicon family (e.g. IncF) if no results are returned
        #if len(idx) == 0: # if exact match failed, search by replicon family
        #    replicon_name = re.findall("Inc[A-Z]{1}|Col", replicon_name)[0]


    return repliconsearchdict

def getLiteratureBasedHostRange(replicon_names,plasmid_lit_db,input_seq=""):
    '''
    # Read literature-based database file
    # Make HR prediction based on the literature evidence on replicon and mob-suite cluster
    # If replicon is not know or mash dist to the nearest reference sequence is > 0.05 ... do not do any prediction
    # Indicate if the plasmid is broad or narrow range class
    :param replicon_names: query replicon names to query literature database
    :param plasmid_lit_db: pandas dataframe representing a database of well curated plasmid mentioned in the literature
    :param input_seq: path to the input plasmid sequence
    :return: report_df: pandas dataframe with all relevant information on the input query including literature reported host range and typical plasmids and PMIDs
             lit_taxids: list of taxonomy ids
    '''


    report_df=pandas.DataFrame()
    repliconsearchdict = findHitsInLiteratureDBbyReplicon(replicon_names,plasmid_lit_db)
    lit_taxids_list=[]
    #report_dict=dict.fromkeys(['apple','ball'],"-")


    for replicon_name in repliconsearchdict.keys():
        idx = repliconsearchdict[replicon_name] #find hits in database based on replicon query
        literature_knowledge = plasmid_lit_db.iloc[idx,:].copy() #get database subset for further analyses

        if literature_knowledge.shape[0] == 0:
            continue
            #raise Exception("Could not extract any information from the literature database matching the search")

        #remove any entries that are map to uncultured bacteria (avoids inflation of the host range)
        select_vector = [True if len(re.findall("uncultured", line)) == 0 else False for line in literature_knowledge["IsolationSpecies"]]
        lit_taxids=list(set(literature_knowledge.loc[select_vector,"IsolationTaxid"]))
        litrank,rankname = getHostRangeRankCovergence(lit_taxids)


        host_range_literature_rank_claim = literature_knowledge["HostRangeRankClaim"].values
        host_range_literature_name_claim = literature_knowledge["HostRangeClaim"].values
        dict_hr_lit_reported = dict(zip(host_range_literature_rank_claim, host_range_literature_name_claim))

        lit_rank_claim_frequency_dict = sorted(Counter(host_range_literature_rank_claim).items(), key=lambda item: item[1], reverse=True)

        host_range_literature_rank_claim = "-"
        host_range_literature_name_claim = "-"
        for key, value in lit_rank_claim_frequency_dict: #account for NaN values and pick the most frequent literature host-range prediction
            if isinstance(key, str):
                host_range_literature_rank_claim = key
                host_range_literature_name_claim = dict_hr_lit_reported[key]
                break

        hostrangeclassdict = Counter(literature_knowledge["HostRangeClass"])
        HostRangeClass = [key for key in hostrangeclassdict.keys()  if  hostrangeclassdict[key] == max(hostrangeclassdict.values())][0]

        LiteratureMinTransferRateRange = "NA"
        LiteratureMaxTransferRateRange = "NA"
        LiteratureMeanTransferRateRange = "NA"

        temp = []
        for r in literature_knowledge["TransferRate"]:
            if isinstance(r,int) or isinstance(r,float):
                temp.append(r)
        literature_knowledge["TransferRate"] = temp

        if any(literature_knowledge["TransferRate"] > 0):
            LiteratureMinTransferRateRange = min([i for i in literature_knowledge["TransferRate"] if i >= 0])
            LiteratureMeanTransferRateRange = mean([i for i in literature_knowledge["TransferRate"] if i >= 0])
            LiteratureMaxTransferRateRange = max([i for i in literature_knowledge["TransferRate"] if i >= 0])
        elif any(literature_knowledge["TransferRate"] == 0):
            LiteratureMinTransferRateRange = "mobilizable/non-conjugative"
            LiteratureMeanTransferRateRange = "mobilizable/non-conjugative"
            LiteratureMaxTransferRateRange = "mobilizable/non-conjugative"



        report_dict = OrderedDict({  "LiteratureQueryReplicon": ",".join(replicon_names),
                         "LiteratureSearchReplicon": replicon_name,
                         "LiteratureFoundPlasmidsNames": ",".join(set(literature_knowledge["Plasmid_Name"])),
                         "LiteratureFoundPlasmidsNumber": len(set(literature_knowledge["Plasmid_Name"])),
                         "LiteratureReportedHostRangePlasmidClass": HostRangeClass,
                         "LiteratureReportedHostPlasmidSpecies": ",".join(set(literature_knowledge["IsolationSpecies"])),
                         "LiteratureReportedPlasmidHostSpeciesNumber": len(set(literature_knowledge["IsolationSpecies"])),
                         "LiteraturePredictedHostRangeTreeRank": litrank,
                         "LiteraturePredictedHostRangeTreeRankSciName": rankname,
                         "LiteratureReportedHostRangeRankInPubs": host_range_literature_rank_claim,
                         "LiteratureReportedHostRangeNameInPubs": host_range_literature_name_claim,
                         "LiteratureMinTransferRateRange": LiteratureMinTransferRateRange,
                         "LiteratureMaxTransferRateRange": LiteratureMaxTransferRateRange,
                         "LiteratureMeanTransferRateRange": LiteratureMeanTransferRateRange,
                         "LiteraturePMIDs": ";".join(sorted(set([str(int(i)) for i in literature_knowledge.loc[:, "PMID"] if pandas.isna(i) == False]))),
                         "LiteraturePublicationsNumber": len(set(literature_knowledge.loc[literature_knowledge["PMID"].isna() == False, "PMID"]))})

        #if input plasmid sequence is provided, do additional closest match based on the sequence similarity and append to the general report

        if input_seq.strip():  #if plasmid sequence is available
            literature_accession_top_hit, literature_mash_dist_top_hit = getClosestLiteratureRefPlasmid(input_seq)
            #print("L",literature_accession_top_hit, literature_mash_dist_top_hit,literature_accession_top_hit.rsplit("."))
            #print(plasmid_lit_db["NCBI_Accession"])
            #print()
            #by accession search
            #idx = [i for i in range(0, len(literature_knowledge)) if len(re.findall(literature_accession_top_hit.rsplit(".")[0],plasmid_lit_db["NCBI_Accession"].values[i])) != 0]
            #remove version "NC_002638.1" number

            literature_accession_top_hit = literature_accession_top_hit.rsplit(".")[0]

            literature_closest_seq_hit_df = pandas.DataFrame()
            for i in range(0, plasmid_lit_db.shape[0]):
                if re.findall(literature_accession_top_hit, str(plasmid_lit_db.iloc[i]["NCBI_Accession"])):
                    #print(literature_accession_top_hit, plasmid_lit_db.iloc[i]["NCBI_Accession"])
                    literature_closest_seq_hit_df = pandas.concat([literature_closest_seq_hit_df, plasmid_lit_db.iloc[[i]]])

            if literature_closest_seq_hit_df.empty:
                raise Exception("Literature top hit search failed! Check mash top hit return from the getClosestLiteratureRefPlasmid()")


            if literature_closest_seq_hit_df.shape[0] != 1:
                logger.warning("Literature top hit dataframe returned more than a single hit ... Expecting a single top hit.Some entries have multiple hits.")

            #additional fields
            report_dict.update(
                                OrderedDict({"LiteratureClosestRefrencePlasmidAcc": literature_accession_top_hit,
                                "LiteratureClosestReferencePlasmidName":literature_closest_seq_hit_df.loc[:,"Plasmid_Name"].values[0],
                                "LiteratureClosestReferencePlasmidSize": int(literature_closest_seq_hit_df.loc[:,"Size"].values[0]),
                                "LiteratureClosestReferenceMashDistance": literature_mash_dist_top_hit,
                                "LiteratureClosestReferenceDonorStrain": literature_closest_seq_hit_df.loc[:,"Donor"].values[0],
                                "LiteratureClosestReferenceRecipientStrain": literature_closest_seq_hit_df.loc[:, "Recipient"].values[0],
                                "LiteratureClosestReferenceTransferRate": literature_closest_seq_hit_df.loc[:, "TransferRate"].values[0],
                                "LiteratureClosestReferenceConjugationTemperature": literature_closest_seq_hit_df.loc[:, "ConjugationTemperature"].values[0]}))


        report_df = pandas.concat([report_df, pandas.DataFrame.from_dict([report_dict])],sort=False) #append results from different replicons (if multiple are present)
        lit_taxids_list = lit_taxids_list + lit_taxids

    report_df.fillna("-", inplace=True)

    return report_df, lit_taxids_list
    #report_table.to_csv(args.outputprefix+'_literature_report.txt',sep="\t",
    #                    float_format='%.1E', index=False, na_rep="NA",mode="w")

def collapseLiteratureReport(df):
    """
    Aggregate data for easier reporting averaging and maximizing fields
    :param: df: dataframe returned from literature report with the following columns (may be more than 1 row if multi-replicons are present)
    #publications	LiteratureClosestRefrencePlasmid	LiteratureClosestPlasmidName	ClosestLiteraturePlasmidSize	LiteratureClosestMashDistance
    :return: single row dataframe
    """

    collapsedlitdf=pandas.DataFrame.from_records([],columns=df.columns,index=[0])
    collapsedlitdf.iloc[0]=df.iloc[0]

    conversiondict={"NarrowHostRange":1,"WideHostRange":2,"BroadHostRange":3}

    #no rank;no rank;superkingdom;phylum;class;order;family;genus;species group;species
    rankconversiondict={"species":1,"genus":2,"family":3,"order":4,"class":5,"phylum":6,"superkingdom":7}
    numrank2nameconversiondict = {1:"species", 2:"genus",3:"family",4:"order",5:"class", 6:"phylum",7:"superkingdom"}

    collapsedlitdf.loc[0,"LiteratureQueryReplicon"] = ",".join(list(set(df["LiteratureQueryReplicon"].values)))
    collapsedlitdf.loc[0, "LiteratureSearchReplicon"] = ",".join(df.loc[:,"LiteratureSearchReplicon"].values)
    collapsedlitdf.loc[0, "LiteratureFoundPlasmidsNames"] = ",".join(df.loc[:,"LiteratureFoundPlasmidsNames"].values)
    collapsedlitdf.loc[0,"LiteratureFoundPlasmidsNumber"]=sum(df.loc[:, "LiteratureFoundPlasmidsNumber"])
    hrclassnum = max([conversiondict[k] for k in set(df.loc[:, "LiteratureReportedHostRangePlasmidClass"].values)])
    collapsedlitdf.loc[0, "LiteratureReportedHostRangePlasmidClass"] = [k for k in conversiondict if conversiondict[k] == hrclassnum][0]

    species_list=[]
    for value in df.loc[:, "LiteratureReportedHostPlasmidSpecies"].values:
        species_list = species_list+value.split(",")

    collapsedlitdf.loc[0,"LiteratureReportedHostPlasmidSpecies"]=",".join(species_list)
    collapsedlitdf.loc[0, "LiteratureReportedPlasmidHostSpeciesNumber"] = len(species_list)
    hrtreeranknum = max([rankconversiondict[k] for k in df.loc[:, "LiteraturePredictedHostRangeTreeRank"].values])
    collapsedlitdf.loc[0,"LiteraturePredictedHostRangeTreeRank"] = [k for k in rankconversiondict if rankconversiondict[k] == hrtreeranknum][0]
    idx=df.loc[:,"LiteraturePredictedHostRangeTreeRank"] == collapsedlitdf.loc[0,"LiteraturePredictedHostRangeTreeRank"]
    collapsedlitdf.loc[0,"LiteraturePredictedHostRangeTreeRankSciName"] = df[idx]["LiteraturePredictedHostRangeTreeRankSciName"].values[0]


    #print(df.loc[:,"LiteratureReportedHostPlasmidSpecies"].values.tolist())
    #print(collapsedlitdf.loc[0,"LiteratureReportedHostPlasmidSpecies"]);exit(1)
    #idx=[]
    #print(df);exit() #KeyError: 'LiteratureReportedHostRangeInPubs'
    #for i in range(0,df.shape[0]):
    #    if isinstance(df.iloc[i]["LiteratureReportedHostRangeInPubs"],int):
    #        idx.append(i)


    if df.shape[0] > 1:
        collapsedlitdf.loc[0,"LiteratureReportedHostRangeRankInPubs"] = numrank2nameconversiondict[max([rankconversiondict[k] for k in df["LiteratureReportedHostRangeRankInPubs"].values if k != "-"])]
        collapsedlitdf.loc[0,"LiteratureReportedHostRangeNameInPubs"] = dict(zip(df["LiteratureReportedHostRangeRankInPubs"], df["LiteratureReportedHostRangeNameInPubs"]))[collapsedlitdf.loc[0,"LiteratureReportedHostRangeRankInPubs"]]
    else:
        collapsedlitdf.loc[0, "LiteratureReportedHostRangeRankInPubs"] = "-"
        collapsedlitdf.loc[0,"LiteratureReportedHostRangeNameInPubs" ] = "-"


    for field in ["LiteratureMinTransferRateRange","LiteratureMaxTransferRateRange", "LiteratureMeanTransferRateRange"]:
        filteredvals=[i for i in df[field].values if isinstance(i,int) or isinstance(i,float)]
        collapsedlitdf.loc[0,field] = mean(filteredvals)

    collapsedlitdf.loc[0,"LiteraturePMIDs"] = ";".join(df["LiteraturePMIDs"])
    collapsedlitdf.loc[0, "LiteraturePublicationsNumber"] = sum(df["LiteraturePublicationsNumber"])

    return(collapsedlitdf)

def getClosestLiteratureRefPlasmid(input_fasta):
    """
    In cases when plasmid sequence is available run additional sequence-based query of the plasmid reference to
    find closest representative plasmid at < 0.05 mash distance. Since plasmids are genomics puzzles of fragements,
    k-mer based approach of distance estimataion might be of benefit in contrast to more deterministic BLAST-based approaches
    :param query single contig sequence in FASTA format
    :return: closest literature-curated plasmid accession number from the literature database and mash distance between query and reference
    """
    reference_db=os.path.dirname(os.path.abspath(__file__))+"/databases/literature_mined_plasmid_seq_db.fasta.msh"
    #input_fasta="/Users/kirill/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/tests/TestData/IncF/ET11_Ecoli_plasmid_529.fasta"
    p = subprocess.Popen(["mash dist -p 1 "+reference_db+" "+input_fasta],stdout=subprocess.PIPE, shell=True)
    stdout = p.communicate()[0].decode("utf-8")

    records = list()
    for line in stdout.split('\n'):
        items=line.rsplit("\t")
        if len(items) != 1:
            records.append((items))

    mash_dist_df = pandas.DataFrame.from_records(records,columns=["Accession","QueryID","mash_dist","p-value","shared_hashes"])
    top_hit_df = mash_dist_df[mash_dist_df.loc[:, "mash_dist"] == min(mash_dist_df.loc[:, "mash_dist"])]
    return top_hit_df.loc[:,"Accession"].values[0], top_hit_df.loc[:, "mash_dist"].values[0]



def mean(numbers):
    if not all([isinstance(n,float) or isinstance(n,int) for n in numbers]):
        mean = "-"
    else:
        mean = float(sum(numbers)) / max(len(numbers), 1)
    return mean

def getHostRangeRankCovergence(taxids):
    """
    :param taxids - list of taxonomy ids allowing to build a tree
    :return rank - taxonomic rank at which phylogenetic tree branches converge
    :return sci_name - taxonomic rank name at which phylogenetic tree branches converge
    Note if no rank is returned due to small diversity amongst tree terminal nodes due to presence of single taxid or low literatuer coverage
    then traverse the tree and report rank at genus level. Of course this literature based host range prediction might be understatement.
    """

    if not isETE3DBTAXAFILEexists():
        logger.info("Did not find taxa.sqlite in {}. Initializaing ete3 taxonomy database".format(ETE3DBTAXAFILE))
        initETE3Database()

    ncbi = NCBITaxa(dbfile=ETE3DBTAXAFILE)
    if not isETE3DBTAXAFILEexists():
        logger.error("Tried ete3 init, but still was not able to find taxa.sqlite file for ete3 lib in {}. Aborting".format(ETE3DBTAXAFILE))
        exit(-1)

    #taxids=[562,573,1288825,439842]
    tree = ncbi.get_topology(taxids)
    tree.annotate_ncbi_taxa(dbfile=ETE3DBTAXAFILE, taxid_attr='name')


    tree_rank=tree.rank
    tree_sci_name=tree.sci_name
    if tree_rank == "no rank":
        ranks_list=list()
        for node in  tree.traverse():
            ranked_taxids = ncbi.get_lineage(taxid=node.taxid)
            ranks_list =  ranks_list+ranked_taxids
        ranks_list_counts_dict = Counter(ranks_list)
        maxchildnodes = max(ranks_list_counts_dict.values())
        for k in ranks_list_counts_dict.keys():
            rank_temp = ncbi.get_rank([k])
            if ranks_list_counts_dict[k] == maxchildnodes and rank_temp[k] == "genus":
                tree_rank = rank_temp[k]
                tree_sci_name = ncbi.get_taxid_translator([k])[k]
    return(tree_rank,tree_sci_name)


#the main function to process
def getRefSeqHostRange(replicon_name_list,  mob_cluster_id_list, relaxase_name_acc_list, relaxase_name_list, matchtype, hr_obs_data):
    """
    Get NCBI RefSeq host range based on the search parameters either/or replicon, mob_cluster_id, relaxase family, relaxase accession
    :param replicon_name_list: a list of replicons to process ['IncP','IncF']
    :param mob_cluster_id_list:  mob-cluster cluster id
    :param relaxase_name_acc_list:  a list of MOB accession ids ['NC_017627_00068' 'NC_011416_00039']
    :param relaxase_name_list: relaxase family names list ['MOBF' 'MOBP']
    :param matchtype: how we will match fields in the database (exactly = entire field, multi = query might represent part of the field)
    :param hr_obs_data: host range observed data
    :return: convergance_rank: effectively the host range approximated by the convergance rank on the phylogenetic tree
             converged_taxonomy_name: the host range name given by the convergence rank
             unique_ref_selected_taxids: a list of reference taxids that matched the query
             ref_taxids_df: pandas dataframe from the RefSeq plasmid database with the selected hits based on the search criteria
             stats_host_range_dict: dictionary with taxonomy rank keys and taxonomy rank names as values for future taxonomy hit statistics calculation
    """

    logger.debug("Inside getRefSeqHostRange()")
    logger.debug("Replicon_name_list:{}\nMob_cluster_id:{}\nRelaxase_name_acc_list:{}\nRelaxase_name_list:{}\nMatch_typer:{}\n".format(
                  replicon_name_list,mob_cluster_id_list,relaxase_name_acc_list,relaxase_name_list,matchtype)) #DEBUG




    logger.info("Loading the plasmid reference database ...")

    #load reference lookup table
    ref_taxids_df = pandas.DataFrame()


    logger.info("Searching  RefSeq plasmid database ...")

    if replicon_name_list is not None:
        for replicon_name in replicon_name_list: #allows to run multi-replicon queries
            n_total_records_before = ref_taxids_df.shape[0]
            logger.debug("Replicon name: {} MOB-Cluster: {}\n".format(replicon_name, mob_cluster_id_list))
            ref_taxids_df = pandas.concat([ref_taxids_df,getHostRangeDBSubset(hr_obs_data, replicon_name, "Ref_rep_type(s)", matchtype="exact")])
            logger.debug("Extracted total records (replicon): {}".format(ref_taxids_df.shape[0] - n_total_records_before))
    if mob_cluster_id_list is not None:
        for mob_cluster_id in mob_cluster_id_list:
            ref_taxids_df = pandas.concat([ref_taxids_df, getHostRangeDBSubset(hr_obs_data, mob_cluster_id, "Ref_cluster_id", matchtype="exact")])#must peform only exact match to DB
        logger.debug("Extracted total records (clusterid): {}".format(ref_taxids_df.shape[0]))
    if relaxase_name_acc_list is not None:
        for relaxase_name_acc in relaxase_name_acc_list:
            ref_taxids_df = pandas.concat([ref_taxids_df,getHostRangeDBSubset(hr_obs_data, relaxase_name_acc, "Ref_relaxase_type_accession(s)",matchtype="exact")])
            logger.debug("Relaxase accession: {}\n".format(relaxase_name_acc))
    if relaxase_name_list is not None:
        for relaxase_name in relaxase_name_list:
            ref_taxids_df = pandas.concat([ref_taxids_df, getHostRangeDBSubset(hr_obs_data, relaxase_name, "Ref_relaxase_type(s)", matchtype="exact")])
    if ref_taxids_df.empty:
        logger.warning("RefSeq Plasmid database found no hits ...")
        logger.warning("Search parameters:\n"+
                      "replicon name list: "+str(replicon_name_list)+"; mob_cluster_id: "+
                       str(mob_cluster_id_list)+"; relaxase_name_acc_list: "+
                       str(relaxase_name_acc_list)+"; relaxase_name_list: "+str(relaxase_name_list)+";")
        #print(ref_taxids_df)
        #raise Exception("Empty dataframe returned from RefSeq plasmid database")


    #select subset of ref data based on criteria (e.g. replicon name)
    # if replicon_name_list != None and mob_cluster_id != None and relaxase_name_acc_list == None and relaxase_name_list == None:
    #
    #     #replicon_name="IncI1"; mob_cluster_id="476"
    #     for replicon_name in replicon_name_list:
    #         n_total_records_before = ref_taxids_df.shape[0]
    #         logging.debug("Replicon name: {} MOB-Cluster: {}\n".format(replicon_name, mob_cluster_id))
    #         ref_taxids_df = pandas.concat([ref_taxids_df,getHostRangeDBSubset(hr_obs_data, replicon_name, "Ref_rep_type(s)", matchtype)])
    #         logging.debug("Extracted records (replicon): {}".format(ref_taxids_df.shape[0] - n_total_records_before))
    #
    #     #merge with cluster_id results too
    #     ref_taxids_df = pandas.concat([ref_taxids_df,getHostRangeDBSubset(hr_obs_data, mob_cluster_id, "Ref_cluster_id", matchtype)])
    #     logging.debug("Extracted total records (clusterid): {}".format(ref_taxids_df.shape[0]))
    #
    # elif replicon_name_list != None and mob_cluster_id == None and relaxase_name_acc_list != None and relaxase_name_list == None:
    #     for replicon_name in replicon_name_list:
    #         logging.debug("Replicon name: {}".format(replicon_name))
    #         ref_taxids_df = pandas.concat([ref_taxids_df, getHostRangeDBSubset(hr_obs_data, replicon_name, "Ref_rep_type(s)", matchtype)])
    #     for relaxase_name_acc in relaxase_name_acc_list:
    #         ref_taxids_df = pandas.concat([ref_taxids_df,getHostRangeDBSubset(hr_obs_data, relaxase_name_acc, "Ref_relaxase_type(s)",matchtype)])
    #         logging.debug("Relaxase accession: {}\n".format(relaxase_name_acc))
    #
    # elif replicon_name_list != None and mob_cluster_id == None and relaxase_name_acc_list == None and relaxase_name_list == None:
    #     for replicon_name in replicon_name_list:
    #         n_total_records_before = ref_taxids_df.shape[0]
    #         logging.debug("Replicon name: {}\n".format(replicon_name, mob_cluster_id))
    #         ref_taxids_df = pandas.concat([ref_taxids_df, getHostRangeDBSubset(hr_obs_data, replicon_name, "Ref_rep_type(s)", matchtype)])
    #         logging.debug("Extracted records (replicon): {}".format(ref_taxids_df.shape[0]-n_total_records_before))
    #
    # elif replicon_name_list == None and mob_cluster_id == None and relaxase_name_list != None:
    #     for relaxase_name in relaxase_name_list:
    #         ref_taxids_df = pandas.concat([ref_taxids_df, getHostRangeDBSubset(hr_obs_data, relaxase_name, "Ref_relaxase_type(s)", matchtype)])
    #
    #     ref_taxids_df = getHostRangeDBSubset(hr_obs_data, relaxase_name_acc, "Ref_relaxase_type(s)", matchtype)
    # elif replicon_name_list == None and mob_cluster_id != None and relaxase_name_acc == None:
    #     print("MOBClusterID: {}".format(mob_cluster_id))
    #     ref_taxids_df = getHostRangeDBSubset(hr_obs_data, mob_cluster_id, "Ref_cluster_id", matchtype)
    # elif relaxase_name_list != None:
    #     print("Relaxase Accession: {}".format(args.relaxase_accession))
    #     ref_taxids_df = getHostRangeDBSubset(hr_obs_data, args.relaxase_accession, "Ref_relaxase_type_accession(s)", matchtype)
    # else:
    #     logging.error("No search paramenters specified (i.e. replicon, mob_cluster_id, relaxase family, relaxase accession #).")
    #     exit("Missing ref plasmid database search parameters (e.g. replicon, relaxases, clusterID)")

    unfiltered_n_ref_hits = ref_taxids_df.shape[0]
    #QC: remove uncultured reference database hits as they thend to inflate host range

    if ref_taxids_df.empty == False:
        select_vector=[True if len(re.findall("uncultured",line)) == 0 else False for line in ref_taxids_df['lineage_names'] ]
        ref_taxids_df=ref_taxids_df[select_vector]


    # Check if the search returned no results. Abort
    n_ref_rep_hits = ref_taxids_df.shape[0]  # number of records on the query replicon/relaxase
    logger.debug("QC removed {} hits (uncultured bacteria filter)".format(unfiltered_n_ref_hits - n_ref_rep_hits))

    if n_ref_rep_hits == 0:
        logger.warning("No plasmid reference database hits were found ...")
        #exit("ERROR: No plasmid reference database hits were found ... Try --loose_match parameter")
        convergance_rank="-"; converged_taxonomy_name="-"; unique_ref_selected_taxids="-";
        ref_taxids_df=pandas.DataFrame(); stats_host_range_dict={}
        return(convergance_rank, converged_taxonomy_name, unique_ref_selected_taxids, ref_taxids_df,stats_host_range_dict)


    #do some stats on the returned records
    ref_taxids = list(ref_taxids_df['taxid']) #might be duplicated taxids
    unique_ref_selected_taxids = set(list(ref_taxids_df['taxid'])) #UNIQUE taxids returned as per query

    logger.debug("Found {} records (with duplicates) in the reference database.({} unique and {} duplicated)".format(n_ref_rep_hits,
                                                                                                                 len(unique_ref_selected_taxids),
                                                                                                                 n_ref_rep_hits-len(unique_ref_selected_taxids)))

    convergance_rank, converged_taxonomy_name = getHostRangeRankCovergence(ref_taxids)
    #ref_taxids_str = ",".join(str(t) for t in sorted(unique_ref_taxids))


    #counts_per_ref_taxid_dict = Counter(ref_taxids_df['taxid']) #Counter({562: 66, 624: 10, 573: 6} #where key is taxid


    # #read selected taxonomy lineages from the reference file

    #create dictionary of ranks and taxnomy names for each reference database extracted lineage
    taxonony_lineages_dict_list = list()
    for key in ref_taxids:
         selection_idx = list(ref_taxids_df['taxid'] == key)
         #taxid2species = ref_taxids_df.loc[selection_idx]['Organism'].iloc[0]
         #print("Taxid:{}\tHits:{}\tOrganism:{}".format(key,value,taxid2species))
         tax_lineage_name_list = (ref_taxids_df.loc[selection_idx]['lineage_names'].iloc[0]).split(";")
         tax_lineage_rank_list = (ref_taxids_df.loc[selection_idx]['lineage_ranks'].iloc[0]).split(";")
         #print(tax_lineage_rank_list);print(tax_lineage_name_list)
         taxonony_lineages_dict_list.append(OrderedDict(zip(tax_lineage_rank_list,tax_lineage_name_list)))
    #
    ranks=['superkingdom','phylum', 'class', 'order', 'family', 'genus', 'species'] #important to have them in hierarchical order
    #create lineage statistics dictionary per rank (e.g. 'superkingdom', ['Bacteria', 'Bacteria' ...],'phylum',[...])
    #print(taxonony_lineages_dict_list)
    stats_host_range_dict=OrderedDict({rank:list() for rank in ranks})
    for lineage_record in taxonony_lineages_dict_list:
         for rank in stats_host_range_dict.keys():
             stats_host_range_dict[rank].append(lineage_record.get(rank))

    #
    #print(stats_host_range_dict)

    # missing_stats={}
    # for key in stats_host_range_dict.keys():
    #     taxa_list = stats_host_range_dict[key]
    #     #print(taxa_list)
    #     if any([i== None for i in taxa_list]):
    #         missing_stats[key]=Counter(taxa_list)[None]/len(taxa_list)
    #     else:
    #         missing_stats[key]=0
    #
    # #print("MISSING RATES PER TAXONONIC RANK:\n", missing_stats)
    #
    # #get predicted host range rank
    # convergance_rank="NA"; converged_taxonomy_name="NA"
    # bottom_up_tax_ranks_list=[i for i in reversed(stats_host_range_dict.keys())]
    # for tax_rank in bottom_up_tax_ranks_list: #rank:list(item names)
    #     value = stats_host_range_dict[tax_rank]
    #     n_unique_items = len(set(value))
    #     unique_items = list(set(value))
    #
    #     if  (any([i == None for i in stats_host_range_dict[tax_rank]]) == True and n_unique_items == 2) or \
    #         (n_unique_items == 1 and any([i != None for i in unique_items])):
    #         convergance_rank = tax_rank
    #         converged_taxonomy_name=list(set([i for i in value if i is not None]))[0] #host range is calculated here
    #         break
    #
    #
    # logging.info("HOST RANGE RESULT: {} (level) {} (rank)\n".format(converged_taxonomy_name, convergance_rank))
    #
    # #print("HOST RANGE:{} RANK LEVEL:{}\n".format(converged_taxonomy_name, convergance_rank))
    #
    # #for rank in ranks:
    # #    print("{}: {}".format(rank,dict(Counter(stats_host_range_dict[rank]))))

    return(convergance_rank, converged_taxonomy_name, unique_ref_selected_taxids, ref_taxids_df,stats_host_range_dict)

def taxonomical_hits_breakdown_stats_per_taxonomical_rank(filename_prefix,stats_host_range_dict,dbtype):
    with open(file=filename_prefix+"_refseqhostrange_phylostats.txt", mode="w", encoding="utf-8") as fp:
        strings2file = ["rank\tsci_name\tdb_hits\tconvergance_rank\tconvergance_sci_name\n"]
        for rank in stats_host_range_dict.keys():
            names = (Counter(stats_host_range_dict[rank]).keys())
            values = (Counter(stats_host_range_dict[rank]).values())
            for couple in zip(names,values):
               strings2file.append("{}\t{}\t{}\n".format(rank,couple[0],couple[1]))
        fp.writelines(strings2file)
        fp.close()
        logger.info("Wrote phylogeny stats into {}".format(filename_prefix + "_"+dbtype+"_hostrange_tree_phylostats.txt"))

def writeOutHostRangeReports(   filename_prefix = None,
                                samplename = None,
                                replicon_name_list = None,
                                mob_cluster_id_list = None,
                                relaxase_name_acc_list = None,
                                relaxase_name_list = None,
                                convergance_rank = None,
                                convergance_taxonomy = None,
                                stats_host_range_dict = None,
                                literature_hr_report = pandas.DataFrame()
                                ):
    ###Should refactor this function to be more abastract, Rank and converaganice should accept both datbases predictions on rank
    dict_molecular_features={"replicons":replicon_name_list,"mob_cluster_ids":mob_cluster_id_list,
                             "relaxase_names":relaxase_name_list, "relaxase_name_accs":relaxase_name_acc_list }

    #collapse data if literature report contains more than one cell

    if literature_hr_report.shape[0] > 1:
        literature_hr_report.to_csv(filename_prefix + '_literature_uncollapsed_report.txt', sep="\t", index=False, na_rep="NA",mode="w")
        literature_hr_report = collapseLiteratureReport(literature_hr_report) #collapse report


    for key in dict_molecular_features.keys():
        if dict_molecular_features[key] != None:
            dict_molecular_features[key] = ",".join([str(x) for x in dict_molecular_features[key]])

    with open(filename_prefix+"_refseqhostrange_report.txt",mode="w", encoding="utf-8") as fp:
        strings2file = ["filename\tquery_replicons\tquery_mob_cluster_ids\tquery_relaxase_names\tquery_relaxase_name_accs\tconvergance_refseq_rank\tconvergance_refseq_sci_name\n"]
        strings2file.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                                  samplename,
                                                  dict_molecular_features["replicons"],
                                                  dict_molecular_features["mob_cluster_ids"],
                                                  dict_molecular_features["relaxase_names"],
                                                  dict_molecular_features["relaxase_name_accs"],
                                                  convergance_rank,convergance_taxonomy))
        fp.writelines(strings2file)
        logger.info("Wrote phylogeny stats into {}".format(filename_prefix+"_refseqhostrange_report.txt"))
    fp.close()

    #if no_header_flag is True: #flag to write header in the report file
    #    with open(file=filename_prefix + "_hostrange_refseqhostrange_report.txt", mode="w") as fp:
    #        fp.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("QueryReplicon(s)", "QueryClusterID", "QueryRelaxase(s)",
    #                                                       "QueryRelaxaseAccession", "RefSeqRank", "RefSeqHostRange"))

    #with open(file=filename_prefix+"_hostrange_refseqhostrange_report.txt", mode="a") as fp:
    #        fp.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(replicons, mob_cluster_id,
    #                                                       relaxases, relaxase_name_acc, convergance_rank,
    #                                                       convergance_taxonomy))
    #fp.close()
    #logging.info("Wrote RefSeq host range report results into {}".format(
    #     filename_prefix + "_hostrange_refseqhostrange_report.txt"))  # replicon,relaxase,cluster,host range



    #decompose resulting phylogenetic tree into frequency counts per reach taxonomic rank (hierachy level)
    #strings2file = ["rank\tsci_name\tRefSeq_db_hits\tQueryReplicon(s)\tQueryClusterID\tQueryRelaxase(s)\tQueryRelaxaseAccession\tRefSeqRank\tRefSeqHostRange\n"]
    taxonomical_hits_breakdown_stats_per_taxonomical_rank(filename_prefix,stats_host_range_dict,dbtype="refseqhostrange")


    # literature report
    if literature_hr_report.shape[0] != 0:
        literature_hr_report.to_csv(filename_prefix + '_literature_report.txt', sep="\t", index=False, na_rep="NA", mode="w")



def getTaxonomyTree(taxids):
    """
    Generate and render taxonomical tree as text or image based on the list of taxonomic ids (taxids)
    :param taxids: list of NCBI Taxonomy ids (e.g. 562 - E.coli)
    :return: tree: object of PhyloTree class from ete3 library
    """

    if not isETE3DBTAXAFILEexists():
        logger.info("Did not find taxa.sqlite in {}. Initializaing ete3 taxonomy database".format(ETE3DBTAXAFILE))
        initETE3Database()

    ncbi = NCBITaxa(dbfile=ETE3DBTAXAFILE)
    if not isETE3DBTAXAFILEexists():
        logger.error("Tried ete3 init, but still was not able to find taxa.sqlite file for ete3 lib in {}. Aborting".format(ETE3DBTAXAFILE))
        exit(-1)

    tree = ncbi.get_topology(taxids)
    #prune tree
    #print(tree.get_ascii(attributes=["rank","name","sci_name"])) ;
    #need to prune tree as extra taxids are returned even at the species level obstructing the freq hit counts
    for node in tree.traverse():
        #print(node.rank, node.taxid)
        if any([t for t in taxids if node.taxid == t]):
            #print("Found species node")
            for c in node.get_children():
                #print("Delete child",c.sci_name,c.taxid)
                c.delete(prevent_nondicotomic=False)

        #if node.is_leaf() and any([ t for t in taxids if t == node.taxid]) == False:
        #   print("Error")
        #if node.rank == "no rank":
        #    node.delete()

    return(tree)

def loadHostRangeDB():
    database_abs_path = os.path.dirname(os.path.abspath(__file__))+"/databases/"+"host_range_ncbirefseq_plasmidDB.csv"
    data_obs_hr = pandas.read_csv(database_abs_path, sep=",", encoding="ISO-8859-1",dtype={'Ref_cluster_id':str,'taxid':str})
    return data_obs_hr


def getHostRangeDBSubset(hr_obs_data, search_name, column_name, matchtype):
    """
    From the NCBI RefSeq Plasmid database get a subset as per search query (e.g. replicon, clusterid, relaxase)
    Use different replicon matching methods (exact, multiple, loose). We recommend loose matching due to notation diversity
    :param hr_obs_data: host range observed data from RefSeq database NCBI
    :param search_name: a search string containing a single item (replicon, relaxase, clusterid name)
    :param column_name: a search string
    :param matchtype: a sting
    :return: dataframe containing selected NCBI plasmid information (accession number, mobility, size, etc.)
    """
    #print(matchtype,column_name,search_name,hr_obs_data.loc[0:5,column_name])
    selection_index = list()

    if matchtype == "exact":
        for item in hr_obs_data[column_name]:
            if str(item) == str(search_name):
                selection_index.append(True)
            else:
                selection_index.append(False)
        return (hr_obs_data.loc[selection_index])
    #elif matchtype == "multiexact":
    #    for item in hr_obs_data[column_name]:
    #        if any([i for i in str(item).split(',') if i == search_name]):
    #            selection_index.append(True)
    #        else:
    #            selection_index.append(False)
    #    return (hr_obs_data.loc[selection_index])
    elif matchtype == "loose":  # will losely match for query cases such as IncF* and complex entires (e.g. IncF/IncH)
        for item in hr_obs_data[column_name]:
            if len(re.findall(str(search_name)+".*", str(item))) != 0:
                selection_index.append(True)
            else:
                selection_index.append(False)
        return (hr_obs_data.loc[selection_index])
    else:
        logging.error("Incorrect match type specified (could be either exact (--exact_match) or multi (--loose_match))")
        exit(1)

#parse cml arguments if run as separate module
def parse_args():
    parser = ArgumentParser(description="Welcome to mob-suite plasmid host range prediction module :-)")
    parser.add_argument('--exact_match',action='store_const', const="True", required=False, help='Single exact search criteria match')
    #parser.add_argument('--multi_match', action='store_const', const="True", required=False, help='Multiple exact search criteria match')
    parser.add_argument('--loose_match', action='store_const', const="True", required=False,help='Wildcard loose search criteria match (e.g. IncF* will match IncFI, IncFII, etc.)')
    parser.add_argument('--replicon_name',  action='store', nargs=1, required=False, help='Replicon name(s)')
    parser.add_argument('--relaxase_name', action='store', nargs=1, required=False, help='Relaxase name')
    parser.add_argument('--relaxase_accession', action='store', nargs=1, required=False, help='Relaxase accession number')
    parser.add_argument('--cluster_id', action='store', nargs=1, required=False, help='MOB-Suite Cluster ID (e.g. 416)')
    parser.add_argument('--host_range_detailed', required=False, help='Complete host range report with phylogeny stats',
                        action='store_true',
                        default=False)
    #parser.add_argument('--render_tree_image', required=False, help='Render host range phylogenetic tree in PNG format (requires running X11 server)',
    #                    action='store_true',
    #                    default=False)
    parser.add_argument('-d', '--database_directory',
                        default=default_database_dir,
                        help='Directory where all download databases are located. Defaults to {}'.format(
                            default_database_dir))
    parser.add_argument('--outdir', action='store', required=True, help='Output files name prefix')
    parser.add_argument('--inputseq', action='store', required=False, help='Single plasmid sequence in FASTA format (optional)')
    parser.add_argument('--debug', required=False, help='Show debug detailed information (optional)', action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    args = parser.parse_args()


    #if multiple replicons, relaxases, relaxases accessions, clusterids are present, then we convert them into a list
    if args.replicon_name:
        args.replicon_name = re.split(",", args.replicon_name[0])
    elif args.relaxase_name:
        args.relaxase_name = re.split(",", args.relaxase_name[0])
    elif args.relaxase_accession:
        args.relaxase_accession = re.split(",", args.relaxase_accession[0])
    elif args.cluster_id:
        args.cluster_id = re.split(",", args.cluster_id[0])


    #CASE1: prohibited characters in the name
    #correct file names that come with the prohibited characters like / or \ (e.g. IncA/C2)
    args.outdir = re.sub("[\\|//]+", "", args.outdir)

    #CASE2: Output file exists. Remove the old version
    #if os.path.exists(args.outdir+".txt"):
    #    os.remove(args.outdir+".txt")

    if args.debug:
        logger.setLevel(logging.DEBUG)
        #logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]", level=)
    #print(parser.parse_args())
    #print(parser.parse_args().output)
    #exit()
    return(args)

#def getHostRange_module_enty_point(replicon_name,  mob_cluster_id, relaxase_name_acc, relaxase_name):
#    global args

def main():
    """
    Main entry point into the host-range module
    :param   inputs   replicon name (e.g. IncI,ColI), MOB-suite plasmid clusterID, single plasmid sequence
    :return: reports  text-based reports on NCBI and literature reported host-range and experimentally reported plasmid range
    """
    global args
    args = parse_args()


    logger.info("Running host range main function")
    logger.info("Parsing command line arguments")


    #determine type of the NCBI reference database match
    if args.exact_match != None:
        matchtype="exact"    # single ORF replicon match
    elif args.loose_match != None:
        matchtype="loose"    # approximate ORF replicon match such as replicon family name (e.g. F)
    #elif args.multi_match != None:
    #    matchtype="multiexact"     # multiple ORFs replicon match
    else:
        matchtype=None

    logger.info("INPUT: Replicon=" + str(args.replicon_name) + " and MOB-SuiteClusterID=" + str(args.cluster_id) + ";")
    logger.info("Started to run the main taxonomy query function per feature")

    # hostrange based on MOB-Suite  RefSeq database
    (rank, host_range, taxids, taxids_df, stats_refseq_host_range_dict) = getRefSeqHostRange(args.replicon_name, args.cluster_id,
                                                                                 args.relaxase_name,
                                                                                 args.relaxase_accession, matchtype,
                                                                                 loadHostRangeDB())
    if len(taxids) == 1:
        rank="NA"
        host_range="NA"
        logger.warning("Only single hit was found in the NCBI RefSeq database. Can not predict host range ... Broaden your search criterion")

    if len(taxids) > 0 and args.host_range_detailed:
        treeRefSeq = getTaxonomyTree(taxids)  # get a phylogenetic tree
        renderTree(tree=treeRefSeq,
                   filename_prefix=args.outdir + "/mob_hostrange_" + args.outdir + "_refseqhostrange_")

    #get literature based host range for each replicon
    lit_report=pandas.DataFrame()
    littaxids = []
    if args.replicon_name != None:
        lit_report, littaxids = getLiteratureBasedHostRange(args.replicon_name, loadliteratureplasmidDB())
    if lit_report.empty == False and args.host_range_detailed:
        treeLiterature = getTaxonomyTree(littaxids)
        renderTree(tree=treeLiterature,
                   filename_prefix=args.outdir+"/mob_hostrange_"+args.outdir+"_literaturehostrange_")


    writeOutHostRangeReports( filename_prefix=args.outdir+"/mob_hostrange_"+args.outdir,
                              samplename="-" ,
                              replicon_name_list = args.replicon_name,
                              mob_cluster_id_list = args.cluster_id,
                              relaxase_name_acc_list = args.relaxase_accession,
                              relaxase_name_list = args.relaxase_name,
                              convergance_rank = rank, convergance_taxonomy = host_range,
                              stats_host_range_dict = stats_refseq_host_range_dict,
                              literature_hr_report=lit_report)

    logger.info("Host Range module run is complete!")

def renderTree(tree,filename_prefix):

    with open(file=filename_prefix+ "asci_tree.txt", mode="w", encoding="utf-8") as fp:
        fp.write(tree.get_ascii(attributes=["rank", "sci_name"]))
    logger.info("Wrote ASCII host range tree into {}".format(filename_prefix + "asci_tree.txt"))
    tree.write(format=2, outfile=filename_prefix + "phylogeny_tree.nwk")
    logger.info("Wrote Newick host range tree into {}".format(filename_prefix + "phylogeny_tree.nwk"))



if __name__ == "__main__":
    # setup the application logging
    main()
    logger.info("Run completed")
    #getTaxidsPerRelaxase()
    #IncI2


#TODO
# Bug #1: mob_recon error. mob_aggregated_report.txt is empty sometimes
# Bug #2: mob_recon error. double counting of identical contigs and inflation of the plasmid size.
# Bug #3: mob_hostrange. NCBI RefSeq results are not printed when run as a standalone module. Output results for NCBI RefSeq prediction too
# Bug #4: mob_hostrange. Output NCBI RefSeq + Literature predictions in a single file
# Feature #1: resolve multi-replicon inputs case for NCBI RefSeq and Literature host-range predictions
# Feature #2: Average size of the plasmids of input Inc group(s)
# Feature #3: Match Inc types in Replicon column via regular expressions as some have a complex multi-replicon structure (e.g. ColE2/ColRNAI).
# Match Col-like plasmids with just Col query
# Feature #4: Match Inc types in step-wise manner (subfamily --> family). If no matches returned by subfamily match (e.g. IncFII) match with family (IncF)
# Feature #5: add stats of the host range taxonomy hits for literature host range hits (similar to refseq ones)
# Feature #6: make column heads easier to report and maintain. refactor that section


# BATCH mode processing
# Add phylogenetic stats for the literature inferred tree. Need to add taxonomy lineages full path in literature database
# Add literature closest transfer rate value if available field in all outputs

# mob_typer run
# cp ~/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/*.py /Users/kirill/miniconda/envs/mob_suite_test/lib/python3.6/site-packages/mob_suite/
# mob_typer -i /Users/kirill/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/tests/TestData/IncF/ET11_Ecoli_plasmid_529.fasta --host_range_detailed -o run_test
# mob_typer -i /Users/kirill/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/tests/TestData/KU295134.fasta  --host_range_detailed -o run_test

# Copy latest database files
# cp host_range_literature_plasmidDB_latest.csv /Users/kirill/miniconda/envs/mob_suite_test/lib/python3.6/site-packages/mob_suite/databases/
# cp host_range_ncbirefseq_plasmidDB_latest.csv /Users/kirill/miniconda/envs/mob_suite_test/lib/python3.6/site-packages/mob_suite/databases/
# cp literature_mined_plasmid_seq_db.fasta.msh /Users/kirill/miniconda/envs/mob_suite_test/lib/python3.6/site-packages/mob_suite/databases/