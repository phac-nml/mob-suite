import pandas, re, logging, os
from argparse import ArgumentParser
from ete3 import NCBITaxa, TreeStyle
from collections import Counter,OrderedDict

#default init arguments
args=ArgumentParser()
args.multi_match = True
args.exact_match = None
args.loose_match = None
args.render_tree = True
args.write_newick = True
args.header = False
args.debug = False
args.outname = ""

OUTDIR = os.getcwd()+"/"

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


def getHostRange(replicon_name_list,  mob_cluster_id, relaxase_name_acc, relaxase_name_list, matchtype):

    if mob_cluster_id != None:
        mob_cluster_id=str(mob_cluster_id) #make sure it is a string

    logging.info("Loading the plasmid reference database ...")

    hr_obs_data = loadHostRangeDB() #load reference lookup table
    ref_taxids_df = pandas.DataFrame()


    logging.info("Searching plasmid database ...")

    #select subset of ref data based on criteria (e.g. replicon name)
    if replicon_name_list != None and mob_cluster_id != None and relaxase_name_acc == None:

        #replicon_name="IncI1"; mob_cluster_id="476"
        for replicon_name in replicon_name_list:
            logging.debug("Replicon name: {} MOB-Cluster: {}\n".format(replicon_name, mob_cluster_id))
            ref_taxids_df = pandas.concat([ref_taxids_df,getHostRangeDBSubset(hr_obs_data, replicon_name,
                                                                              "Ref_rep_type(s)", matchtype)])

        #merge with cluster_id results too
        ref_taxids_df = pandas.concat([ref_taxids_df, \
                                       getHostRangeDBSubset(hr_obs_data, mob_cluster_id,
                                                            "Ref_culster_id.at.dist-0.05", matchtype)])

    elif replicon_name_list != None and mob_cluster_id == None and relaxase_name_acc != None:
        for replicon_name in replicon_name_list:
            logging.debug("Replicon name: {}".format(replicon_name))
            ref_taxids_df = pandas.concat(
                [ref_taxids_df, getHostRangeDBSubset(hr_obs_data, replicon_name, "Ref_rep_type(s)", matchtype)])

        logging.debug("Relaxase accession: {}\n".format(relaxase_name_acc))
        ref_taxids_df = pandas.concat([getHostRangeDBSubset(hr_obs_data, replicon_name, "Ref_rep_type(s)", matchtype), \
                                       getHostRangeDBSubset(hr_obs_data, relaxase_name_acc, "Ref_relaxase_type(s)", matchtype)])
    elif replicon_name_list != None and mob_cluster_id == None and relaxase_name_acc == None:
        for replicon_name in replicon_name_list:
            logging.debug("Replicon name: {}\n".format(replicon_name, mob_cluster_id))
            ref_taxids_df = pandas.concat(
                [ref_taxids_df, getHostRangeDBSubset(hr_obs_data, replicon_name, "Ref_rep_type(s)", matchtype)])
    elif replicon_name_list == None and mob_cluster_id == None and relaxase_name_list != None:
        for relaxase_name in relaxase_name_list:
            ref_taxids_df = pandas.concat(
                [ref_taxids_df, getHostRangeDBSubset(hr_obs_data, relaxase_name, "Ref_relaxase_type(s)", matchtype)])
        print("Relaxase family: {}\n".format(relaxase_name))
        ref_taxids_df = getHostRangeDBSubset(hr_obs_data, relaxase_name_acc, "Ref_relaxase_type(s)", matchtype)
    elif replicon_name_list == None and mob_cluster_id != None and relaxase_name_acc == None:
        print("MOBClusterID: {}".format(mob_cluster_id))
        ref_taxids_df = getHostRangeDBSubset(hr_obs_data, mob_cluster_id, "Ref_culster_id.at.dist-0.05", matchtype)
    elif relaxase_name_list != None:
        print("Relaxase Accession: {}".format(args.relaxase_accession))
        ref_taxids_df = getHostRangeDBSubset(hr_obs_data, args.relaxase_accession, "Ref_relaxase_type_accession(s)", matchtype)
    else:
        logging.error("No search paramenters specified (i.e. replicon, mob_cluster_id, relaxase family, relaxase accession #).")
        exit("Missing ref plasmid database search parameters (e.g. replicon, relaxases, clusterID)")

    #QC: remove uncultured reference database hits as they thend to inflate host range
    select_vector=[True if len(re.findall("uncultured",line)) == 0 else False for line in ref_taxids_df["lineage_names"] ]
    ref_taxids_df=ref_taxids_df[select_vector]

    # Check if the search returned no results. Abort
    n_ref_rep_hits = ref_taxids_df.shape[0]  # number of records on the query replicon/relaxase
    if n_ref_rep_hits == 0:
        logging.error("ERROR: No hits were found ...")
        exit("ERROR: No plasmid reference database hits were found ...")


    #do some stats on the returned records
    ref_taxids = list(ref_taxids_df['taxid']) #might be duplicated taxids
    unique_ref_selected_taxids = set(list(ref_taxids_df['taxid'])) #UNIQUE taxids returned as per query

    logging.debug("Found {} records (with duplicates) in the reference database.({} unique and {} duplicated)".format(n_ref_rep_hits, \
                                                                                                                 len(unique_ref_selected_taxids),
                                                                                                                 n_ref_rep_hits-len(unique_ref_selected_taxids)))
    #ref_taxids_str = ",".join(str(t) for t in sorted(unique_ref_taxids))


    #counts_per_ref_taxid_dict = Counter(ref_taxids_df['taxid']) #Counter({562: 66, 624: 10, 573: 6} #where key is taxid


    #read selected taxonomy lineages from the reference file
    taxonony_lineages_dict_list=list()

    #create dictionary of ranks and taxnomy names for each reference database extracted lineage
    for key in ref_taxids:
        selection_idx = list(ref_taxids_df['taxid'] == key)
        #taxid2species = ref_taxids_df.loc[selection_idx]['Organism'].iloc[0]
        #print("Taxid:{}\tHits:{}\tOrganism:{}".format(key,value,taxid2species))
        tax_lineage_name_list = (ref_taxids_df.loc[selection_idx]['lineage_names'].iloc[0]).split(";")
        tax_lineage_rank_list = (ref_taxids_df.loc[selection_idx]['lineage_ranks'].iloc[0]).split(";")
        #print(tax_lineage_rank_list);print(tax_lineage_name_list)
        taxonony_lineages_dict_list.append(OrderedDict(zip(tax_lineage_rank_list,tax_lineage_name_list)))

    ranks=['superkingdom','phylum', 'class', 'order', 'family', 'genus', 'species'] #important to have them in hierarchical order

    #create lineage statistics dictionary per rank (e.g. 'superkingdom', ['Bacteria', 'Bacteria' ...],'phylum',[...])
    stats_host_range_dict=OrderedDict({rank:list() for rank in ranks})

    for lineage_record in taxonony_lineages_dict_list:
        for rank in stats_host_range_dict.keys():
            stats_host_range_dict[rank].append(lineage_record.get(rank))


    missing_stats={}
    for key in stats_host_range_dict.keys():
        taxa_list = stats_host_range_dict[key]
        #print(taxa_list)
        if any([i== None for i in taxa_list]):
            missing_stats[key]=Counter(taxa_list)[None]/len(taxa_list)
        else:
            missing_stats[key]=0

    #print("MISSING RATES PER TAXONONIC RANK:\n", missing_stats)

    convergance_rank="NA"; converged_taxonomy_name="NA"
    bottom_up_tax_ranks_list=[i for i in reversed(stats_host_range_dict.keys())]
    for tax_rank in bottom_up_tax_ranks_list: #rank:list(item names)
        value = stats_host_range_dict[tax_rank]
        n_unique_items = len(set(value))
        unique_items = list(set(value))

        if  (any([i == None for i in stats_host_range_dict[tax_rank]]) == True and n_unique_items == 2) or \
            (n_unique_items == 1 and any([i != None for i in unique_items])):
            convergance_rank = tax_rank
            converged_taxonomy_name=list(set([i for i in value if i is not None]))[0] #host range is calculated here
            break


    logging.info("HOST RANGE:{} RANK LEVEL:{}\n".format(converged_taxonomy_name, convergance_rank))

    print("HOST RANGE:{} RANK LEVEL:{}\n".format(converged_taxonomy_name, convergance_rank))

    #for rank in ranks:
    #    print("{}: {}".format(rank,dict(Counter(stats_host_range_dict[rank]))))

    return(convergance_rank, converged_taxonomy_name, unique_ref_selected_taxids, ref_taxids_df,stats_host_range_dict)



def writeOutHostRangeResults(   filename = args.outname,\
                                replicon_name_list = None, \
                                mob_cluster_id = None, \
                                relaxase_name_acc = None, \
                                relaxase_name_list = None, \
                                convergance_rank = None, \
                                convergance_taxonomy = None, stats_host_range_dict = None,\
                                header_flag = None, treeObject=None):


    if replicon_name_list != None:
        replicons = ",".join(replicon_name_list)
    else:
        replicons = replicon_name_list
    if relaxase_name_list  != None:
        relaxases = ",".join(relaxase_name_list )
    else:
        relaxases = relaxase_name_list

    if header_flag is not None:
        with open(file=OUTDIR + filename + "_hostrange_report.txt", mode="w") as fp:
            fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Replicon(s)", "ClusterID", "Relaxase(s)", \
                                                           "RelaxaseAccession", "Rank", "HostRange", "GenusCountStats"))

    with open(file=OUTDIR+filename+"_hostrange_report.txt", mode="a") as fp:
            fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(replicons, mob_cluster_id, \
                                                           relaxases, relaxase_name_acc, convergance_rank, \
                                                           convergance_taxonomy,dict(Counter(stats_host_range_dict["genus"]))))
    fp.close()

    with open(file=OUTDIR + filename + "_asci_tree.txt", mode="w") as fp:
        fp.write(treeObject.get_ascii(attributes=["rank", "sci_name"]))
    fp.close()

    #decompose resulting phylogenetic tree into frequency counts per reach taxonomic rank (hierachy level)
    strings2file = ["rank\tname\tref_db_hits\n"]
    with open(file=OUTDIR+filename+"_phyloStats.txt", mode="w") as fp:
        for rank in stats_host_range_dict.keys():
            names = (Counter(stats_host_range_dict[rank]).keys())
            values = (Counter(stats_host_range_dict[rank]).values())
            for couple in zip(names,values):
               strings2file.append("{}\t{}\t{}\n".format(rank,couple[0],couple[1]))
        fp.writelines(strings2file)
        fp.close()

    logging.info("Wrote host range report results into {}".format(OUTDIR+filename+"_hostrange_report.txt")) #replicon,relaxase,cluster,host range
    logging.info("Wrote phylogeny stats into {}".format(OUTDIR + filename + "_phyloStats.txt"))
    logging.info("Wrote ASCII host range tree into {}".format(OUTDIR + filename + "_asci_tree.txt"))


def getTaxonomyTree(taxids, ref_taxids_df, filename=args.outname):
    ncbi = NCBITaxa()

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

    #exit()


    ts = TreeStyle()
    ts.show_leaf_name = False; ts.show_scale=False; ts.show_branch_length = False;

    #print(dir(ts))
    #print(tree.get_ascii(attributes=["sci_name"]))

    #pandas.DataFrame.to_csv(ref_taxids_df, "/Users/kirill/WORK/MOBSuiteHostRange2018/selected_df.csv")

    for node in tree.traverse():
        #print(node.features) #{'rank', 'common_name', 'support', 'sci_name', 'named_lineage', 'taxid', 'lineage', 'name', 'dist'}
        #print(node.sci_name)
        #print(ref_taxids_df)
        nhits = len(ref_taxids_df.loc[ref_taxids_df["taxid"] == node.taxid])
        #print("{}:{}".format(node.taxid,nhits))


        node.img_style['size'] = nhits
        node.img_style['fgcolor'] = "red"
        node.img_style['hz_line_color'] = "red"
        node.img_style['vt_line_color'] = "red"
        node.img_style['draw_descendants']=True
        #print(node.img_style)
        node.name = re.sub("_","",node.sci_name) + "|taxid"+str(node.taxid)+"|"+str(nhits)+"hits" #Newick does not like spaces in labels
        #node.name=node.sci_name
        #print(len(node.name))
        #print(node.name, node.img_style)
        #node.set_style()
        #print(node.features)
    #exit()
    #tree.show(tree_style=ts)

    if args.render_tree == True:
        tree.render(OUTDIR+filename+"_phylogeny_tree.png",  dpi=2800, w=2000, tree_style=ts)

    if args.write_newick == True:
        tree.write(format=2, outfile=OUTDIR+filename+"_phylogeny_tree.nwk")
    return(tree)

def loadHostRangeDB():
    database_abs_path = os.path.dirname(__file__)+"/databases/"+"host_range_plasmidDB.csv"
    data_obs_hr = pandas.read_csv(database_abs_path, sep=",")
    return data_obs_hr


def getHostRangeDBSubset(hr_obs_data, search_name, column_name, matchtype):
    selection_index = list()
    if matchtype == "exact":
        for item in hr_obs_data[column_name]:
            if str(item) == str(search_name):
                selection_index.append(True)
            else:
                selection_index.append(False)
        return (hr_obs_data.loc[selection_index])
    elif matchtype == "multi":
        for item in hr_obs_data[column_name]:
            if any([i for i in str(item).split(',') if i == search_name]):
                selection_index.append(True)
            else:
                selection_index.append(False)
        return (hr_obs_data.loc[selection_index])
    elif matchtype == "loose":  # will losely match for cases such as IncF*
        for item in hr_obs_data[column_name]:
            if len(re.findall(search_name+".*", item)) != 0:
                selection_index.append(True)
            else:
                selection_index.append(False)
        return (hr_obs_data.loc[selection_index])
    else:
        exit("No match type specified (could be either exact (--exact_match) or multi (--multi_match))")


def parse_args():
    parser = ArgumentParser(description="Welcome to mob-suite Host Range module :-)")
    parser.add_argument('--exact_match',action='store_const', const="True", required=False, help='Exact single replicon/relaxase match')
    parser.add_argument('--multi_match', action='store_const', const="True", required=False, help='Exact multi replicon/relaxase match')
    parser.add_argument('--loose_match', action='store_const', const="True", required=False,help='Wildcard loose replicon/relaxase match')
    parser.add_argument('--replicon_name',  action='store', nargs=1, required=False, help='Replicon name')
    parser.add_argument('--relaxase_name', action='store', nargs=1, required=False, help='Relaxase name')
    parser.add_argument('--relaxase_accession', action='store', required=False, help='Relaxase accession number')
    parser.add_argument('--cluster_id', action='store', required=False, help='MOB-Suite Cluster ID (e.g. 416)')
    parser.add_argument('--render_tree',action='store_true', default=False,  required=False, help='render taxanomic tree')
    parser.add_argument('--write_newick', action='store_true', required=False, help='write a newick tree')
    parser.add_argument('--header', action='store_const', const="header", required=False, help='Print header in the output')
    parser.add_argument('--outname', action='store', required=True, help='output files name prefix')
    parser.add_argument('--debug', required=False, help='Show debug information', action='store_true')


    args = parser.parse_args()
    print(args);exit()

    #CASE1: prohibited characters in the name
    #correct file names that come with the prohibited characters like / or \ (e.g. IncA/C2)
    args.outname = re.sub("[\\|//]+", "", args.outname)

    #CASE2: Output file exists. Remove the old version
    if os.path.exists(OUTDIR+args.outname+".txt"):
        os.remove(OUTDIR+args.outname+".txt")

    if args.debug:
        logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]", level=logging.DEBUG)
    #print(parser.parse_args())
    #print(parser.parse_args().output)
    #exit()
    return(args)

def getHostRange_module_enty_point(replicon_name,  mob_cluster_id, relaxase_name_acc, relaxase_name):
    global args

def main():
    global args
    args = parse_args()

    logging.info("Running host range main function")
    # get all command line arguments
    logging.info("Parsing command line arguments")


    #determine type of the reference database match
    if args.exact_match != None:
        matchtype="exact"
    elif args.loose_match != None:
        matchtype="loose"
    elif args.multi_match != None:
        matchtype="multi"
    else:
        matchtype=None

    logging.info("INPUT: Replicon=" + str(args.replicon_name) + " and MOB-SuiteClusterID=" + str(args.cluster_id) + ";")

    logging.info("Started to run the main taxonomy query function per feature")

    (rank, host_range, taxids, taxids_df,stats_host_range) = getHostRange(args.replicon_name, args.cluster_id,
                                      args.relaxase_name, args.relaxase_accession, matchtype)

    tree = getTaxonomyTree(taxids, taxids_df, args.outname)  # get phylogenetic tree

    writeOutHostRangeResults(
                        filename=args.outname,\
                        replicon_name_list = args.replicon_name, \
                        mob_cluster_id = args.cluster_id, \
                        relaxase_name_acc = args.relaxase_accession, \
                        relaxase_name_list = args.relaxase_name,\
                        convergance_rank = rank, convergance_taxonomy = host_range, \
                        stats_host_range_dict = stats_host_range, header_flag = args.header, treeObject=tree)

if __name__ == "__main__":
    # setup the application logging
    main()
    logging.info("Run completed")
    #getTaxidsPerRelaxase()
    #IncI2


