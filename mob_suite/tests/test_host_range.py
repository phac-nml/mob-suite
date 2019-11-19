from mob_suite.mob_host_range import getRefSeqHostRange, getLiteratureBasedHostRange, loadliteratureplasmidDB, \
    loadHostRangeDB, getTaxonomyTree, renderTree
import ete3, numpy, os, logging
import locale

logger=logging.getLogger()
logger.setLevel(logging.DEBUG)


#test the main function
def test_loadliteratureHostRangeDB():
    lit_database_data = loadliteratureplasmidDB()
    print(type(lit_database_data["Size"].dtype))
    assert isinstance(lit_database_data["PMID"].values[0], object), "Incorrect value type for PMID column"
    assert isinstance(lit_database_data["Size"].values[0], numpy.float64), "Incorrect value type for plasmid Size column"
    assert isinstance(lit_database_data["TransferRate"].values[0], numpy.float64), "Incorrect value type for TransferRate column"

def test_HostRangeDB():
    refseq_plasmid_database_data = loadHostRangeDB()
    assert isinstance(refseq_plasmid_database_data["Ref_cluster_id"].values[0],object), "Incorrect value type for PMID column"
    assert isinstance(refseq_plasmid_database_data["lineage_taxid"].values[0], object), "Incorrect value type for lineage_taxid column"
    assert isinstance(refseq_plasmid_database_data["lineage_ranks"].values[0], object), "Incorrect value type for lineage_ranks column"
    assert isinstance(refseq_plasmid_database_data["lineage_names"].values[0], object), "Incorrect value type for lineage_names column"
    assert isinstance(refseq_plasmid_database_data["taxid"].values[0], object), "Incorrect value type for lineage_ranks column"

def test_getRefSeqHostRange_replicon():
    """
    Check the RefSeq host range prediction functionality using IncI1 as a query replicon.
    Test for host range correct rank, rank scientific name, hits dataframe columns and hits stats dictionary presence
    """
    ref_taxids_df_columns = ["Ref_accession", "Ref_total_length", "Ref_gc", "Ref_rep_type(s)",
                             "Ref_rep_type_accession(s)",
                             "Ref_relaxase_type(s)", "Ref_relaxase_type_accession(s)", "Ref_mpf_type",
                             "Ref_mpf_type_accession(s)",
                             "Ref_orit_type(s)", "Ref_orit_accession(s)", "Ref_PredictedMobility", "Ref_cluster_id",
                             "Description", "Organism", "Conformation", "Molecule", "plasmid_id", "taxid",
                             "lineage_taxid",
                             "lineage_ranks", "lineage_names"]

    for matchtype in ["exact","loose"]:
        (returnhrrank,returnhrsciname,unique_ref_selected_taxids,ref_taxids_df,stats_host_range_dict) = getRefSeqHostRange(
                           replicon_name_list = ["IncI1"],
                           mob_cluster_id_list= None,
                           relaxase_name_acc_list=None,
                           relaxase_name_list=None,
                           matchtype=matchtype,
                           hr_obs_data=loadHostRangeDB())

        assert returnhrrank == "order", "Rank did not match. The output rank="+returnhrrank+" and expected=order. Match type="+matchtype
        assert returnhrsciname == "Enterobacterales", "host range sci name did not match. The output host range="+returnhrsciname+" and expected=Enterobacterales. Match type="+matchtype
        assert unique_ref_selected_taxids != [], "Could not extract any hits from the RefSeq Host range database. Taxid list is empty. Check the getRefSeqHostRange()"
        assert ref_taxids_df.columns.tolist() == ref_taxids_df_columns, "RefSeq plasmid database column names do not match reference list anymore"
        assert stats_host_range_dict != {}, "Empty RefSeq host-range hits dictionary. Check the getRefSeqHostRange()"

    #test the no hits returned by the RefSeq host range prediction
    #with pytest.raises(Exception) as e:
    (returnhrrank, returnhrsciname, unique_ref_selected_taxids, ref_taxids_df, stats_host_range_dict) = getRefSeqHostRange(
            replicon_name_list=None,
            mob_cluster_id_list=None,
            relaxase_name_acc_list=None,
            relaxase_name_list=None,
            matchtype="exact",
            hr_obs_data=loadHostRangeDB())
    assert returnhrrank == "-", "Wrong value. Returned "+returnhrrank
    assert returnhrsciname == "-", "Wrong value. Returned "+returnhrsciname
    assert stats_host_range_dict == {},"Wrong value. Returned "+stats_host_range_dict
    assert ref_taxids_df.empty == True


    #assert str(e.value) == "Empty dataframe returned from RefSeq plasmid database", "Error message: "+str(e.value)

    #test when there are hits for RefSeq host-range search only (e.g. ColRNAI_rep_cluster_1857)
    getRefSeqHostRange(
        replicon_name_list=["ColRNAI_rep_cluster_1857"],
        mob_cluster_id_list=None,
        relaxase_name_acc_list=None,
        relaxase_name_list=None,
        matchtype="exact",
        hr_obs_data=loadHostRangeDB())

def test_getRefSeqHostRange_multireplicon():
    """
    Test multi-replicon case when multiple replicons are queried simultaneously. The prediction results are aggregated together together
    :return:
    """
    for matchtype in ["exact","loose"]:
        (returnhrrank,returnhrsciname,unique_ref_selected_taxids,ref_taxids_df,stats_host_range_dict) = getRefSeqHostRange(
                           replicon_name_list = ["IncI1","IncI2"],
                           mob_cluster_id_list = None,
                           relaxase_name_acc_list = None,
                           relaxase_name_list = None,
                           matchtype=matchtype,
                           hr_obs_data=loadHostRangeDB())
        assert returnhrrank == "order", "Rank did not match. The output rank="+returnhrrank+" and expected=order. Match type="+matchtype
        assert returnhrsciname == "Enterobacterales", "host range sci name did not match. The output host range="+returnhrsciname+" and expected=Enterobacterales. Match type="+matchtype
        assert unique_ref_selected_taxids != [], "Could not extract any hits from the RefSeq Host range database. Taxid list is empty. Check the getRefSeqHostRange()"
        assert stats_host_range_dict != {}, "Empty RefSeq host-range hits dictionary. Check the getRefSeqHostRange()"

    for matchtype in ["exact", "loose"]:
        (returnhrrank, returnhrsciname, unique_ref_selected_taxids, ref_taxids_df,
         stats_host_range_dict) = getRefSeqHostRange(
            replicon_name_list=["IncFIIA", "IncFII"],
            mob_cluster_id_list=None,
            relaxase_name_acc_list=None,
            relaxase_name_list=None,
            matchtype=matchtype,
            hr_obs_data=loadHostRangeDB())


def test_getHostrange_clusterID():
    (convergance_rank, converged_taxonomy_name,
     unique_ref_selected_taxids, ref_taxids_df,
     stats_host_range_dict) = getRefSeqHostRange(replicon_name_list=None,
                       mob_cluster_id_list=[476],  #IncI1 plasmid cluster
                       relaxase_name_acc_list=None,
                       relaxase_name_list=None,
                       matchtype="exact",
                       hr_obs_data=loadHostRangeDB())
    assert convergance_rank == "family", "The host range rank is incorrect. Reported "+ convergance_rank
    assert converged_taxonomy_name == "Enterobacteriaceae", "The host range rank scientific name is incorrect. Reported " + converged_taxonomy_name
    assert unique_ref_selected_taxids != None, "Empty taxonomy id list."+unique_ref_selected_taxids+" Check search criteria"

def test_getHostRange_cluster_and_replicon():
    (convergance_rank, converged_taxonomy_name,
     unique_ref_selected_taxids, ref_taxids_df,
     stats_host_range_dict) = getRefSeqHostRange(replicon_name_list=['rep_cluster_101'],
                                                 mob_cluster_id_list=[11171],  # Non-existing cluster
                                                 relaxase_name_acc_list=None,
                                                 relaxase_name_list=None,
                                                 matchtype="loose_match",
                                                 hr_obs_data=loadHostRangeDB())
    print(convergance_rank, converged_taxonomy_name,
     unique_ref_selected_taxids, ref_taxids_df,
     stats_host_range_dict)



def test_getHostrange_clusterIDs():
    """
    mob_hostrange --loose_match --cluster_id "369,592,370,3558"  --host_range_detailed --outdir "run_test"
    :return:
    """
    (convergance_rank, converged_taxonomy_name,
     unique_ref_selected_taxids, ref_taxids_df,
     stats_host_range_dict) = getRefSeqHostRange(replicon_name_list=None,
                       mob_cluster_id_list=[369,592,370,3558],  #IncFII plasmid cluster
                       relaxase_name_acc_list=None,
                       relaxase_name_list=None,
                       matchtype="exact",
                       hr_obs_data=loadHostRangeDB())
    assert convergance_rank == "order", "The host range rank is incorrect. Reported "+ convergance_rank
    assert converged_taxonomy_name == "Enterobacterales", "The host range rank scientific name is incorrect. Reported " + converged_taxonomy_name
    assert unique_ref_selected_taxids != None, "Empty taxonomy id list."+unique_ref_selected_taxids+" Check search criteria"


def test_getHostrange_relaxase_acc():
    (convergance_rank, converged_taxonomy_name,
     unique_ref_selected_taxids, ref_taxids_df,
     stats_host_range_dict) = getRefSeqHostRange(replicon_name_list=None,
                       mob_cluster_id_list=None,
                       relaxase_name_acc_list=["NC_001735_00050"],
                       relaxase_name_list=None,
                       matchtype="multiexact",
                       hr_obs_data=loadHostRangeDB())
    assert convergance_rank == "superkingdom", "The host range rank is incorrect. Reported "+ convergance_rank
    assert converged_taxonomy_name == "Bacteria", "The host range rank scientific name is incorrect. Reported " + converged_taxonomy_name
    assert unique_ref_selected_taxids != None, "Empty taxonomy id list."+unique_ref_selected_taxids+" Check search criteria"

    (convergance_rank, converged_taxonomy_name,
     unique_ref_selected_taxids, ref_taxids_df,
     stats_host_range_dict) = getRefSeqHostRange(replicon_name_list=None,
                        mob_cluster_id_list = None,
                        relaxase_name_acc_list = ["NC_009739_00028","NC_016978_00054"],
                        relaxase_name_list = None,
                        matchtype = "multiexact",
                        hr_obs_data = loadHostRangeDB())
    assert convergance_rank == "phylum", "The host range rank is incorrect. Reported " + convergance_rank
    assert converged_taxonomy_name == "Proteobacteria", "The host range rank scientific name is incorrect. Reported " + converged_taxonomy_name
    assert unique_ref_selected_taxids != None, "Empty taxonomy id list." + unique_ref_selected_taxids + " Check search criteria"


def test_getRefSeqHostRange_relaxase_name():
    (convergance_rank, converged_taxonomy_name,
     unique_ref_selected_taxids, ref_taxids_df,
     stats_host_range_dict) = getRefSeqHostRange(replicon_name_list = None,
                                                 mob_cluster_id_list = None,
                                                 relaxase_name_acc_list=None,
                                                 relaxase_name_list=["MOBC"],
                                                 matchtype="multiexact",
                                                 hr_obs_data=loadHostRangeDB())

    assert convergance_rank == "superkingdom", "The host range rank is incorrect. Reported " + convergance_rank
    assert converged_taxonomy_name == "Bacteria", "The host range rank scientific name is incorrect. Reported " + converged_taxonomy_name
    assert unique_ref_selected_taxids != None, "Empty taxonomy id list." + unique_ref_selected_taxids + " Check search criteria"

    (convergance_rank, converged_taxonomy_name,
     unique_ref_selected_taxids, ref_taxids_df,
     stats_host_range_dict) = getRefSeqHostRange(replicon_name_list=None,
                                                 mob_cluster_id_list = None,
                                                 relaxase_name_acc_list=None,
                                                 relaxase_name_list=["MOBV","MOBC"],
                                                 matchtype="multiexact",
                                                 hr_obs_data=loadHostRangeDB())
    assert convergance_rank == "superkingdom", "The host range rank is incorrect. Reported " + convergance_rank
    assert converged_taxonomy_name == "Bacteria", "The host range rank scientific name is incorrect. Reported " + converged_taxonomy_name
    assert unique_ref_selected_taxids != None, "Empty taxonomy id list." + unique_ref_selected_taxids + " Check search criteria"

def test_literature_hostrange_single_replion():
    """
    replicon_names,plasmid_lit_db,input_seq=""
    Test the literature-based plasmid host range prediction and transfer rate
    """
    lit_report_columns = ["LiteratureQueryReplicon", "LiteratureSearchReplicon","LiteratureFoundPlasmidsNames",
    "LiteratureFoundPlasmidsNumber","LiteratureReportedHostRangePlasmidClass","LiteratureReportedHostPlasmidSpecies",
    "LiteratureReportedPlasmidHostSpeciesNumber","LiteraturePredictedHostRangeTreeRank",
    "LiteraturePredictedHostRangeTreeRankSciName","LiteratureReportedHostRangeRankInPubs",
    "LiteratureReportedHostRangeNameInPubs","LiteratureMinTransferRateRange","LiteratureMaxTransferRateRange",
    "LiteratureMeanTransferRateRange","LiteraturePMIDs","LiteraturePublicationsNumber"]

    report,littaxids = getLiteratureBasedHostRange(replicon_names =["IncF"],
                                plasmid_lit_db = loadliteratureplasmidDB(),input_seq="")

    #check the column names in literature report
    for column in report.columns.values:
        assert any([c == column for c in lit_report_columns]), "The column name "+column+" is not it the literature report dataframe"

    assert report.empty == False, "Literature host range prediction is empty. Check your search criteria/criterium"
    #assert report.columns.tolist() == lit_report_columns, "Literature host range report does not match all default columns"
    #assert report.shape[0] == 1, "Literature host range report dimension is incorrect. The expect dimension is a single row (1,15)"
    assert report.loc[0,"LiteraturePredictedHostRangeTreeRank"] == "family", "Literature hit-based host range host range prediction is incorrect. Expected family, got"+report.loc[0,"LiteraturePredictedDBHostRangeTreeRank"]
    assert report.loc[0,"LiteraturePredictedHostRangeTreeRankSciName"] == "Enterobacteriaceae"
    assert report.loc[0,"LiteratureReportedHostRangeRankInPubs"] == "family","Literature reported host range host range prediction is incorrect. Expected family, got"+report.loc[0,"LiteraturePredictedDBHostRangeTreeRank"]

    #test when there are no hits for a given replicon in literature database
    report, littaxids = getLiteratureBasedHostRange(replicon_names=["NA"],
                                                    plasmid_lit_db=loadliteratureplasmidDB(), input_seq="")
    assert report.empty, "The literature results df should be empty but it is not"

def test_literature_hostrange_multi_replion():
    """
    Test the literature-based plasmid host range prediction and transfer rate with multiple query replicons belonging to the same Inc family
    """

    report,littaxids = getLiteratureBasedHostRange(replicon_names =["IncFI","IncFII"],
                                plasmid_lit_db = loadliteratureplasmidDB(),input_seq="")

    assert report.empty == False, "Literature host range prediction is empty. Check your search criteria/criterium"
    assert report.shape[0] == 2, "Literature host range report dimension is incorrect. The expect dimension is double row"
    assert all( report["LiteraturePredictedHostRangeTreeRank"].values == ["genus","family"]), "LiteraturePredictedHostRange seems to be incorrect"
    assert all(report.loc[0,"LiteraturePredictedHostRangeTreeRankSciName"] == ["Salmonella","Enterobacteriaceae"]), "LiteraturePredictedHostRange expressed in sci names seems to be incorrect"

    #test when there are no hits for a given replicon in literature database
    report, littaxids = getLiteratureBasedHostRange(replicon_names=["ColRNAI_rep_cluster_1857"],
                                                    plasmid_lit_db=loadliteratureplasmidDB(), input_seq="")
    assert report.empty, "The literature results df should be empty but it is not"

def test_getTaxonomyTree_onRefSeqDB_and_Render():
    convergance_rank, converged_taxonomy_name, unique_ref_selected_taxids, \
    ref_taxids_df, stats_host_range_dict = getRefSeqHostRange(replicon_name_list=["IncP"],matchtype="exact",
                                               mob_cluster_id_list = None, relaxase_name_acc_list=None,
                                               relaxase_name_list= None, hr_obs_data=loadHostRangeDB())
    tree = getTaxonomyTree(taxids=unique_ref_selected_taxids)
    assert isinstance(tree, ete3.PhyloTree), "Class mismatch output by getTaxonomyTree(). Output class is " + type(tree).__name__

    if os.path.exists("run_test") == False:
        os.mkdir("run_test")
    renderTree(tree=tree, filename_prefix="run_test/run_test_refseqhostrange")

def test_getTaxonomyTree_onLiteratureDB_and_Render():
    report_df, lit_taxids = getLiteratureBasedHostRange(replicon_names=["IncP"],plasmid_lit_db=loadliteratureplasmidDB(),)
    tree = getTaxonomyTree(taxids=lit_taxids)
    assert isinstance(tree,ete3.PhyloTree), "Class mismatch output by getTaxonomyTree(). Output class is "+type(tree).__name__
    renderTree(tree=tree, filename_prefix="run_test/run_test_literaturehostrange")



def test_nonASCII_characters_in_host_range():
    "Check for potential printing errors such as UnicodeEncodeError: 'ascii' codec can't encode characters in position 523-524: ordinal not in range(128)"
    dbplasmids = loadHostRangeDB()

    try:
        print(dbplasmids)
    except UnicodeError as e:
        print(e)
        logger.error("Database has non-ascii characters that might not display well on all locales")
        logger.error("Check that host_range_ncbirefseq_plasmidDB.csv has only ASCII characters")
        print("Your environment:")
        os.system("env")
        exit(-1)

    dblit = loadliteratureplasmidDB()

    assert all([isinstance(value,float) for value in dblit["TransferRate"].values]) == True
    try:
        print(dblit)
    except UnicodeError as e:
        print(e)
        logger.error("Database has non-ascii characters that might not display well on all locales")
        logger.error("Check that host_range_literature_plasmid.csv has only ASCII characters")
        print("Your environment:")
        os.system("env")
        exit(-1)

    print("Your environment:")
    os.system("env")

    #set locale via export LC_ALL = POSIX or C.UTF-8




