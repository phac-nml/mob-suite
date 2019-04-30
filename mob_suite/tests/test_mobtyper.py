from mob_suite.mob_host_range import getRefSeqHostRange, getLiteratureBasedHostRange, loadliteratureplasmidDB, \
    loadHostRangeDB, getTaxonomyTree, renderTree
import mob_suite.mob_typer
import pytest,os,sys
import ete3, pandas,argparse


#test the main function
def test_getHostrange_replicon():
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

    for matchtype in ["exact","multiexact","loose"]:
        (returnhrrank,returnhrsciname,unique_ref_selected_taxids,ref_taxids_df,stats_host_range_dict) = getRefSeqHostRange(
                           replicon_name_list = ["IncI1"],
                           mob_cluster_id= None,
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
    with pytest.raises(Exception) as e:
        getRefSeqHostRange(
            replicon_name_list=None,
            mob_cluster_id=None,
            relaxase_name_acc_list=None,
            relaxase_name_list=None,
            matchtype="exact",
            hr_obs_data=loadHostRangeDB())
    assert str(e.value) == "Empty dataframe returned from RefSeq plasmid database", "Error message: "+str(e.value)

    #test when there are hits for RefSeq host-range search only (e.g. ColRNAI_rep_cluster_1857)
    getRefSeqHostRange(
        replicon_name_list=["ColRNAI_rep_cluster_1857"],
        mob_cluster_id=None,
        relaxase_name_acc_list=None,
        relaxase_name_list=None,
        matchtype="exact",
        hr_obs_data=loadHostRangeDB())

def test_getHostrange_clusterID():
    (convergance_rank, converged_taxonomy_name,
     unique_ref_selected_taxids, ref_taxids_df,
     stats_host_range_dict) = getRefSeqHostRange(replicon_name_list=None,
                       mob_cluster_id=476,  #IncI1 plasmid cluster
                       relaxase_name_acc_list=None,
                       relaxase_name_list=None,
                       matchtype="multiexact",
                       hr_obs_data=loadHostRangeDB())
    assert convergance_rank == "family", "The host range rank is incorrect. Reported "+ convergance_rank
    assert converged_taxonomy_name == "Enterobacteriaceae", "The host range rank scientific name is incorrect. Reported " + converged_taxonomy_name
    assert unique_ref_selected_taxids != None, "Empty taxonomy id list."+unique_ref_selected_taxids+" Check search criteria"

def test_getHostrange_relaxase_acc():
    (convergance_rank, converged_taxonomy_name,
     unique_ref_selected_taxids, ref_taxids_df,
     stats_host_range_dict) = getRefSeqHostRange(replicon_name_list=None,
                       mob_cluster_id=None,
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
                        mob_cluster_id = None,
                        relaxase_name_acc_list = ["NC_009739_00028","NC_016978_00054"],
                        relaxase_name_list = None,
                        matchtype = "multiexact",
                        hr_obs_data = loadHostRangeDB())
    assert convergance_rank == "phylum", "The host range rank is incorrect. Reported " + convergance_rank
    assert converged_taxonomy_name == "Proteobacteria", "The host range rank scientific name is incorrect. Reported " + converged_taxonomy_name
    assert unique_ref_selected_taxids != None, "Empty taxonomy id list." + unique_ref_selected_taxids + " Check search criteria"


def test_getHostrange_relaxase_name():
    (convergance_rank, converged_taxonomy_name,
     unique_ref_selected_taxids, ref_taxids_df,
     stats_host_range_dict) = getRefSeqHostRange(replicon_name_list=None,
                                                 mob_cluster_id=None,
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
                                                 mob_cluster_id=None,
                                                 relaxase_name_acc_list=None,
                                                 relaxase_name_list=["MOBV","MOBC"],
                                                 matchtype="multiexact",
                                                 hr_obs_data=loadHostRangeDB())
    assert convergance_rank == "superkingdom", "The host range rank is incorrect. Reported " + convergance_rank
    assert converged_taxonomy_name == "Bacteria", "The host range rank scientific name is incorrect. Reported " + converged_taxonomy_name
    assert unique_ref_selected_taxids != None, "Empty taxonomy id list." + unique_ref_selected_taxids + " Check search criteria"

def test_literature_hostrange_replion():
    """
    replicon_names,plasmid_lit_db,input_seq=""
    Test the literature-based plasmid host range prediction and transfer rate
    """
    lit_report_columns = ["LiteratureQueryReplicon",
                          "LiteratureSearchReplicon",
                          "LiteratureFoundPlasmidsNames",
                          "LiteratureFoundPlasmidsNumber",
                          "LiteratureReportedHostRangePlasmidClass",
                          "LiteratureReportedHostPlasmidSpecies",
                          "LiteratureReportedPlasmidHostSpeciesNumber",
                          "LiteraturePredictedDBHostRangeTreeRank",
                          "LiteraturePredictedDBHostRangeTreeRankSciName",
                          "LiteratureReportedHostRangeInPubs",
                          "LiteratureMinTransferRateRange",
                          "LiteratureMaxTransferRateRange",
                          "LiteratureMeanTransferRateRange",
                          "LiteraturePMIDs", "LiteraturePublicationsNumber"]
    report,littaxids = getLiteratureBasedHostRange(replicon_names =["IncF"],
                                plasmid_lit_db = loadliteratureplasmidDB(),input_seq="")
    assert report.empty == False, "Literature host range prediction is empty. Check your search criteria/criterium"
    assert report.columns.tolist() == lit_report_columns, "Literature host range report does not match all default columns"
    assert report.shape[0] == 1, "Literature host range report dimension is incorrect. The expect dimension is a single row (1,15)"
    assert report.loc[0,"LiteraturePredictedDBHostRangeTreeRank"] == "family", "Literature hit-based host range host range prediction is incorrect. Expected family, got"+report.loc[0,"LiteraturePredictedDBHostRangeTreeRank"]
    assert report.loc[0,"LiteraturePredictedDBHostRangeTreeRankSciName"] == "Enterobacteriaceae"
    assert report.loc[0,"LiteratureReportedHostRangeInPubs"] == "family","Literature reported host range host range prediction is incorrect. Expected family, got"+report.loc[0,"LiteraturePredictedDBHostRangeTreeRank"]

    #test when there are no hits for a given replicon in literature database
    report, littaxids = getLiteratureBasedHostRange(replicon_names=["ColRNAI_rep_cluster_1857"],
                                                    plasmid_lit_db=loadliteratureplasmidDB(), input_seq="")
    assert report.empty, "The literature results df should be empty but it is not"



def test_getTaxonomyTree_onRefSeqDB_and_Render():
    convergance_rank, converged_taxonomy_name, unique_ref_selected_taxids, \
    ref_taxids_df, stats_host_range_dict = getRefSeqHostRange(replicon_name_list=["IncP"],matchtype="exact",
                                               mob_cluster_id=None, relaxase_name_acc_list=None,
                                               relaxase_name_list= None, hr_obs_data=loadHostRangeDB())
    tree = getTaxonomyTree(taxids=unique_ref_selected_taxids)
    assert isinstance(tree, ete3.PhyloTree), "Class mismatch output by getTaxonomyTree(). Output class is " + type(tree).__name__


    renderTree(tree=tree, taxids=unique_ref_selected_taxids, filename_prefix="./run_test/run_test_refseqhostrange")

def test_getTaxonomyTree_onLiteratureDB_and_Render():
    report_df, lit_taxids = getLiteratureBasedHostRange(replicon_names=["IncP"],plasmid_lit_db=loadliteratureplasmidDB(),)
    tree = getTaxonomyTree(taxids=lit_taxids)
    assert isinstance(tree,ete3.PhyloTree), "Class mismatch output by getTaxonomyTree(). Output class is "+type(tree).__name__
    renderTree(tree=tree, taxids=lit_taxids, filename_prefix="./run_test/run_test_literaturehostrange")




#test the entire mob-typer + mob_host_range modules. AB040415 has multiple replicons (IncFIB,IncFII)
def test_mob_typer_host_range_multi_replicon():
    args = [
        "--infile", os.path.dirname(__file__) + "/TestData/AB040415.fasta",
        "--outdir", "run_test",
        "--host_range_detailed"
    ]
    sys.argv[1:] = args
    mob_suite.mob_typer.main()
    results_df = pandas.read_csv("./run_test/mobtyper_AB040415_report.txt", sep="\t")
    assert results_df["RefSeqHRrank"].values[0] == "order"
    assert results_df["RefSeqHRSciName"].values[0] == "Enterobacterales"

    args = [
        "--infile", os.path.dirname(__file__) + "/TestData/AY603981.fasta",
        "--outdir", "run_test",
        "--host_range"
    ]
    sys.argv[1:] = args
    mob_suite.mob_typer.main()
    results_df = pandas.read_csv("./run_test/mobtyper_AY603981_report.txt", sep="\t")
    assert results_df["RefSeqHRrank"].values[0] == "-"
    assert results_df["RefSeqHRSciName"].values[0] == "-"
    assert results_df["PredictedMobility"].values[0] == "Mobilizable"

    args = [
        "--infile", os.path.dirname(__file__) + "/TestData/AB011548.fasta",
        "--outdir", "run_test",
        "--host_range"
    ]
    sys.argv[1:] = args
    mob_suite.mob_typer.main()
    results_df = pandas.read_csv("./run_test/mobtyper_AY603981_report.txt", sep="\t")
    assert results_df["RefSeqHRrank"].values[0] == "-"
    assert results_df["RefSeqHRSciName"].values[0] == "-"
    assert results_df["PredictedMobility"].values[0] == "Mobilizable"

    args = [
        "--infile", os.path.dirname(__file__) + "/TestData/AB011548.fasta",
        "--outdir", "run_test",
        "--host_range"
    ]
    sys.argv[1:] = args
    mob_suite.mob_typer.main()
    results_df = pandas.read_csv("./run_test/mobtyper_AB011548_report.txt", sep="\t")
    assert results_df["RefSeqHRrank"].values[0] == "superkingdom"
    assert results_df["RefSeqHRSciName"].values[0] == "Bacteria"
    assert results_df["PredictedMobility"].values[0] == "Mobilizable"
    assert results_df["LitPredDBHRRank"].values[0]  == "-"
    assert results_df["LitPredDBHRRankSciName"].values[0] == "-"



#    assert any([len(re.findall("order\tEnterobacterales",out)) >= 1 for out in output]) == True, "Something went wrong with host range prediction";

# mob_typer -i /Users/kirill/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/tests/TestData/IncF/ET11_Ecoli_plasmid_529.fasta -o run_test --host_range_detailed
# mob_typer -i /Users/kirill/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/tests/TestData/IncF/ET5_Ecoli_plasmid_973.fasta -o run_test --host_range
# mob_typer -i /Users/kirill/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/tests/TestData/IncF/ET4_Ecoli_plasmid_969.fasta -o run_test --host_range
# mob_typer -i /Users/kirill/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/tests/TestData/AB040415.fasta -o run_test --host_range
# mob_typer -i /Users/kirill/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/tests/TestData/AY603981.fasta -o run_test --host_range #no replicons
# mob_typer -i /Users/kirill/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/tests/TestData/AB011548.fasta -o run_test --host_range_detailed