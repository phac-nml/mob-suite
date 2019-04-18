import sys, os,re,pandas
from mob_suite.mob_host_range import getRefSeqHostRange, getLiteratureBasedHostRange, loadliteratureplasmidDB, loadHostRangeDB
import mob_suite.mob_typer


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
        assert ref_taxids_df.columns.tolist() == ref_taxids_df_columns, "RefSeq plasmid database column names do not match"
        assert stats_host_range_dict != {}, "Empty RefSeq host-range hits dictionary. Check the getRefSeqHostRange()"


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
    report = getLiteratureBasedHostRange(replicon_names =["IncF"],
                                plasmid_lit_db = loadliteratureplasmidDB(),
                                input_seq="")
    assert report.empty == False, "Literature host range prediction is empty. Check your search criteria/criterium"
    assert report.columns.tolist() == lit_report_columns, "Literature host range report does not match all default columns"
    assert report.shape[0] == 1, "Literature host range report dimension is incorrect. The expect dimension is a single row (1,15)"
    assert report.loc[0,"LiteraturePredictedDBHostRangeTreeRank"] == "family", "Literature hit-based host range host range prediction is incorrect. Expected family, got"+report.loc[0,"LiteraturePredictedDBHostRangeTreeRank"]
    assert report.loc[0,"LiteraturePredictedDBHostRangeTreeRankSciName"] == "Enterobacteriaceae"
    assert report.loc[0,"LiteratureReportedHostRangeInPubs"] == "family","Literature reported host range host range prediction is incorrect. Expected family, got"+report.loc[0,"LiteraturePredictedDBHostRangeTreeRank"]




def test_getHostrange_clusterID(capsys):
    getRefSeqHostRange(None,267,None,None,"multi")
    out, err = capsys.readouterr()
    assert "RANGE:Enterobacteriaceae" in out , "host range prediction was incorrect"


#test the entire mob-typer + mob_host_range modules. AB040415 has multiple replicons (IncFIB,IncFII)
def test_mob_typer_host_range_multi_replicon():
   args=[
          "--infile", os.path.dirname(__file__)+"/TestData/AB040415.fasta",
          "--outdir", "test_AB040415",
          "--host_range"
          ]

   sys.argv[1:] = args
   mob_suite.mob_typer.main()

   with open(file=os.path.dirname(__file__)+"/test_AB040415/mobtyper_AB040415.fasta_report.txt") as fp_in:
       output = fp_in.readlines()

   assert any([len(re.findall("order\tEnterobacterales",out)) >= 1 for out in output]) == True, "Something went wrong with host range prediction";

# mob_typer -i /Users/kirill/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/tests/TestData/IncF/ET11_Ecoli_plasmid_529.fasta -o run_test --host_range
# mob_typer -i /Users/kirill/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/tests/TestData/IncF/ET5_Ecoli_plasmid_973.fasta -o run_test --host_range
# mob_typer -i /Users/kirill/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/tests/TestData/IncF/ET4_Ecoli_plasmid_969.fasta -o run_test --host_range