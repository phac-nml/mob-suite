import mob_suite.mob_typer
import os,sys
import pandas


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
# cp  *.py /Users/kirill/miniconda/envs/mob_suite_test/lib/python3.6/site-packages/mob_suite/