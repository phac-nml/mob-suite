import mob_suite.mob_typer
import os,sys
import pandas
import logging

TEST_ROOT = os.path.dirname(__file__)
LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
logging.basicConfig(format=LOG_FORMAT)
logger=logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def test_mob_typer_IncX1():
    create_output_dir()
    args = [
        "--infile", TEST_ROOT + "/TestData/IncX/IncX1.fasta",
        "--outdir", TEST_ROOT+"/run_test",
        "--host_range_detailed",
        "--debug"
    ]
    sys.argv[1:] = args
    mob_suite.mob_typer.main()
    results_df = pandas.read_csv(os.path.join(TEST_ROOT,"run_test/mobtyper_IncX1.fasta_report.txt"), sep="\t")
    print(results_df)
    assert results_df["LitPMIDs"].values[0] == "21625636;22470007"

def create_output_dir():
    if not os.path.exists(os.path.join(TEST_ROOT, "run_test")):
        os.mkdir(os.path.join(TEST_ROOT, "run_test"))
#test the entire mob-typer + mob_host_range modules. AB040415 has multiple replicons (IncFIB,IncFII)
#IncFIB,IncFII multi-plasmids
def test_mob_typer_host_range_multi_replicon():
    create_output_dir()

    logger.info("Testing mob_typer on IncF {} plasmid".format("AB040415.fasta"))
    logger.info("Current working directory:{}".format(os.getcwd()))
    logger.info("List diretory of the input files: {}".format(os.listdir(os.path.dirname(__file__)+"/TestData/")))

    args = [
        "--infile", os.path.dirname(__file__) + "/TestData/AB040415.fasta",
        "--outdir", TEST_ROOT+"/run_test",
        "--host_range_detailed"
    ]
    sys.argv[1:] = args
    mob_suite.mob_typer.main()
    results_df = pandas.read_csv(os.path.join(TEST_ROOT,"run_test/mobtyper_AB040415.fasta_report.txt"), sep="\t")

    assert results_df["NCBI-HR-rank"].values[0] == "order"
    assert results_df["NCBI-HR-Name"].values[0] == "Enterobacterales"
    assert results_df["LitRepHRPlasmClass"].values[0]  == "NarrowHostRange"
    assert results_df["LitPredDBHRRank"].values[0]  == "family"
    assert results_df["LitPredDBHRRankSciName"].values[0]  == "Enterobacteriaceae"
    assert results_df["LitRepHRRankInPubs"].values[0]  == "family"
    assert results_df["LitRepHRNameInPubs"].values[0]  == "Enterobacteriaceae"
    assert results_df["LitPMIDsNumber"].values[0]  == 5


    args = [
        "--infile", os.path.dirname(__file__) + "/TestData/AB011548.fasta",
        "--outdir", TEST_ROOT+"/run_test",
        "--host_range"
    ]
    sys.argv[1:] = args
    mob_suite.mob_typer.main()

    results_df = pandas.read_csv(os.path.join(TEST_ROOT,"run_test/mobtyper_AB011548.fasta_report.txt"), sep="\t")
    assert results_df["NCBI-HR-rank"].values[0] == "superkingdom"
    assert results_df["NCBI-HR-Name"].values[0] == "Bacteria"
    assert results_df["PredictedMobility"].values[0] == "Mobilizable"
    assert results_df["LitPredDBHRRank"].values[0]  == "-"
    assert results_df["LitPredDBHRRankSciName"].values[0] == "-"

def test_mob_typer_host_range_multi_replicon_KU295134():
    create_output_dir()
    #KU295134.fasta with IncFII and IncN replicons with closest literature reference NC_011385
    #suitable to check the multi-replicon case and how host range data collapse is functioning
    logger.info("Testing mob_typer on IncFII and IncF {} plasmid".format("KU295134.fasta"))
    logger.info("Current working directory:{}".format(os.getcwd()))
    logger.info("List diretory of the input files: {}".format(os.listdir(os.path.dirname(__file__) + "/TestData/")))

    args = [
        "--infile", os.path.dirname(__file__) + "/TestData/KU295134.fasta",
        "--outdir", TEST_ROOT+"/run_test",
        "--host_range_detailed"
    ]
    sys.argv[1:] = args
    mob_suite.mob_typer.main()
    results_df = pandas.read_csv(os.path.join(TEST_ROOT,"run_test/mobtyper_KU295134.fasta_report.txt"), sep="\t")


    assert results_df["NCBI-HR-rank"].values[0] == "class"
    assert results_df["NCBI-HR-Name"].values[0] == "Gammaproteobacteria"
    assert results_df["PredictedMobility"].values[0] == "Conjugative"

def test_mob_typer_host_range_no_replicon_data():
    create_output_dir()
    logger.info("Testing mob_typer on Mobilizable {} plasmid from Pseudomonas".format("AY603981.fasta"))
    logger.info("Current working directory:{}".format(os.getcwd()))
    logger.info("List diretory of the input files: {}".format(os.listdir(os.path.dirname(__file__) + "/TestData/")))

    args = [
        "--infile", os.path.dirname(__file__) + "/TestData/AY603981.fasta",
        "--outdir", TEST_ROOT+"/run_test",
        "--host_range_detailed"
    ]
    sys.argv[1:] = args
    mob_suite.mob_typer.main()
    results_df = pandas.read_csv(os.path.join(TEST_ROOT,"run_test/mobtyper_AY603981.fasta_report.txt"), sep="\t")

    assert results_df["NCBI-HR-rank"].values[0] == "genus"
    assert results_df["NCBI-HR-Name"].values[0] == "Pseudomonas"
    assert results_df["PredictedMobility"].values[0] == "Mobilizable"

def test_mob_typer_broad_host_range_IncF():
    create_output_dir()
    logger.info("Testing mob_typer on IncF {} plasmid".format("ET4_Ecoli_plasmid_969.fasta"))
    logger.info("Current working directory:{}".format(os.getcwd()))
    logger.info("List diretory of the input files: {}".format(os.listdir(os.path.dirname(__file__) + "/TestData/")))

    args = [
        "--infile", os.path.dirname(__file__) + "/TestData/IncF/ET4_Ecoli_plasmid_969.fasta",
        "--outdir", TEST_ROOT+"/run_test",
        "--host_range_detailed"
    ]
    sys.argv[1:] = args
    mob_suite.mob_typer.main()
    results_df = pandas.read_csv(os.path.join(TEST_ROOT,"run_test/mobtyper_ET4_Ecoli_plasmid_969.fasta_report.txt"), sep="\t")

    assert results_df["NCBI-HR-rank"].values[0] == "order"
    assert results_df["NCBI-HR-Name"].values[0] == "Enterobacterales"
    assert results_df["PredictedMobility"].values[0] == "Conjugative"


def test_mean_and_multireplicon_frame():
    logger.info("Testing multireplicon case with multiple transfer rate means")
    args = [
        "--infile", os.path.dirname(__file__) + "/TestData/pCAV1453-208.fasta",
        "--outdir", TEST_ROOT + "/run_test",
        "--host_range_detailed"
    ]

    sys.argv[1:] = args
    mob_suite.mob_typer.main()