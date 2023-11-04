import sys, os, logging
import mob_suite.mob_recon
#test all mob_recon functions including aggregation of results
TEST_ROOT = os.path.dirname(__file__)
ROOT_MODULE = os.path.dirname(mob_suite.__file__)
logger=logging.getLogger()
LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
logging.basicConfig(format=LOG_FORMAT, level=logging.DEBUG)

if not os.path.exists(os.path.join(TEST_ROOT,"run_test")):
    os.mkdir(os.path.join(TEST_ROOT,"run_test"))


def test_mob_recon_with_mob_typer_report():
    if os.path.exists("run_test") == False:
        os.mkdir("run_test")
    #IncFIB,IncFII multi-plasmids
    args = [
        "--infile", os.path.dirname(__file__) + "/TestData/Pseudomonas/test_contigs.fasta",
        "--outdir", os.path.dirname(__file__)+"/run_test/mob_recon",
        "--debug",
        "--force"
    ]
    sys.argv[1:] = args
    mob_suite.mob_recon.main()

    mobtyper_results_file = os.path.join(TEST_ROOT,"run_test/mob_recon/mobtyper_results.txt")
    assert sum(1 for line in open(mobtyper_results_file)) == 2 , "Results file is empty, something went wrong"

def test_mob_recon_no_plasmid_biomarkers():

    args = [
        "--infile", ROOT_MODULE + "/example/assembly_no_biomarkers.fasta",
        "--outdir", TEST_ROOT + "/run_test/mob_recon_no_plasmid_markers",
        "--debug",
        "--force"
    ]

    sys.argv[1:] = args
    mob_suite.mob_recon.main()
    assert os.path.exists(os.path.join(TEST_ROOT,"run_test/mob_recon_no_plasmid_markers/biomarker_report.txt")) == False

def test_mob_recon_typical_run():
    args = [
        "--infile", ROOT_MODULE + "/example/SRR3703080_illumina_unicycler.fasta",
        "--outdir", TEST_ROOT+"/run_test/mob_recon_SRR3703080",
        "--debug",
        "--force"
    ]
    sys.argv[1:] = args
    mob_suite.mob_recon.main()

    assert os.path.exists(os.path.join(TEST_ROOT,"run_test/mob_recon_SRR3703080/biomarkers.blast.txt")) == True, "File does not exist"
    assert os.path.exists(os.path.join(TEST_ROOT,"run_test/mob_recon_SRR3703080/chromosome.fasta")) == True, "File does not exist"
    assert os.path.exists(os.path.join(TEST_ROOT,"run_test/mob_recon_SRR3703080/mobtyper_results.txt")) == True, "File does not exist"
    assert os.path.exists(os.path.join(TEST_ROOT,"run_test/mob_recon_SRR3703080/contig_report.txt")) == True, "File does not exits"
    assert os.path.exists(os.path.join(TEST_ROOT,"run_test/mob_recon_SRR3703080/plasmid_AA474.fasta")) == True, "File does not exists"

    