import sys, os, logging
import mob_suite.mob_recon
#test all mob_recon functions including aggregation of results
TEST_ROOT = os.path.dirname(__file__)
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
