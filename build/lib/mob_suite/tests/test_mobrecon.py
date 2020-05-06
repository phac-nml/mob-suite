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
        "--run_typer",
        "--outdir", os.path.dirname(__file__)+"/run_test/mob_recon",
        "--debug"
    ]
    sys.argv[1:] = args
    mob_suite.mob_recon.main()

    mobtyper_results_file = os.path.join(TEST_ROOT,"run_test/mob_recon/mobtyper_aggregate_report.txt")
    assert sum(1 for line in open(mobtyper_results_file)) == 4 , "Results file is empty, something went wrong"

def test_run_mob_typer():
    """
    Test if mob_recon can call mob_typer and successfuly concatenate mob_typer results
    :return:
    """
    plasmid_files=[TEST_ROOT+"/run_test/mob_recon/plasmid_novel_0.fasta",
                   TEST_ROOT+"/run_test/mob_recon/plasmid_novel_1.fasta"]
    out_dir=TEST_ROOT+"/run_test/mob_recon"
    num_threads=1
    #database_dir="/Users/kirill/WORK/MOBSuiteHostRange2018/Source/mob-suite/mob_suite/databases"

    mobtyper_results="Tab-file header will go here. Line 1\n"
    for file in plasmid_files:
        file=os.path.join(TEST_ROOT,file)
        mobtyper_results = mobtyper_results + "{}\n".format(mob_suite.mob_recon.run_mob_typer(plasmid_file_abs_path=file,
                                                                                            outdir=out_dir,
                                                                                            num_threads=int(num_threads)
                                                                                            ))
    mobtyper_results_file = os.path.join(out_dir, 'mobtyper_aggregate_report.txt')
    fh = open(mobtyper_results_file, 'w')
    fh.write(mobtyper_results)
    print("mob_typer aggregate results",mobtyper_results)
    fh.close()
    assert sum(1 for line in open(mobtyper_results_file)) == 3, "Results file is empty, something went wrong"
