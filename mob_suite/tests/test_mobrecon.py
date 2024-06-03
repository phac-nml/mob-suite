import sys, os, logging, json, pandas, pytest, random, hashlib
import mob_suite.mob_recon
import mob_suite.utils
from mob_suite.constants import MOB_CLUSTER_INFO_HEADER, default_database_dir
from Bio import SeqIO

#test all mob_recon functions including aggregation of results
TEST_ROOT = os.path.dirname(__file__)
ROOT_MODULE = os.path.dirname(mob_suite.__file__)
logger=logging.getLogger()
LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
logging.basicConfig(format=LOG_FORMAT, level=logging.DEBUG)

def check_if_output_dir_exists_and_create():
    if not os.path.exists(os.path.join(TEST_ROOT,"run_test")):
        os.mkdir(os.path.join(TEST_ROOT,"run_test"))


def test_mob_recon_with_mob_typer_report():
    check_if_output_dir_exists_and_create()
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
    check_if_output_dir_exists_and_create()
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
    check_if_output_dir_exists_and_create()
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


def test_contig_order_on_final_results(tmpdir):
    check_if_output_dir_exists_and_create()
    in_fasta = ROOT_MODULE + "/example/SRR3703080_illumina_unicycler.fasta"
    records=SeqIO.index(in_fasta, "fasta")
    for iter in range(0,2):
        contigs_list = [r for r in records.keys()]
        random.shuffle(contigs_list)
        with open(os.path.join(tmpdir, f"shuffled_{iter}.fasta"), "w") as output_handle:
            SeqIO.write([records[k] for k in contigs_list], output_handle, "fasta")
        args = [
            "--infile", os.path.join(tmpdir, f"shuffled_{iter}.fasta"),
            "--outdir", TEST_ROOT+f"/run_test/mob_recon_SRR3703080_shuffled_{iter}",
            "--debug",
            "--force"
        ]
        sys.argv[1:] = args
        mob_suite.mob_recon.main()

    df1 = pandas.read_csv(TEST_ROOT+f"/run_test/mob_recon_SRR3703080_shuffled_0/contig_report.txt", sep="\t")[1:] #remove sample_id column
    df2 = pandas.read_csv(TEST_ROOT+f"/run_test/mob_recon_SRR3703080_shuffled_1/contig_report.txt", sep="\t")[1:] #remove sample_id column
    assert all(df1==df2), "Two contig_reports.txt are not equal"

    df1 = pandas.read_csv(TEST_ROOT+f"/run_test/mob_recon_SRR3703080_shuffled_0/mobtyper_results.txt", sep="\t")[1:]
    df2 = pandas.read_csv(TEST_ROOT+f"/run_test/mob_recon_SRR3703080_shuffled_1/mobtyper_results.txt", sep="\t")[1:]
    assert all(df1==df2), "Two mobtyper_results.txt are not equal"

    df1 = pandas.read_csv(TEST_ROOT+f"/run_test/mob_recon_SRR3703080_shuffled_0/biomarkers.blast.txt", sep="\t")[1:]
    df2 = pandas.read_csv(TEST_ROOT+f"/run_test/mob_recon_SRR3703080_shuffled_1/biomarkers.blast.txt", sep="\t")[1:]
    assert all(df1==df2), "Two biomarkers.blast.txt are not equal"

    hash_1 = hashlib.md5(open(TEST_ROOT+f"/run_test/mob_recon_SRR3703080_shuffled_0/plasmid_AA474.fasta",'rb').read()).hexdigest()
    hash_2 = hashlib.md5(open(TEST_ROOT+f"/run_test/mob_recon_SRR3703080_shuffled_1/plasmid_AA474.fasta",'rb').read()).hexdigest()
    assert hash_1==hash_2, "Two plasmid_AA474.fasta files are not equal"