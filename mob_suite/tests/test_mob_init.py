import os, sys, logging, subprocess, time, mob_suite.mob_init

TEST_ROOT = os.path.dirname(__file__)

logger = logging.getLogger()
LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
logging.basicConfig(format=LOG_FORMAT, level=logging.DEBUG)

try:
    os.mkdir(TEST_ROOT + "/run_test")
except:
    pass


# "repetitive.dna.fas.nin"
def test_download_databases_with_input_dir():
    if os.path.exists(TEST_ROOT + "/run_test") == False:
        os.mkdir(TEST_ROOT + "/run_test")

    database_dir = TEST_ROOT + "/run_test/databases_custom_dir"
    args = [
        "-d", database_dir,
    ]
    sys.argv[1:] = args
    mob_suite.mob_init.main()


def test_concurrent_init():
    try:
        p1 = subprocess.Popen(['mob_init', '-v', '-d', TEST_ROOT + '/run_test/databases'],
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              shell=False
                              )
    except:
        p1 = subprocess.Popen([sys.executable, '../mob_init.py', '-v', '-d', TEST_ROOT + '/run_test/databases'],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              )
    print("Strated process 1  ...")
    time.sleep(10)

    try:
        p2 = subprocess.Popen(['mob_init', '-v', '-d', TEST_ROOT + '/run_test/databases'],
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              shell=False)
    except:
        p2 = subprocess.Popen([sys.executable, '../mob_init.py', '-v', '-d', TEST_ROOT + '/run_test/databases'],
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              shell=False)
    print("Strated process 2  ...")
    time.sleep(100);
    p2.kill()

    p1stdout, p1stderr = p1.communicate()
    p2stdout, p2stderr = p2.communicate()
    print("Process 1 stdout steam", p1stdout.decode(), "Process 1 stderr steam", p1stderr.decode())
    print("Process 2 stdout steam", p2stdout.decode(), "Process 2 stderr steam", p2stderr.decode())
    assert "MOB init completed successfully" in p1stderr.decode()
    assert "Lock file found" in p2stderr.decode()


def check_file_list(default_database_dir):
    file_ref_list = ["repetitive.dna.fas.nin",
                     "ncbi_plasmid_full_seqs.fas.msh",
                     "repetitive.dna.fas",
                     "mob.proteins.faa",
                     "rep.dna.fas",
                     "clusters.txt",
                     "repetitive.dna.fas.nsq",
                     "ncbi_plasmid_full_seqs.fas.nhr",
                     "ncbi_plasmid_full_seqs.fas.nin",
                     "mpf.proteins.faa",
                     "orit.fas",
                     "ncbi_plasmid_full_seqs.fas.nsq",
                     "repetitive.dna.fas.nhr",
                     "ncbi_plasmid_full_seqs.fas"
                     ]

    for file in file_ref_list:
        file_path = default_database_dir + "/" + file
        assert os.path.exists(file_path), "File " + file + "does not exist. Check downloaded database and retry."
