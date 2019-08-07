import os,sys, logging,  subprocess,time

logger=logging.getLogger()
LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
logging.basicConfig(format=LOG_FORMAT, level=logging.DEBUG)

#"repetitive.dna.fas.nin"
def check_file_hash(default_database_dir):
    file_ref_list = ["repetitive.dna.fas.nin",
                       "ncbi_plasmid_full_seqs.fas.msh",
                       "repetitive.dna.fas",
                       "mob.proteins.faa",
                       "rep.dna.fas",
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
        file_path=default_database_dir+"/"+file
        assert os.path.exists(file_path), "File "+file+"does not exist. Check downloaded database and retry."
        cleanup(file_path)

def cleanup(full_path_file):
    os.remove(full_path_file)

def test_download_databases_with_input_dir():

    if os.path.exists("run_test") == False:
        os.mkdir("run_test")

    database_dir = "run_test/databases"
    args = [
        "-d", database_dir,
    ]
    sys.argv[1:] = args

def test_concurrent_init():


    p1 = subprocess.Popen(['python', '../mob_init.py', '-v', '-d' 'run_test/databases'],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           shell=False
                           )
    print("Strated process 1  ...")
    time.sleep(10)

    p2 = subprocess.Popen(['python', '../mob_init.py', '-v', '-d' 'run_test/databases'],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           shell=False)
    print("Strated process 2  ...")

    p1stdout, p1stderr = p1.communicate()
    p2stdout, p2stderr = p2.communicate()
    print("Process 1", p1stdout.decode())
    print("Process 2", p2stdout.decode())
    assert "completed successfully" in p1stdout.decode()
    assert "completed successfully" in p2stdout.decode()





