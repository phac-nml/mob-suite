import mob_suite.mob_init
import os,sys, hashlib,shutil

#"repetitive.dna.fas.nin"
def check_file_hashesh(default_database_dir):
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
    default_database_dir = os.path.dirname(__file__) + "/run_test/databases"
    args = [
        "-d", default_database_dir,
    ]
    sys.argv[1:] = args

    mob_suite.mob_init.main()
    check_file_hashesh(default_database_dir)


def test_download_databases_with_custom_url():
    default_database_dir = os.path.dirname(__file__) + "/run_test/databases"

    test_db_urls = [    "https://ndownloader.figshare.com/articles/5841882/versions/1",
                        "https://share.corefacility.ca/index.php/s/oeufkw5HyKz0X5I/download"]
    for test_db_url in test_db_urls:

        args = [
            "-d", default_database_dir,
            "-u", test_db_url
        ]
        sys.argv[1:] = args

        print("Testing "+test_db_url+"...")
        mob_suite.mob_init.main()
        check_file_hashesh(default_database_dir)

