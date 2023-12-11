import os, time, random, sys
from statistics import mean
from Bio import SeqIO
import pandas as pd
import logging
import itertools
import mob_suite.mob_typer



SAMPLES_WITH_ISSUES_DF = pd.DataFrame(columns=['plasmid_id','issue_type','marker_type','marker_name','marker_accession'])
TEST_ROOT = os.path.dirname(__file__)
PACKAGE_ROOT = os.path.dirname(mob_suite.__file__)
LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
NRANDOMSAMPLES2TEST = 300
COMPLEX_CASES_LIST = ['CP011291', 'CP045773', 'KU932024', 'NC_002134', 'CP033122', 'CP031575', 'CP021102',
                      'CP033949', 'CP025336', 'CP025336', 'CP028546']

logging.basicConfig(format=LOG_FORMAT)
logger=logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def clean_missing_dashes(L):
    return [i for i in L if i != '-']


def test_mean_and_multireplicon_frame(mode = None):
    if os.path.exists("run_test") == False:
        os.mkdir("run_test")

    if mode == None:
        return 0
    logger.info("Testing MOB-typer for biomarker and mobtyper reports consistency in terms of biomarker names (not accessions)")
    logger.info(f"Selected testing mode: '{mode}'")
    elapsed_time_list=list()
    args = [
        "--num_threads", "4",
        "--infile", TEST_ROOT + "/run_test/tmp_plasmid.fasta",
        "--out_file", TEST_ROOT + "/run_test/mobtyper_report.txt",
        "--biomarker_report_file", TEST_ROOT + "/run_test/biomarker_report.txt"
    ]
    sys.argv[1:] = args
    
    with open(os.path.join(PACKAGE_ROOT + "/databases/ncbi_plasmid_full_seqs.fas")) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        if mode == "random sampling":
            records=random.sample(records, NRANDOMSAMPLES2TEST)
        n_records = len([rec.id for rec in records])
        for n, record in enumerate(records):
            if mode == "list of sample ids":
                if record.id not in COMPLEX_CASES_LIST:
                    continue
            start_time = time.time()
            logger.info(f"*** Sample #{n+1} out of {n_records} -> record ID {record.id} ***")
            SeqIO.write(record, TEST_ROOT + "/run_test/tmp_plasmid.fasta", "fasta")
            mob_suite.mob_typer.main()
            mobtyper_report_df = pd.read_csv(TEST_ROOT + "/run_test/mobtyper_report.txt",sep="\t")
            if os.path.getsize( TEST_ROOT + "/run_test/biomarker_report.txt") == 1:
                logger.warning(f'Empty biomarker report for {record.id}')
                continue

            # check biomarker report in sync with mobtyper report    
            mobtyper_biomarker_df = pd.read_csv(TEST_ROOT + "/run_test/biomarker_report.txt",sep="\t")   
            mobtyper_marker_types = clean_missing_dashes(list(itertools.chain(* mobtyper_report_df[['orit_type(s)','rep_type(s)', 'relaxase_type(s)', 'mpf_type']].stack().str.split(','))))
            mobtyper_marker_accessions=clean_missing_dashes(list(itertools.chain(* mobtyper_report_df[['rep_type_accession(s)','relaxase_type_accession(s)', 'mpf_type_accession(s)', 'orit_accession(s)']].stack().str.split(','))))
            
            #for accession, biomarker_name in mobtyper_biomarker_df.qseqid.str.split('|'):
            for index, row in mobtyper_biomarker_df.iterrows():
                accession, biomarker_name = row.qseqid.split('|')       
                match_result = mobtyper_report_df.apply(lambda item: item.astype(str).str.contains(accession) , axis=0).any().any(axis=None)
                if match_result == False:
                    logger.info(f"Issue with {record.id}. Marker with accession {biomarker_name}:{accession} is not found in mobtyper report")
                    SAMPLES_WITH_ISSUES_DF.loc[SAMPLES_WITH_ISSUES_DF.shape[0],:]=[
                        record.id,'missing in mobtyper report',row.biomarker, biomarker_name, accession]
                if biomarker_name in mobtyper_marker_types:
                    mobtyper_marker_types.remove(biomarker_name)
                if accession in mobtyper_marker_accessions:
                    mobtyper_marker_accessions.remove(accession)   
            if len(mobtyper_marker_types) != 0:
                logger.error(f'Following biomarker types are missing in biomarkers report:{",".join(mobtyper_marker_types)}')
                SAMPLES_WITH_ISSUES_DF.loc[SAMPLES_WITH_ISSUES_DF.shape[0],:]=[
                    record.id,'missing in biomarkers report',row.biomarker,",".join(mobtyper_marker_types),",".join(mobtyper_marker_accessions)]
        
            elapsed_time=time.time() - start_time
            elapsed_time_list.append(elapsed_time)
            logger.info("Elapsed {0:.0f} seconds for {1} sample typing".format(elapsed_time, record.id))
            logger.info("*** Estimated time to completion {0:.2f} hours (mean time per test {1:.0f} s) ***".format((n_records-n)*mean(elapsed_time_list)/3600, mean(elapsed_time_list) ))
           

    SAMPLES_WITH_ISSUES_DF.to_csv(TEST_ROOT+'/plasmid_accessions_with_issues.txt', sep="\t", index=False)          
    n_issues2solve = len(SAMPLES_WITH_ISSUES_DF['plasmid_id'].unique()) 
    logger.info(f"Analysis completed with {n_issues2solve} issues ...")
    if len(SAMPLES_WITH_ISSUES_DF) == 0:
        logger.info("PASSED ALL TEST, CONGRATULATIONS! 100% consistency!")
    else:
        logger.info("Oh well, more work to do with {0:.1f}% success rate".format((NRANDOMSAMPLES2TEST-n_issues2solve)/n_records*100))
        logger.info("*** Problematic samples listed in {} ***".format(TEST_ROOT+'/plasmid_accessions_with_issues.txt'))       


    assert len(SAMPLES_WITH_ISSUES_DF) == 0, "Test failed with {0:.1f}% success rate".format((NRANDOMSAMPLES2TEST-n_issues2solve)/n_records*100)


test_mean_and_multireplicon_frame(mode="list of sample ids")
test_mean_and_multireplicon_frame(mode="random sampling")