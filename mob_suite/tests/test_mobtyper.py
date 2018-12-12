import sys, os,re
from mob_suite.mob_host_range import getHostRange, mob_host_range_main
import mob_suite.mob_typer


def test_getHostrange(capsys):
    getHostRange(["IncI1"],None,None,None,"exact")
    out, err = capsys.readouterr()
    assert "Enterobacterales" in out , "host range prediction was incorrect"

def test_mob_typer_host_range_multi_replicon():
   args=[
          "--infile", os.getcwd()+"/TestData/AB040415.fasta",
          "--outdir", "test_AB040415",
          "--host_range"
          ]
   sys.argv[1:] = args
   mob_suite.mob_typer.main()

   with open(file=os.getcwd()+"/test_AB040415/mobtyper_AB040415.fasta_report.txt") as fp_in:
       output = fp_in.readlines()

   assert any([len(re.findall("order\tEnterobacterales",out)) >= 1 for out in output]) == True, "Something went wrong with host range prediction";