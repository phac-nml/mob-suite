from mob_suite.wrappers import detectCircularity


fasta_file = '/Users/jrobertson/Desktop/test.fasta'
out_dir = '/Users/jrobertson/Desktop/test_circularize'


c = detectCircularity()
c.run(fasta_file,out_dir)
