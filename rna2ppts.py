import argparse
import subprocess
import time

import rna2ppts

parser = argparse.ArgumentParser()
parser.add_argument("-r1", help="First RNA-FastQFile")
parser.add_argument("-r2", help="Second RNA-FastQFile")
parser.add_argument("-s", "--starindex", help="Path to starindex")
parser.add_argument("-f", "--fasta", help="Path to fasta reference file")
parser.add_argument("-o", "--outputdir",help="output dir")
parser.add_argument("-a", "--annotation",help="path to annotation")
parser.add_argument("-v", "--vcf",help="path to vcf")
parser.add_argument("--verbose", help="increase output verbosity",action="store_true")
parser.add_argument("--createindex", help="creates star index req. -s -f -a",action="store_true")
args=vars(parser.parse_args())
#print args

star_exec="STAR"
spladder_exec="/software/spladder/python/spladder.py"

	
t_start=time.time()			  
#run_star1stpass(args)
#run_starreindex(args)
run_star2ndpass(args)
run_spladder(args)
#do peptide stuff

t_end=time.time()
print "all done in %i seconds"%(t_end-t_start)
