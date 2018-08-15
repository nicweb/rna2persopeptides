import argparse
import subprocess
import time
import os

from rna2ppts import wrapper
from rna2ppts import converter 
from rna2ppts import cleaner

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
parser.add_argument("--runall", help="run complete",action="store_true")
parser.add_argument("--clean", help="clean only",action="store_true")
parser.add_argument("--convert", help="convert only",action="store_true")
args=vars(parser.parse_args())
#print args

	
t_start=time.time()	

if (args['createindex']):
	print "calling creating index subroutine"
	run_createindex(args)

def convert_path(path):
	#dir = os.path.dirname(__file__)
	dir = os.getcwd()
	if os.path.isabs(path):
		print "input path is absolute"
	else:
		path = os.path.join(dir, path)
		print "absolute path is %s" % path
	return path
	
for i in ['r1','r2','starindex','fasta','outputdir','annotation','vcf']:
			if args[i] != None:
				args[i]=convert_path(args[i])
				

if args['runall']:
	wrapper.run_star1stpass(args)
	wrapper.run_starreindex(args)
	wrapper.run_star2ndpass(args)
	wrapper.run_spladder(args)
#do peptide stuff
if args['runall'] or args['convert']:
	converter.calc_proteins('spladdrout',args['fasta'],args['vcf'])
if args['runall'] or args['clean']:
	cleaner.sequence_cleaner("spladdrout/predected_genes.fa")
t_end=time.time()
print "all done in %i seconds"%(t_end-t_start)
