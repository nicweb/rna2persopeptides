import argparse
import subprocess
import time

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

def run_external_cmd(command):
	#bad error handling at the moment
	#infos only printed when process finishes
	print "Executing Command %s" %' '.join(command)
	t0=time.clock()
	try:
		subprocess.check_output(command, stderr=subprocess.STDOUT)
	except subprocess.CalledProcessError as e:
		print e.output
	t1=time.clock()
	print "Execution took %i seconds"%(t1-t0)

def createindex(args):
#create directory and create star index
	if not (args['starindex'] and args['annotation'] and args['fasta']):
		print "not enough arguments for creating"
		exit(1)
	print "creating directory"
	subprocess.check_output(["mkdir",args['starindex']])
	print "Creating Star Index"
	starcreateindex_cmd=["STAR", \
					"--runMode", "genomeGenerate", \
					"--runThreadN", "4", \
					"--genomeDir", args['starindex'], \
					"--genomeFastaFiles",args['fasta'], \
					"--sjdbGTFfile",args['annotation'], \
					"--sjdbGTFtagExonParentTranscript","Parent", \
					"-sjdbOverhang","100"]
	run_external_cmd(starcreateindex_cmd)


if (args['createindex'] != ''):
	print "calling creating index subroutine"
	createindex(args)


def star1stpass(args):
	star1stpass_cmd=["STAR", \
           "--genomeDir",args['starindex'], \
           "--readFilesIn",args['r1'],args['r2'], \
           "--runThreadN","4", \
           "--outFilterMultimapScoreRange","1", \
           "--outFilterMultimapNmax","20", \
           "--outFilterMismatchNmax","10", \
           "--alignIntronMax","500000", \
           "--alignMatesGapMax","1000000", \
           "--sjdbScore","2", \
           "--alignSJDBoverhangMin","1", \
           "--genomeLoad","NoSharedMemory", \
           "--readFilesCommand","cat", \
           "--outFilterMatchNminOverLread","0.33", \
           "--outFilterScoreMinOverLread","0.33", \
           "--sjdbOverhang","100", \
           "--outSAMstrandField","intronMotif", \
           "--outSAMtype","None", \
           "--outSAMmode","None", \
           "--outFileNamePrefix",args['outputdir']]
	print "calling 1.pass Star Mapping subroutine"
	run_external_cmd(star1stpass_cmd)
		   
def star2ndpass(args):
	star2ndpass_cmd=["STAR", \
             "--genomeDir","GENOME_TMP", \
             "--readFilesIn",args['r1'],args['r2'], \
             "--runThreadN","4", \
             "--outFilterMultimapScoreRange","1", \
             "--outFilterMultimapNmax","20", \
             "--outFilterMismatchNmax","10", \
             "--alignIntronMax","500000", \
             "--alignMatesGapMax","1000000", \
             "--sjdbScore","2", \
             "--alignSJDBoverhangMin","1", \
             "--genomeLoad","NoSharedMemory", \
             "--limitBAMsortRAM","70000000000", \
             "--readFilesCommand","cat", \
             "--outFilterMatchNminOverLread","0.33", \
             "--outFilterScoreMinOverLread","0.33", \
             "--sjdbOverhang","100", \
             "--outSAMstrandField","intronMotif", \
             "--outSAMattributes","NH HI NM MD AS XS", \
             "--outSAMunmapped","Within", \
             "--outSAMtype","BAM","SortedByCoordinate", \
             "--outSAMheaderHD","@HD VN:1.4", \
             "--outSAMattrRGline","ID:SM", \
             "--outFileNamePrefix",args['outputdir']]
	print "calling 2. pass star mapping subroutine"
	run_external_cmd(star2ndpass_cmd)

def starreindex(args):
	starreindex_cmd=["STAR", \
             "--runMode","genomeGenerate", \
             "--genomeDir","GENOME_TMP", \
             "--genomeFastaFiles",args['fasta'], \
             "--sjdbOverhang","100", \
             "--runThreadN","4", \
             "--sjdbFileChrStartEnd","SJ.out.tab", \
             "--outFileNamePrefix",args['outputdir']]
	print "calling star redindexing"
	run_external_cmd(starreindex_cmd)
	
def run_spadder():
	run_spladder_cmd=["python", \
              "~/tools/spladder/python/spladder.py", \
              "-b","Aligned.sortedByCoord.out.bam", \
              "-o","spladdrout", \
              "-a",args['annotation']]
	print "calling spladdr subroutine"
	run_external_cmd(run_spladder_cmd)
	
			  
run_star1stpass()
run_starreindex()
run_star2ndpass()
run_spladder()

