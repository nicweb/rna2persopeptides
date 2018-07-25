import argparse
import subprocess

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
print args

def run_external_cmd(command):
    print "Executing Command %s" %''.join(command)
	try:
        subprocess.check_output(command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print e.output

def createindex(args):
#create directory and create star index
	if not (args['starindex'] and args['annotation'] and args['fasta']):
		print "not enough arguments for creating"
        exit(1)
    print "creating directory"
    subprocess.check_output(["mkdir",args['starindex']])
    print "Creating Star Index"
    starcreateindex=["STAR","--runMode","genomeGenerate","--runThreadN", "4","--genomeDir", args['starindex'],"--genomeFastaFiles",args['fasta'],"--sjdbGTFfile",args['annotation'],"--sjdbGTFtagExonParentTranscript","Parent","-sjdbOverhang","100"]
	run_external_cmd(starindex)


if (args['createindex'] != ''):
	createindex(args)


def star1stpass():
	star1pass=["STAR",
           "--genomeDir","$1",
           "--readFilesIn","$2", "$3"
           "--runThreadN","4",
           "--outFilterMultimapScoreRange","1",
           "--outFilterMultimapNmax","20",
           "--outFilterMismatchNmax","10",
           "--alignIntronMax","500000",
           "--alignMatesGapMax","1000000",
           "--sjdbScore","2",
           "--alignSJDBoverhangMin","1",
           "--genomeLoad","NoSharedMemory",
           "--readFilesCommand","cat",
           "--outFilterMatchNminOverLread","0.33",
           "--outFilterScoreMinOverLread","0.33",
           "--sjdbOverhang","100",
           "--outSAMstrandField","intronMotif",
           "--outSAMtype","None",
           "--outSAMmode","None",
           "--outFileNamePrefix","/path"]
		   
def star2ndpass():
	star2ndpass=["STAR",
             "--genomeDir","GENOME_TMP",
             "--readFilesIn","$1","$2",
             "--runThreadN","4",
             "--outFilterMultimapScoreRange","1",
             "--outFilterMultimapNmax","20",
             "--outFilterMismatchNmax","10",
             "--alignIntronMax","500000",
             "--alignMatesGapMax","1000000",
             "--sjdbScore","2",
             "--alignSJDBoverhangMin","1",
             "--genomeLoad","NoSharedMemory",
             "--limitBAMsortRAM","70000000000",
             "--readFilesCommand","cat",
             "--outFilterMatchNminOverLread","0.33",
             "--outFilterScoreMinOverLread","0.33",
             "--sjdbOverhang","100",
             "--outSAMstrandField","intronMotif",
             "--outSAMattributes","NH HI NM MD AS XS",
             "--outSAMunmapped","Within",
             "--outSAMtype","BAM","SortedByCoordinate",
             "--outSAMheaderHD","@HD VN:1.4",
             "--outSAMattrRGline","ID:SM",
             "--outFileNamePrefix","/path"]

def starreindex():
	starreindex=["STAR",
             "--runMode","genomeGenerate",
             "--genomeDir","GENOME_TMP",
             "--genomeFastaFiles","/mnt/test_data/refs/GRCh37.primary_assembly.genome.fa",
             "--sjdbOverhang","100",
             "--runThreadN","4",
             "--sjdbFileChrStartEnd","SJ.out.tab",
             "--outFileNamePrefix","/path"]

def run_spadder():
	run_spladder=["python",
              "~/tools/spladder/python/spladder.py",
              "-b","Aligned.sortedByCoord.out.bam",
              "-o","spladdrout",
              "-a",args['annotation']
	run_external_cmd(run_spladder)
	
			  
run_star1stpass()
run_starreindex()
run_star2ndpass()
run_spladder()

