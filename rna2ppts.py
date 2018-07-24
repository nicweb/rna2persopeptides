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
if (args['createindex'] != ''):
    if not (args['starindex'] and args['annotation'] and args['fasta']):
        print "not enough arguments for creating"
        exit(1)
    print "creating directory"
    subprocess.check_output(["mkdir",args['starindex']])
    print "Creating Star Index"
    starcreateindex=["STAR","--runMode","genomeGenerate","--runThreadN", "4","--genomeDir", args['starindex'],"--genomeFastaFiles",args['fasta'],"--sjdbGTFfile",args['annotation'],"--sjdbGTFtagExonParentTranscript","Parent","-sjdbOverhang","100"]
    try:
        subprocess.check_output(starcreateindex, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print e.output
