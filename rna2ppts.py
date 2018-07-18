import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-r1", help="First RNA-FastQFile")
parser.add_argument("-r2", help="Second RNA-FastQFile")
parser.add_argument("-s", "--starindex", help="Second RNA-FastQFile")
parser.add_argument("-f", "--fasta", help="Second RNA-FastQFile")
parser.add_argument("-o", "--outputdir",help="Second RNA-FastQFile")
parser.add_argument("-a", "--annotation",help="Second RNA-FastQFile")
parser.add_argument("-v", "--vcf",help="Second RNA-FastQFile")
parser.add_argument("--verbose", help="increase output verbosity",
                            action="store_true")

parser.parse_args()

