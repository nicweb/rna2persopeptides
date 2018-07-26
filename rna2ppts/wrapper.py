def run_external_cmd(command):
	#bad error handling at the moment
	#infos only printed when process finishes
	print "Executing Command %s" %' '.join(command)
	t0=time.time()
	try:
		subprocess.check_output(command, stderr=subprocess.STDOUT)
	except subprocess.CalledProcessError as e:
		print "ERROR executing Job"
		print e.output
		exit(1)
	t1=time.time()
	print "Execution took %i seconds"%(t1-t0)

def run_createindex(args):
#create directory and create star index
	if not (args['starindex'] and args['annotation'] and args['fasta']):
		print "not enough arguments for creating"
		exit(1)
	print "creating directory"
	subprocess.check_output(["mkdir",args['starindex']])
	print "Creating Star Index"
	starcreateindex_cmd=[star_exec, \
					"--runMode", "genomeGenerate", \
					"--runThreadN", "4", \
					"--genomeDir", args['starindex'], \
					"--genomeFastaFiles",args['fasta'], \
					"--sjdbGTFfile",args['annotation'], \
					"--sjdbGTFtagExonParentTranscript","Parent", \
					"-sjdbOverhang","100"]
	run_external_cmd(starcreateindex_cmd)


if (args['createindex']):
	print "calling creating index subroutine"
	run_createindex(args)


def run_star1stpass(args):
	if not (args['starindex'] and args['r1'] and args['r2'] and args['fasta'] and args['outputdir']):
		print "not enough arguments for 1. star run"
		exit(1)
	star1stpass_cmd=[star_exec, \
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
           "--outSAMmode","None"]#, \
       #    "--outFileNamePrefix",args['outputdir']]
	print "calling 1.pass Star Mapping subroutine"
	run_external_cmd(star1stpass_cmd)
		   
def run_star2ndpass(args):
	star2ndpass_cmd=[star_exec, \
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
             "--outSAMattributes","NH","HI","NM","MD","AS","XS", \
             "--outSAMunmapped","Within", \
             "--outSAMtype","BAM","SortedByCoordinate", \
             "--outSAMheaderHD","@HD VN:1.4", \
             "--outSAMattrRGline","ID:SM"]#, \
           #  "--outFileNamePrefix",args['outputdir']]
	print "calling 2. pass star mapping subroutine"
	run_external_cmd(star2ndpass_cmd)

def run_starreindex(args):
	#TODO: Add annotation???
	print "creating directory"
	subprocess.check_output(["mkdir","GENOME_TMP"])
	starreindex_cmd=[star_exec, \
             "--runMode","genomeGenerate", \
             "--genomeDir","GENOME_TMP", \
             "--genomeFastaFiles",args['fasta'], \
             "--sjdbOverhang","100", \
             "--runThreadN","4", \
             "--sjdbFileChrStartEnd","SJ.out.tab"]#, \
         #    "--outFileNamePrefix",args['outputdir']]
	print "calling star redindexing"
	run_external_cmd(starreindex_cmd)
	
def run_spladder(args):
	print "indexing file..."
	subprocess.check_output(["samtools","index","Aligned.sortedByCoord.out.bam"])
	run_spladder_cmd=["python", \
              spladder_exec, \
              "-b","Aligned.sortedByCoord.out.bam", \
              "-o","spladdrout", \
              "-a",args['annotation']]
	print "calling spladder subroutine"
	run_external_cmd(run_spladder_cmd)