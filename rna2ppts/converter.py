from BCBio import GFF
import operator
from Bio.SeqRecord import SeqRecord
import glob, os
from Bio import SeqIO
import time
import pysam

import gffutils
import pyfaidx
from pysam import VariantFile
from Bio.Seq import Seq

#sourcebase= '/mnt/test_data/hepavac34/rna_benign/spladdrout'

#reffile = '/mnt/test_data/refs/GRCh37.primary_assembly.genome.fa'
#vcffile= '/tmp/vcf.vcf.gz' #must be indexed!

fileending = 'gff3'



#db = gffutils.create_db('/tmp/ensg228794.gff', ':memory:')


outfile = 'predicted_exons.fa'
outfile2 = 'predicted_genes.fa'  



#Getting all gff-files in base directory
def get_gff_files(sourcebase):
    print "Checking for GFF Files in directory..."
    gfffiles = []
    os.chdir(sourcebase)
    for file in glob.glob("*.%s" %fileending):
        gfffiles.append(file)
    print "Found %i file(s) in %s ending with %s." %(len(gfffiles), sourcebase, fileending)
    return gfffiles

def get_gene_sequences(gene,db,vcfrecords,fasta):
	sequences = []
	geneseq = []
	num_mod=0
    #print gene
	for i in db.children(gene,featuretype='exon'):
		if vcfrecords != '':
			events=check_feature_for_vcfevent(i)
		else:
			events=[]
		if len(events) == 0:
			rec = i.sequence(fasta)#.translate(to_stop=True)
			sequences.append(SeqRecord(Seq(rec),i.id,i.attributes['Parent'][0],""))
			geneseq.append(rec)
		else:
			num_mod+=1
        #    sequences.append(Seq.translate(i.sequence(fasta),to_stop=True))
			rec = get_modified_sequence(i.sequence(fasta),events,i.start)
			sequences.append(SeqRecord(Seq(rec),i.id,i.attributes['Parent'][0],""))
			geneseq.append(rec)
            #sequences.append(get_modified_sequence(i.sequence(fasta),events,i.start))
            #return get_modified_sequence(i.sequence(fasta),events)
        #if i.id!='exon_7569':
        #    print "this is the record"
        #    print sequences[-1]
        
	return sequences,[SeqRecord(Seq(''.join(geneseq)),gene.id,"")], num_mod 

def check_feature_for_vcfevent(feature):
    vcfevents = []
    if feature.chrom not in vcfrecords.header.contigs:
        print "Chrom not found"
        return vcfevents
    for rec in vcfrecords.fetch(feature.chrom,feature.start,feature.stop):
        vcfevents.append(rec)
       #a=rec
    if len(vcfevents) == 0:
        return vcfevents
        print "no events in range detected (%s - %i on %i)"%(feature.chrom,feature.start,feature.stop)
    else:
        print "%i events in range detected (%s - %i on %i)"%(len(vcfevents),feature.chrom,feature.start,feature.stop)
        
    return vcfevents

def get_modified_sequence(in_seq,events,offset):
    print "Personalisizing sequence len:%i with %i events"%(len(in_seq),len(events))
    lengthchange=0
    for event in events:
        eventpos=event.start-offset
        print "Eventstart %i, offset: %i -> pos:%i"%(event.start,offset,eventpos)
        print "Char at  %i: %s (region +-3: %s)" %(eventpos,in_seq[eventpos],in_seq[eventpos-3:eventpos+3])
        alleles=event.alleles
        l = list(in_seq)
        
        if is_snp(alleles):
          #  print "Found SNP substituting %s with %s" %(l[eventpos],alleles[1])
            if l[eventpos]!= alleles[0]:
                print "warning SNP missmatch"
            l[eventpos]=alleles[1]
        if is_del(alleles):
            l[eventpos]=alleles[1] #TODO falsch!
            lengthchange-=len(alleles[1])-len(alleles[0])
            print "Found deletion: deleting %s from sequences"%alleles[0]
            print alleles
            print event
        if is_insert(alleles):
            l.insert(eventpos,alleles[1])
            lengthchange+=len(alleles[0])-len(alleles[1])
            print "found insert: inserting %s into sequence at pos %i"%alleles[1],eventpos
            print event
    new_seq=''.join(l)
    if (len(in_seq)-len(new_seq)) != lengthchange:
        print "Warning: length changed %i (should be %i)" %(len(in_seq)-len(new_seq),lengthchange)
    return new_seq 

def is_snp(e):
    if len(e[0]) == 1 and len(e[1]) == 1:
        return True
    return False

def is_insert(e):
    if len(e[0])<len(e[1]):
        return True
    return False

def is_del(e):
    if len(e[0])>len(e[1]):
        return True
    return False

def translate_records(records):
	new_records = []
	for rec in records:
		for i in range(0,3):
			new_records.append(SeqRecord(rec.seq[i:].translate(to_stop=True),"%s %s TL-Windows: %i Strand: +"%(rec.name,rec.id,i),rec.name,""))
			new_records.append(SeqRecord(rec.seq[i:].reverse_complement().translate(to_stop=True),"%s %s TL-Windows: %i Strand: - "%(rec.name,rec.id,i),rec.name,""))
	return new_records

def calc_proteins(sourcebase,reffile,vcffile):
	gfffiles = get_gff_files(sourcebase)
	if vcffile != None:
		print "VCFFile Provided"
		vcfrecords=VariantFile(vcffile)
	else:
		no_vcf=True
		vcfrecords=''
	print "DEBUG: reffile for pyfaidx is: %s"%reffile	
	fasta = pyfaidx.Fasta(reffile)
	
	
	for infile in gfffiles:
		total_genes=0
		total_exons=0
		total_genes_wrote=0
		total_exons_wrote=0
		total_mod=0
		db = gffutils.create_db(infile, ':memory:')
		for gene in db.features_of_type('gene'):
			total_genes+=1
			seqs,geneseq,num_mod=get_gene_sequences(gene,db,vcfrecords,fasta)
			total_exons+=len(seqs)
			total_mod+=num_mod
            #a=seqs
           # return seqs
            #return geneseq
			total_exons_wrote+=write_records(translate_records(seqs),'exons')
			total_genes_wrote+=write_records(translate_records(geneseq),'genes')
		print "processed (and wrote): %s with %i(%i) genes and %i(%i) exons (%i modified))" %(infile,total_genes,total_genes_wrote,total_exons,total_exons_wrote,total_mod)
        
    #ref_recs = load_reffile(reffile)
    
     #   records, records_until_stop = get_records(infile,ref_recs)
      #  write_records(records,records_until_stop,os.path.basename(infile))
        

def write_records(records,type):
    #translate records
    
    if type == 'exons':
        with open(outfile, "a") as out_handle:
            #print(records[1])
            return SeqIO.write(records, out_handle, "fasta")
    if type == 'genes':
        with open(outfile2, "a") as out_handle:
            #print(records[1])
            return SeqIO.write(records, out_handle, "fasta")