import sys, requests
from Bio import SeqIO
import pickle
import os.path

server = "http://rest.ensembl.org"
ext = "/lookup/id"
headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

ids_lookup = {}

tmp_mappingfile="/mnt/data_leon/id_mapping.pickle"

def load_mapping(mappingfile):
	if os.path.isfile(mappingfile):
		pickle_off = open(mappingfile,"rb")
		emp = pickle.load(pickle_off)
	else:
		print "no valid mapping file found"
		emp ={}
	print "Loaded %i mappings from file"%len(emp)	
	return emp

def save_mapping(mappingfile,mapping)
	pickling_on = open(mappingfile"wb")
	pickle.dump(mapping, pickling_on)
	pickling_on.close()
	print "Saved %i mappings to file"%len(mapping)

def request_ids(query):
	q={}
	q['ids'] = query
	r = requests.post(server+ext, headers=headers, data=q)
		if not r.ok:
			r.raise_for_status()
			#sys.exit()
 	decoded = r.json()
	reply={}
	for key in decoded.keys:
		reply[key]=decoded[key]['description']
	return reply

def get_gene_description(ids):
	mapping = load_mapping(tmp_mappingfile)
	queries =[]
	for id in ids:
		if id not in mapping and id not in queries:
			queries.append(id)
			if len(queries) = 100:
				#r = requests.post(server+ext, headers=headers, data='{ "ids" : ["ENSG00000157764", "ENSG00000248378" ] }')
				reply = request_ids(queries)
				queries = []
				print "Queried %i and got %i ids back"%(len(queries),len(reply))
				mapping.update(reply)
	if len(queries) > 0:
		reply=request_ids(queries)
		mapping.update(reply)
	return mapping
	save_mapping(tmp_mappingfile,mapping)
	


def sequence_cleaner(fasta_file, min_length=0, por_n=100):
    # Create our hash table to add the sequences
    sequences={}
	sequence_ids={}
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq).upper()
        # Check if the current sequence is according to the user parameters
        if (len(sequence) >= min_length and
            (float(sequence.count("N"))/float(len(sequence)))*100 <= por_n):
        # If the sequence passed in the test "is it clean?" and it isn't in the
        # hash table, the sequence and its id are going to be in the hash
            sequence_ids[sequence] = seq_record.id
			if sequence not in sequences:
                sequences[sequence] = seq_record.description
				
				
       # If it is already in the hash table, we're just gonna concatenate the ID
       # of the current sequence to another one that is already in the hash table
            else:
                sequences[sequence] += "|" + seq_record.description


    # Write the clean sequences
	ids=[]
	for k in sequence_ids.keys():
		if sequence_ids[k] not in ids:
			ids.append(sequence_ids[k])
	mapping=get_gene_descritpion(ids)
	
	for k in sequences.keys():
		sequences[k] += "|" + mapping[sequences_ids[k]]
	
	
	
    # Create a file in the same directory where you ran this script
    with open("clear_" + fasta_file, "w+") as output_file:
        # Just read the hash table and write on the file as a fasta format
        for sequence in sequences:
            output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")

    print("CLEAN!!!\nPlease check clear_" + fasta_file)


userParameters = sys.argv[1:]

try:
    if len(userParameters) == 1:
        sequence_cleaner(userParameters[0])
    elif len(userParameters) == 2:
        sequence_cleaner(userParameters[0], float(userParameters[1]))
    elif len(userParameters) == 3:
        sequence_cleaner(userParameters[0], float(userParameters[1]),
                         float(userParameters[2]))
    else:
        print("There is a problem!")
except:
    print("There is a problem!")
