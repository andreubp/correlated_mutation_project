from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def exec_blast(infile, config_file):
	db, evalue = parse_blast_conf(config_file)
	fasta_string = SeqIO.read(infile, format="fasta")
	result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string.seq)
	save_file = open("blast_out.xml", "w")
	save_file.write(result_handle.read())
	save_file.close()
	result_handle.close()

def parse_blast_conf(config_file):
	op_config = open(config_file, "r")
	
	for line in op_config:
		if line.startswith("blast"):
			line = line.split("\t")
			db = line[1]
			evalue = line[2]
	return(db, evalue)
