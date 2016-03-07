from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Blast import NCBIXML

def exec_blast(infile, config_file):
	db, evalue = parse_blast_conf(config_file)
	blastp_cline = blastp(query=infile, db=db, evalue = evalue, outfmt=5, out="blast_out.xml")
	stdout, stderr = blastp_cline()
	
def parse_blast_conf(config_file):
	op_config = open(config_file, "r")
	
	for line in op_config:
		if line.startswith("blast"):
			line = line.split("\t")
			db = line[1]
			evalue = line[2]
	return(db, evalue)
