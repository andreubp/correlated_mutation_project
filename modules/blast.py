from Bio.Blast.Applications import NcbiblastpCommandline as blastp

def exec_blast(infile, config_file):
	db, evalue = parse_blast_conf(config_file)
	blastp_cline = blastp(query=infile, db=db, evalue = evalue, outfmt=5, out="")
	
def parse_blast_conf(config_file):
	
	return(db, evalue)
