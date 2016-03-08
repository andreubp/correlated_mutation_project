from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline
import os


def exec_blast(infile, config_file, out_name):
	db, evalue = parse_blast_conf(config_file)
	fasta_string = SeqIO.read(infile, format="fasta")#
	result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string.seq)
	stdout, stderr = result_handle()
	name = out_name + ".xml"
	print (name)
	save_file = open(name, "w")
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

def parse_blast_conf(config_file):
	op_config = open(config_file, "r")

	for line in op_config:
		if line.startswith("clustalW"):
			line = line.split("\t")
			db = line[1]
			evalue = line[2]
	return (options)


def parse_seq_XML(blast_xml):
	hits = set()
	for record in NCBIXML.parse(blast_xml):
		for align in record.alignments:
			hit_id = align.hit_id
		blast_xml.close()
        
       
	for hit_id in hits:
		indexs = [i for i in xrange(len(hit_id)) if hit_id.find("|", n) == n]
		gi_id = hit_id[indexs[0]:indexs[1]]
		print(gi_id)


def clustalW(infil, config_file):
	clustalw2= r"/Applications/clustalw2"
	cline = ClustalwCommandline(clustalw2, infile=infil)
	print (cline)
	assert os.path.isfile(clustalw2), "Clustal W executable missing"
	stdout, stderr = cline()
