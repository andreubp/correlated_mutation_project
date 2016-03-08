import os
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio import SeqIO, Entrez, AlignIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalwCommandline



def exec_blast(infile, config_file, out_name):
	db, evalue = parse_blast_conf(config_file)
	fasta_string = SeqIO.read(infile, format="fasta")#
	result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string.seq)
	name = out_name + ".xml"
	save_file = open(name, "w")
	save_file.write(result_handle.read())
	save_file.close()
	result_handle.close()
	return (name)

def parse_blast_conf(config_file):
	op_config = open(config_file, "r")

	for line in op_config:
		if line.startswith("blast"):
			line = line.split("\t")
			db = line[1]
			evalue = line[2]
	return(db, evalue)

def parse_seq_XML(blast_xml, output):
	print ("PARSE XML...")
	blast_xml_op = open (blast_xml, 'r')
	hits = set()
	outfile = output +".mfa"
	op_outfile = open(outfile, 'w')
	for record in NCBIXML.parse(blast_xml_op):
		print ("ROUND...")
		for align in record.alignments:
			hits.add(align.hit_id,)
#			for hsp in align.hsps:
#				print (hsp.expect)
#				print (hsp.query_end)
#				print (hsp.query_start)
			hit_id=align.hit_id.split("|")
#			if hsp.expect >= "0.05":
			efetch = Entrez.efetch(db="protein", id=hit_id[1], rettype="fasta")
			for line in efetch:
				id_info = line
				break
			sequence = ""
			for line in efetch:
				sequence += line.rstrip()
			organism = id_info[id_info.find("[") + 1:id_info.find("]")]
			organism = organism.split()
			if len(organism) != 1:
				species = str(organism[0] + "_" + organism[1])
			else:
				species = str(organism[0] + "_" + "sp.")
			op_outfile.write("> "+ species + "|"+ hit_id[1] + "| \n" + sequence + "\n")
	op_outfile.close()
	return (outfile)
	blast_xml_op.close()

def clustalW(infil):
	clustalw2= r"/Applications/clustalw2"
	cline = ClustalwCommandline(clustalw2, infile=infil, align="input", seqnos="ON", outorder="input", type="PROTEIN")
	print (cline)
	assert os.path.isfile(clustalw2), "Clustal W executable missing"
	stdout, stderr = cline()

def read_clustaw(clustalw_file):
	clustalw = open(clustalw_file, 'r')
	align = AlignIO.read(clustalw, "clustal")
	for record in align:
		print (record.seq, record.id)
