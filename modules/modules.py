import os
import math
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
	clustalw2= r"/usr/bin/clustalw"
	cline = ClustalwCommandline(clustalw2, infile=infil, align="input", seqnos="ON", outorder="input", type="PROTEIN")
	assert os.path.isfile(clustalw2), "Clustal W executable missing"
	stdout, stderr = cline()

def read_clustaw(clustalw_file):
	clustalw = open(clustalw_file, 'r')
	align = AlignIO.read(clustalw, "clustal")
	transposed = transpose_alignment(align)
	return transposed

def transpose_alignment(align):
	index = range(len(align[0]))
	transposed = list()
	for i in index:
		transposed.append(''.join([seq[i] for seq in align]))
	return transposed


def mutual_information(transposed):
	mi = dict()
	length = range(len(transposed))
	for i in length:
		#H(i) = -sum_x(P(x)ln2(P(x)))
		entropy_i = entropy(transposed[i])	

		for j in length:
			entropy_j = entropy(transposed[j])
			joint = joint_entropy(transposed[i], transposed[j])
			mi = entropy_i + entropy_j - joint
			print(i, transposed[i], j, transposed[j], entropy_i, entropy_j, joint, mi)
			


def entropy(column_string):
	frequencies = dict()
	total = len(column_string)
	entropy = 0
	for i in column_string:
		if i in frequencies:
			frequencies[i] +=1
		else:
			frequencies[i] = 1

	for key in frequencies:
		frequencies[key] /= total
		entropy += frequencies[key]*math.log(frequencies[key], 2)
	return -entropy


def joint_entropy(column_i, column_j):
	freq_ij = dict()
	total = len(column_i)
	entropy = 0
	for index in range(total):
		i = column_i[index]
		j = column_j[index]
		if i+j in freq_ij:
			freq_ij[i+j] +=1
		else:
			freq_ij[i+j] = 1

	for key in freq_ij:
		freq_ij[key] /= total
		entropy += freq_ij[key]*math.log(freq_ij[key], 2)
	return -entropy


