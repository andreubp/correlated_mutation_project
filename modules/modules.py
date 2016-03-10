#!/usr/local/bin/python3

"""
This is the documentation for correlated_mutation_project. The functions below are
used in the execution of the different paths of main python. The main dependencies
of these module are the sys, os, math and Bio packages but all of them should be
present in a correct python instalation.
This module only contain documented functions used in the main script.
The authors of this module are Andreu Bofill and Marina Reixachs.
"""

                        ########################
                        #        Authors:      #
                        #     Andreu Bofill    #
                        #    Marina Reixachs   #
                        ########################

########################
#        Modules       #
########################
import os, math, numpy as np
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio import SeqIO, Entrez, AlignIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalwCommandline
import plotly.plotly as py, plotly.graph_objs as go, plotly.tools as tls
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages

def exec_blast(infile, config_file, out_name):
	"""
	From a sequence input, and given a configuration file, execute Blast software from NCBI web server and print the result in a xml output file
	"""
	db, evalue = parse_blast_conf(config_file)
	fasta_string = SeqIO.read(infile, format="fasta")
	result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string.seq)
	output= out_name + ".xml"
	save_file = open(output, "w")
	save_file.write(result_handle.read())
	save_file.close()
	result_handle.close()
	return (output)

def parse_blast_conf(config_file):
	"""
	This method is used to read the configuration file to extract the necesary information for each part of the program
	"""
	op_config = open(config_file, "r")

	for line in op_config:
		if line.startswith("blast"):
			line = line.split("\t")
			db = line[1]
			evalue = line[2]
			return(db, evalue)

def parse_seq_XML(blast_xml, output):
	"""
	Read the blast_xml file generated before and extract the sequence and the id of each sequence in Blast and save them to
	multiple fasta file. It will allow ClustalW to generate a Multiple Sequence Alignment from all these sequence extracted.
	"""
	blast_xml_op = open (blast_xml, 'r')
	hits = set()
	outfile = output +".mfa"
	op_outfile = open(outfile, 'w')
	for record in NCBIXML.parse(blast_xml_op):
		for align in record.alignments:
			hits.add(align.hit_id,)
#			for hsp in align.hsps: #				print (hsp.expect) #				print (hsp.query_end) #				print (hsp.query_start)
			hit_id=align.hit_id.split("|")
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
			sentence = "> "+ species + "|"+ hit_id[1] + "| \n" + sequence + "\n"
			op_outfile.write(sentence)
	op_outfile.close()
	return (outfile)
	blast_xml_op.close()

def clustalW(infil):
	"""
	This method run ClustalW software and extract a multiple sequence alignment (MSA) from a multiple fasta file. We
	need to especify the path of the clustalW program in our computers in our configuration file.  The MSA is saved
	in a .aln file.
	"""
	clustalw2= r"/Applications/clustalw2"
	cline = ClustalwCommandline(clustalw2, infile=infil, align="input", seqnos="ON", outorder="input", type="PROTEIN")
	assert os.path.isfile(clustalw2), "Clustal W executable missing"
	stdout, stderr = cline()

def read_clustaw(clustalw_file):
	"""
	Read the MSA generated from ClustalW and save it to a new variable. This method calls another method, transpose_alignemnt
	mentioned below. It allows us to work with this alignemnt in a easier way.
	"""
	clustalw = open(clustalw_file, 'r')
	align = AlignIO.read(clustalw, "clustal")
	transposed = transpose_alignment(align)
	return transposed

def transpose_alignment(align):
	"""
	From a MSA, transpose all columns and rows, such that the columns in the alignment are saved as elements in a list. So
	finally we have a list of columns as strings.
	"""
	index = range(len(align[0]))
	transposed = list()
	for i in index:
		transposed.append(''.join([seq[i] for seq in align]))
	return transposed

def mutual_information(transposed):
	"""
	Calculates MI scores between all positions in a single protein or two different proteins
	from the transposed MSA columns and returns a list with the scores for all possible pair
	of positions.
	MI = sum_i(sum_j( H(i) + H(j) - H(i,j) ))
	"""
	mi = []
	length = range(len(transposed))
	for i in length:
		entropy_i = entropy(transposed[i])
		mi_list = []
		for j in length:
			entropy_j = entropy(transposed[j])
			joint = joint_entropy(transposed[i], transposed[j])
			mi_calc = entropy_i + entropy_j - joint
			mi_list.append(mi_calc)
		mi.append(mi_list)
	return mi

def entropy(column_string):
	"""
	Calculates the entropy of a single column(position) in the transposed MSA.
	H(i) = -sum_x(P(x)ln2(P(x)))
	"""
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
	"""
	Calculates the joint entropy for two columns(positions) in the transposed MSA.
	H(i,j) = -sum_x(sum_y( P(x,y)ln2(P(x,y)) ))
	"""
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

def plot_heatmap(mi):
	"""
	Given a list with the MI scores for all possible pairs of residues in the protein(s) sequence(s)
	plots a heatmap using matplotlib with the MI scores for each pair and saves it in PDF format.
	The axis represent the positions in the sequence and the legend includes the color scale for MI values.
	"""
	fig = plt.figure()
	data = np.array(mi)
	fig, ax = plt.subplots()
	heatmap = ax.pcolor(data, cmap=plt.cm.jet)

	ax.tick_params(direction='out')

	majorLocator   = MultipleLocator(20)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator   = MultipleLocator(1)

	ax.xaxis.set_major_locator(majorLocator)
	ax.xaxis.set_major_formatter(majorFormatter)
	ax.xaxis.set_minor_locator(minorLocator)

	ax.yaxis.set_major_locator(majorLocator)
	ax.yaxis.set_major_formatter(majorFormatter)
	ax.yaxis.set_minor_locator(minorLocator)

	ax.invert_yaxis()
	ax.xaxis.tick_top()

	ax.set_xlim(0, len(mi))
	ax.set_ylim(0, len(mi))

	plt.xticks(rotation=90)

	cb = plt.colorbar(heatmap)
	cb.set_label('MI value')

	pdf = PdfPages('heatmap.pdf')
	pdf.savefig(fig)
	pdf.close()




def plotly_heatmap(mi):
	"""
	Given a list with the MI scores for all possible pairs of residues in the protein(s) sequence(s)
	creates a plotly heatmap.
	"""
	tls.set_credentials_file(username="mars13", api_key="llj6ors56n")
	data = [ go.Heatmap(
			z=mi,
			x=[i for i in range(len(mi))],
			y= [i for i in range(len(mi))]) ]
	plot_url = py.plot(data, filename = 'mi_heatmap')
