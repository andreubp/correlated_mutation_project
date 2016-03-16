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
	db, evalue = parse_config(config_file, "blast")
	fasta_string = SeqIO.read(infile, format="fasta")
	result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string.seq)
	output= out_name + ".xml"
	save_file = open(output, "w")
	save_file.write(result_handle.read())
	save_file.close()
	result_handle.close()
	return (output)

def parse_config (config_file, option):
	"""
	This method is used to read the configuration file to extract the necesary information for each part of the program
	"""
	op_config = open(config_file, "r")
	if option == "blast":
		for line in op_config:
			if line.startswith("blast"):
				line = line.split("\t")
				db = line[1]
				evalue = line[2]
				return(db, evalue)
	elif option == "clustalw":
		for line in op_config:
			if line.startswith ("clustalw"):
				line = line.split("\t")
				clustal_path = line[1]
				return (clustal_path)

	elif option == "plotly":
		for line in op_config:
			if line.startswith("plotly"):
				line = line.split("\t")
				username = line[1]
				api_key = line[2]
				return (username, api_key)

class BlastResult(object):
	def __init__(self, hit_id, organism, sequence, evalue, coverage):
		self.hit = hit_id
		self.evalue = evalue
		self.coverage = coverage
		self.species = organism
		self.sequence = sequence

	def trim_coverage(self, cut_off):
		if self.coverage >= cut_off:
			return True
		else:
			return False


def parse_blast_XML(blast_xml):
	"""
	Read the blast_xml file generated before and extract the sequence and the id of each sequence in Blast and save them to
	multiple fasta file. It will allow ClustalW to generate a Multiple Sequence Alignment from all these sequence extracted.
	"""
	blast_xml_op = open (blast_xml, 'r')
	for record in NCBIXML.parse(blast_xml_op):
		for align in record.alignments:
			hit_id = align.hit_id.split("|")
			prev_eval = 1
			coverage = align.length / 390 ######arreglar per posar longitud sequencia
			for hsp in align.hsps:
				if hsp.expect < prev_eval:
					prev_eval = hsp.expect
			efetch = Entrez.efetch(db="protein", id=hit_id, rettype="fasta")
			for line in efetch:
				line = line.rstrip()
				if line.startswith(">"):
					id_info = line
					sequence = ""
				else:
					sequence += line
			sequence += line

			organism = id_info[id_info.find("[") + 1:id_info.find("]")]
			organism = organism.split()
			if len(organism) != 1:
				species = str(organism[0] + "_" + organism[1])

			yield BlastResult(hit_id[1], species, sequence, prev_eval, coverage)


def get_sequences(blast_xml, output, blast_xml_2 = False):
	species = set()
	final_results = []
	for result in parse_blast_XML(blast_xml):
		if result.species not in species:
			final_results.append(result)
			species.add(result.species)
		else:
			[ result for element in final_results if element.species == result.species and result.evalue < element.evalue]

	if blast_xml_2 == False:
		outfile = output +".mfa"
		op_outfile = open(outfile, 'w')
		for element in final_results:
			sentence = "> "+ element.species + "|"+ element.hit + "| \n" + element.sequence + "\n"
			op_outfile.write(sentence)
		op_outfile.close()
		return (outfile)

	else:
		species_2 = set()
		final_results_2 = []
		for result in parse_blast_XML(blast_xml_2):
			if result.species not in species_2:
				final_results_2.append(result)
				species_2.add(result.species)
			else:
				[ result for element in final_results_2 if element.species == result.species and result.evalue < element.evalue]

		final_species = species.intersection(species_2)

		filtered_results = [element for element in final_results if element.species in final_species]
		filtered_results_2 = [element for element in final_results_2 if element.species in final_species]

		outfile1 = output +"_1.mfa"
		op_outfile1 = open(outfile1, 'w')
		for element in filtered_results:
			sentence = "> "+ element.species + "|"+ element.hit + "| \n" + element.sequence + "\n"
			op_outfile1.write(sentence)
		op_outfile1.close()

		outfile2 = output +"_2.mfa"
		op_outfile2 = open(outfile2, 'w')
		for element in filtered_results_2:
			sentence = "> "+ element.species + "|"+ element.hit + "| \n" + element.sequence + "\n"
			op_outfile2.write(sentence)
		op_outfile2.close()

		return (outfile1, outfile2)


def clustalW(infil, config_file):
	"""
	This method run ClustalW software and extract a multiple sequence alignment (MSA) from a multiple fasta file. We
	need to especify the path of the clustalW program in our computers in our configuration file.  The MSA is saved
	in a .aln file.
	"""
	clustalw_path= parse_config(config_file, "clustalw")
#	clustalw2= r clustalw_path
	clustalw2= r'/usr/bin/clustalw2'
	cline = ClustalwCommandline(clustalw2, infile=infil, align="input", seqnos="ON", outorder="input", type="PROTEIN")
	assert os.path.isfile(clustalw2), "Clustal W executable missing"
	stdout, stderr = cline()

def read_clustalw(clustalw_file):
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
	transposed_gap = list()
	transposed = list()
	for i in index:
		transposed.append(''.join([seq[i] for seq in align]))
	for column in transposed:
		gap = 0
		for el in column:
			if el == "-":
				gap += 1
		if gap <= (len(column)/2):
			transposed_gap.append(column)
	return transposed_gap

def mutual_information(transposed, transposed_2 = False):
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
		if transposed_2 == False:
			for j in length:
				entropy_j = entropy(transposed[j])
				joint = joint_entropy(transposed[i], transposed[j])
				mi_calc = entropy_i + entropy_j - joint
				mi_list.append(mi_calc)
			mi.append(mi_list)

		else:
			length_2 = range(len(transposed_2))
			for j in length_2:
				entropy_j = entropy(transposed_2[j])
				joint = joint_entropy(transposed[i], transposed_2[j])
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

	###check which seq belongs to each axe
	ax.set_xlabel('Seq 2')
	ax.set_ylabel('Seq 1')

	ax.set_xlim(0, len(mi[1]))
	ax.set_ylim(len(mi), 0)

	plt.xticks(rotation=90)

	cb = plt.colorbar(heatmap)
	cb.set_label('MI value')

	#pdf = PdfPages('heatmap.pdf')
	#pdf.savefig(fig)
	fig.savefig('heatmap.png')
	#pdf.close()

def plotly_heatmap(mi):
	"""
	Given a list with the MI scores for all possible pairs of residues in the protein(s) sequence(s)
	creates a plotly heatmap.
	"""
	(username, api_key)= parse_config(config_file, "plotly")
	tls.set_credentials_file(username=username, api_key=api_key)
	data = [ go.Heatmap(
			z=mi,
			x=[i for i in range(len(mi))],
			y= [i for i in range(len(mi))]) ]
	plot_url = py.plot(data, filename = 'mi_heatmap')
