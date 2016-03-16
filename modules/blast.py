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
import os, sys, math
from Bio import SeqIO, Entrez, AlignIO
from Bio.Blast import NCBIWWW, NCBIXML
from modules.parse_config import *


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

def exec_blast(infile, config_file, out_name):
	"""
	From a sequence input, and given a configuration file, execute Blast software from NCBI web server and print the result in a xml output file
	"""
	db, evalue = parse_config(config_file, "blast")
	fasta_string = SeqIO.read(infile, format="fasta")
	result_handle = NCBIWWW.qblast("blastp", db, fasta_string.seq)
	output= out_name + ".xml"
	save_file = open(output, "w")
	save_file.write(result_handle.read())
	save_file.close()
	result_handle.close()
	return (output)

def parse_blast_XML(blast_xml, config_file):
	"""
	Read the blast_xml file generated before and extract the sequence and the id of each sequence in Blast and save them to
	multiple fasta file. It will allow ClustalW to generate a Multiple Sequence Alignment from all these sequence extracted.
	"""
	blast_xml_op = open (blast_xml, 'r')
	Entrez.email = parse_config(config_file, "email")
	db, evalue = parse_config(config_file, "blast")

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

			if prev_eval <= float(evalue):
				yield BlastResult(hit_id[1], species, sequence, prev_eval, coverage)

def get_sequences(input1, blast_xml, output, config_file, blast_xml_2 = False, input2=False):
	species = set()
	final_results = []
	results_id = 0
	for result in parse_blast_XML(blast_xml, config_file):
		results_id += 1
		if result.species not in species:
			final_results.append(result)
			species.add(result.species)
		else:
			[ result for element in final_results if element.species == result.species and result.evalue < element.evalue]

	if blast_xml_2 == False:
		outfile = output +".mfa"
		op_outfile = open(outfile, 'w')
		infile = open(input1, 'r')
		for line in infile:
			op_outfile.write(line)
		op_outfile.write("\n")

		for element in final_results:
			sentence = "> "+ element.species + "|"+ element.hit + "| \n" + element.sequence + "\n"
			op_outfile.write(sentence)
		op_outfile.close()
		print("%s hits found. After filtering, we have %s hits." % (results_id,len(final_results)), file=sys.stderr)
		return (outfile)

	else:
		species_2 = set()
		final_results_2 = []
		results_id_2 = 0
		for result in parse_blast_XML(blast_xml_2, config_file):
			results_id_2 += 1
			if result.species not in species_2:
				final_results_2.append(result)
				species_2.add(result.species)
			else:
				[ result for element in final_results_2 if element.species == result.species and result.evalue < element.evalue]

		final_species = species.intersection(species_2)

		filtered_results = [element for element in final_results if element.species in final_species]
		filtered_results_2 = [element for element in final_results_2 if element.species in final_species]

		print("With the first protein, we found %s hits, and after filtering, %s." % (results_id,len(filtered_results)), file=sys.stderr)
		print("With the Second protein, we found %s hits, and after filtering, %s." % (results_id_2,len(filtered_results_2)), file=sys.stderr)

		outfile1 = output +"_1.mfa"
		op_outfile1 = open(outfile1, 'w')
		infile1 = open(input1, 'r')
		for line in infile1:
			op_outfile1.write(line)
		op_outfile1.write("\n")

		for element in filtered_results:
			sentence = "> "+ element.species + "|"+ element.hit + "| \n" + element.sequence + "\n"
			op_outfile1.write(sentence)
		op_outfile1.close()

		outfile2 = output +"_2.mfa"
		op_outfile2 = open(outfile2, 'w')
		infile2 = open(input2, 'r')
		for line in infile2:
			op_outfile2.write(line)
		op_outfile2.write("\n")

		for element in filtered_results_2:
			sentence = "> "+ element.species + "|"+ element.hit + "| \n" + element.sequence + "\n"
			op_outfile2.write(sentence)
		op_outfile2.close()

		return (outfile1, outfile2)
