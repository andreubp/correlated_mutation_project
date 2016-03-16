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
import os
from Bio import  AlignIO
from Bio.Align.Applications import ClustalwCommandline
from modules.parse_config import *

def clustalW(infil, config_file, prefix):
	"""
	This method run ClustalW software and extract a multiple sequence alignment (MSA) from a multiple fasta file. We
	need to especify the path of the clustalW program in our computers in our configuration file.  The MSA is saved
	in a .aln file.
	"""
	#clustalw_path = parse_config(config_file, "clustalw")
	clustalw_path=r"/Applications/clustalw2"
	cline = ClustalwCommandline(clustalw_path, infile=infil, align="input", outfile=prefix, seqnos="ON", outorder="input", type="PROTEIN")
	assert os.path.isfile(clustalw_path), "Clustal W executable missing"
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
		if column[0] != "-":
			transposed_gap.append(column)
	"""	for el in column:
			if el == "-":
				gap += 1
		if gap <= (len(column)/2):
			transposed_gap.append(column)
	"""
	return transposed_gap
