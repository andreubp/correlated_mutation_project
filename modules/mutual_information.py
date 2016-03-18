#!/usr/local/bin/python3

"""
This is the documentation for correlated_mutation_project. The functions below are
used in the execution of the different paths of main python. The main dependencies
of these module are the sys, os, math, numpy, potly and matplotlib packages but 
all of them should be present in a correct python instalation.
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
import plotly.plotly as py, plotly.graph_objs as go, plotly.tools as tls
import matplotlib.pyplot as plt
import numpy as np
from modules.parse_config import *



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
	if len(column_i) != len(column_j):
		raise IndexError("The two MSA should have the same number of related sequences (same species)")
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
	if entropy != 0.0:
		return -entropy
	else: 
		return entropy

def write_mi_output(mi, output, filtered):
	"""
	From the given list of list with the mutual information values writes a tab separated format output such that:

		position_i	position_j	mi_value

	When the input is a single protein both position columns will refer to the same protein. If two different proteins
	are given as input, the first position column will refer to positions in the first protein and the second column to
	positions in the second protein.
	If a cut off value is specified it will only retrieve pairs of positions with MI values higher than the cut off.
	"""

	op_outfile = open(output, "w")
	op_outfile.write("###correlated mutations output file\n###Contains mutual information values for each pair of sequences\n\n")
	
	if filtered == 0.0:
		op_outfile.write("#pos1\tpos2\tMI\n")
		for index1 in range(len(mi)):
			position = mi[index1]
			for index2 in range(len(position)):
				op_outfile.write("%d\t%d\t%.5f\n" %(index1 + 1, index2 + 1, position[index2]))
	else:
		op_outfile.write("###FILTERING CUT OFF: %.2f\n" %filtered)
		op_outfile.write("#pos1\tpos2\tMI\n")
		for index1 in range(len(mi)):
			position = mi[index1]
			for index2 in range(len(position)):
				if position[index2] >= filtered:
					op_outfile.write("%d\t%d\t%.5f\n" %(index1 + 1, index2 + 1, position[index2]))
	op_outfile.close()

def plot_heatmap(mi, output):
	"""
	Given a list with the MI scores for all possible pairs of residues in the protein(s) sequence(s)
	plots a heatmap using matplotlib with the MI scores for each pair and saves it in PDF format.
	The axis represent the positions in the sequence and the legend includes the color scale for MI values.
	"""
	fig = plt.figure()
	data = np.array(mi)
	fig, ax = plt.subplots()
	heatmap = ax.pcolor(data, cmap=plt.cm.jet)

	ax.invert_yaxis()
	ax.xaxis.tick_top()

	ax.set_xlabel('Seq 2')
	ax.set_ylabel('Seq 1')
	ax.xaxis.set_label_position('top')

	ax.set_xlim(0, len(mi[0]))
	ax.set_ylim(len(mi), 0)

	xmajor_ticks = np.arange(0, len(mi[1]), 20)
	xminor_ticks = np.arange(0, len(mi[1]), 1)

	ymajor_ticks = np.arange(0, len(mi), 20)
	yminor_ticks = np.arange(0, len(mi), 1)

	ax.tick_params(axis = 'both', which = 'major', labelsize = 5)
	ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)

	ax.set_xticks(xmajor_ticks)
	ax.set_xticks(xminor_ticks, minor = True)
	ax.set_yticks(ymajor_ticks)
	ax.set_yticks(yminor_ticks, minor = True)

	ax.tick_params(which = 'both', direction = 'out')

	plt.xticks(rotation=90)

	cb = plt.colorbar(heatmap)
	cb.set_label('MI value')

	fig.savefig(output, dpi = 700)

def plotly_heatmap(mi, output, config_file):
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
	plot_url = py.plot(data, filename = output)
