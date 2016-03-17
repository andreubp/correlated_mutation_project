#!/usr/local/bin/python3

                    #####################################
                    #    correlated_mutation_project    #
                    #####################################

###############################################################################
#              version2: This is the last stable version of main.             #
# It includes a userfriendly paths to drive the function of the software.     #
# It includes the possibility of calculing the mutual information scores from #
# one or two sequence and extrapoled then the correlated mutation between the #
# different positions of the sequence/s. Finally it allows to plot a Heat-Map #
# about correlated mutation of each position on sequence.                     #
# This program accepts one or two inputs with FASTA format or raw sequence    #
# that allow to do the first step of the program, that's a Blast              #
###############################################################################

                        ########################
                        #        Authors:      #
                        #     Andreu Bofill    #
                        #    Marina Reixachs   #
                        ########################

########################
#        Modules       #
########################

import argparse
import sys
import os
from modules.blast import *
from modules.clustalw import *
from modules.mutual_information import *
from modules.parse_config import *

parser = argparse.ArgumentParser(description="Correlated mutations")

parser.add_argument('-i1', '--input1',
			dest='infile1',
			action='store',
			default=None,
			required=False,
			help='Input file name')

parser.add_argument('-i2', '--input2',
			dest='infile2',
			action='store',
			default=None,
			required=False,
			help='Input file name')

parser.add_argument('-mfa1', '--multifasta1',
			dest='multifasta_1',
			action='store',
			default=None,
			required=False,
			help='Multifasta Input file prepared to make a clustalW alignment')

parser.add_argument('-mfa2', '--multifasta2',
			dest='multifasta_2',
			action='store',
			default=None,
			required=False,
			help='Multifasta Input file prepared to make a clustalW alignment')

parser.add_argument('-c', '--config',
			dest='params',
			action='store',
			default='parameters.config',
			required=False,
			help='Parameters configuration file')

parser.add_argument('-o', '--output',
			dest='outfile',
			action='store',
			default=None,
			required=False,
			help='Prefix for output files')

parser.add_argument('-p', '--plot',
			dest='heatmap',
			action='store',
			default="png",
			required=False,
			help='Heatmap program options: png or plotly')

parser.add_argument('-f', '--filter',
			dest='filtered',
			action='store',
			default=0.0,
			required=False,
			type=float,
			help='If specified returns only positions with MI values higher than cutt-off value')
args = parser.parse_args()

def checkDependencies():
	"""
	Checks the import of the necessary python modules
	"""
	try:
		import Bio
		del Bio
	except ImportError:
		raise ImportError("ERROR: Unable to import Biopython")

	try:
		import numpy
		del numpy
	except ImportError:
		raise ImportError("ERROR: Unable to import Numpy")

	if args.heatmap == "plotly":
		try:
			import plotly
			del plotly
		except ImportError:
			raise ImportError("ERROR: Unable to import plotly")

	elif args.heatmap == "png":
		try:
			import matplotlib
			del matplotlib
		except ImportError:
			raise ImportError("ERROR: Unable to import Matplotlib")
	else:
		sys.stderr.write("ERROR: The format specified for the heatmap is no available. Try 'png' or 'plotly'\n")
		exit(1)

def clustalw_f(multifasta1,prefix_output):
	clustalW(multifasta1, args.params, prefix_output+".aln")
	module= read_clustalw(prefix_output+".aln")
	return (module)

def mi_f(prefix_output, prefix):
	root = parse_config(args.params, "root")
	module= read_clustalw(prefix_output+".aln")
	sys.stderr.write("Generating Mutual Information table...\n")
	mi = mutual_information(module)
	sys.stderr.write("You could see it in a few seconds in the next file:\t %s\n\n" %(prefix_output+"_mi.tsv"))
	write_mi_output(mi, prefix_output + "_mi.tsv", args.filtered)
	sys.stderr.write("Plotting results...\n")
	if args.heatmap == "plotly":
		plotly_heatmap(mi, prefix, args.params)
	elif args.heatmap == "png":
		plot_heatmap(mi,prefix_output+'.png')
	sys.stderr.write("\tThe plot has just been created. \n\n")
	sys.stderr.write("All the files generated in this program are saved on folder:\t%s\n" %(root))

def mi_two_f(module1, module2, prefix_output, prefix):
	root = parse_config(args.params, "root")
	sys.stderr.write("Generating Mutual Information table...\n")
	mi = mutual_information(transposed= module1, transposed_2 = module2)
	sys.stderr.write("You could see it in a few seconds in the next file: \t %s\n\n" %(prefix_output+"_mi.tsv"))
	write_mi_output(mi, prefix_output + "_mi.tsv", args.filtered)
	sys.stderr.write("Plotting results...\n")
	if args.heatmap == "plotly":
		plotly_heatmap(mi, prefix, args.params)
	elif args.heatmap == "png":
		plot_heatmap(mi,prefix_output+'.png')
	sys.stderr.write("\tThe plot has just been created. \n\n")
	sys.stderr.write("All the files generated in this program are saved on folder:\t%s\n" %(root))

def runCoevolution():
	root = parse_config(args.params, "root")
	if not args.infile2:
		if args.infile1:
			prefix_2 = None
			s = "only ONE protein sequence"
			if not args.outfile:
				outfile= input_name(args.infile1)
			else:
				outfile = args.outfile
		if args.multifasta_1:
			if not args.multifasta_2:
				prefix_2 = None
				s = "only ONE multifasta file"
				if not args.outfile:
					outfile= input_name(args.multifasta_1)
				else:
					outfile = args.outfile
			else:
				s = "TWO multifasta file"
				if not args.outfile:
					outfile= input_name(args.multifasta_1)+"__"+input_name(args.multifasta_2)
				else:
					outfile = args.outfile
				prefix_output_2 = root.rstrip() +outfile+"_2"
				prefix_2 = outfile+"_2"
		prefix_output = root.rstrip() + outfile + "_1"
		prefix = outfile+"_1"
		if prefix_2 != None:
			sys.stderr.write("You have %s into the input with these next prefixes:\t%s\t & \t %s\n\n" %(s, prefix, prefix_2))
		else:
			sys.stderr.write("You have %s into the input with this next prefix:\t%s\n\n" %(s, prefix))

		if args.infile1:
			if not args.multifasta_1:
				sys.stderr.write("Executing Blast...\n")
				file1= exec_blast(args.infile1, args.params, prefix_output)
				sys.stderr.write("\tBlast finished correctly.\n")
				multifasta1 = get_sequences(args.infile1, file1, prefix_output,args.params)
				sys.stderr.write("Running ClustalW...\n")
				clustalw_f(multifasta1, prefix_output)
				sys.stderr.write("\tClustalW finished correctly.\n\n")

				mi_f(prefix_output, prefix)
			else:
				print ("FUCK YOU!") #### RAISE ERROR BLBLABLA
		else:
			if args.multifasta_2:
				sys.stderr.write("Running ClustalW for the %s...\n" %(prefix))
				module1=clustalw_f(args.multifasta_1, prefix_output)
				sys.stderr.write("\tClustalW finished correctly.\n")
				sys.stderr.write("Running ClustalW for the %s...\n" %(prefix_2))
				module2=clustalw_f(args.multifasta_2, prefix_output_2)
				sys.stderr.write("\tClustalW finished correctly.\n\n")
				mi_two_f(module1, module2, prefix_output, prefix)
			else:
				multifasta1 = args.multifasta_1
				clustalw_f(multifasta1, prefix_output)
				mi_f(prefix_output, prefix)
	else:
		if not args.outfile:
			outfile = input_name(args.infile1) + "__" + input_name(args.infile2)
		else:
			outfile = args.outfile
		prefix_output_mfa = root.rstrip() + outfile
		prefix_output = root.rstrip() + outfile+"_1"
		prefix_output_2 = root.rstrip() +outfile+"_2"
		prefix = outfile+"_1"
		prefix_2 =outfile+"_2"
		sys.stderr.write("You have TWO protein sequence into the input and they are saved in all this program with these next prefixes:\t%s\t & \t%s\n\n" %(prefix, prefix_2))
		sys.stderr.write("Executing Blast for the %s...\n" %(prefix))
		file1= exec_blast(args.infile1, args.params, prefix_output)
		sys.stderr.write("\tBlast finished correctly.\n")
		sys.stderr.write("Executing Blast for the %s...\n" %(prefix_2))
		file2= exec_blast(args.infile2, args.params, prefix_output_2)
		sys.stderr.write("\tBlast finished correctly.\n")
		multifasta1, multifasta2 = get_sequences(args.infile1, file1, prefix_output_mfa, args.params, blast_xml_2 = file2, input2 = args.infile2)
		sys.stderr.write("Running ClustalW for the %s...\n" %(prefix))
		module1=clustalw_f(multifasta1, prefix_output)
		sys.stderr.write("\tClustalW finished correctly.\n\n")
		sys.stderr.write("Running ClustalW for the %s...\n" %(prefix_2))
		module2=clustalw_f(multifasta2, prefix_output_2)
		sys.stderr.write("\tClustalW finished correctly.\n\n")
		mi_two_f(module1, module2, prefix_output, prefix)

if __name__ == "__main__":
	checkDependencies()
	sys.stderr.write("\n\t\t CORRELATED MUTATIONS TOOL\n\n")
	sys.stderr.write("Dependencies OK.\n")
	runCoevolution()
	sys.stderr.write("\n\tThe program is done. See you!\n")
