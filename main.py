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

parser = argparse.ArgumentParser(description="Correlated mutations")

parser.add_argument('-i1', '--input1',
			dest='infile1',
			action='store',
			default=None,
			required=True,
			help='Input file name')

parser.add_argument('-i2', '--input2',
			dest='infile2',
			action='store',
			default=None,
			required=False,
			help='Input file name')


parser.add_argument('-p', '--param',
			dest='params',
			action='store',
			default='parameters.config',
			required=False,
			help='Parameters configuration file')

parser.add_argument('-o', '--output',
			dest='outfile',
			action='store',
			default=None,
			required=True,
			help='Prefix for output files')

args = parser.parse_args()

if args.infile2:
	prefix_output = args.outfile+"1"
	prefix_output_2 = args.outfile+"2"

	#file1= exec_blast(args.infile1, args.params, prefix_output)
	#file2= exec_blast(args.infile2, args.params, prefix_output_2)
	#multifasta1, multifasta2= get_sequences(blast_xml = file1, output = args.outfile, blast_xml_2 = file2)
	#align = clustalW(multifasta1, args.params)
	#align2 = clustalW(multifasta2, args.params)
	transposed = read_clustalw("prova_doble_1.aln")
	transposed_2 = read_clustalw("prova_doble_2.aln")
	mi = mutual_information(transposed= transposed, transposed_2 = transposed_2)
	plot_heatmap(mi)

else:
	prefix_output = args.outfile
	print("Executing Blast...", file=sys.stderr)
	file1= exec_blast(args.infile1, args.params, prefix_output)
	print("Blast finished correctly.", file=sys.stderr)
	multifasta1 = get_sequences(args.infile1, file1, prefix_output,args.params)

	print("Runnning ClustalW...", file=sys.stderr)
	clustalW(prefix_output+".mfa", args.params)
	print("ClustalW finished correctly", file=sys.stderr)
	module= read_clustalw(prefix_output+".aln")
	print("Generating Mutual Information table", file=sys.stderr)
	mi = mutual_information(module)
	print("Ploting results...", file=sys.stderr)
	plot_heatmap(mi)
	#plotly_heatmap(mi)
	print("The program is done. See you!", file=sys.stderr)
