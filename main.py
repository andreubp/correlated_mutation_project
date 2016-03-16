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

parser.add_argument('-hm', '--heatmap',
			dest='heatmap',
			action='store',
			default="png",
			required=False,
			help='Heatmap program options: png or plotly')

args = parser.parse_args()

if args.infile2:
	root = parse_config(args.params, "root")
	prefix_output = root.rstrip() + args.outfile+"_1"
	prefix_output_2 = root.rstrip() +args.outfile+"_2"
	prefix = args.outfile+"_1"
	prefix_2 =args.outfile+"_2"
	sys.stderr.write("You have TWO protein sequence into the input and they are saved in all this program with these next prefixes:\n\t%s\n\t%s\n" %(prefix, prefix_2))

	sys.stderr.write("Executing Blast for the %s...\n" %(prefix))
	file1= exec_blast(args.infile1, args.params, prefix_output)
	sys.stderr.write("\tBlast finished correctly.\n")
	sys.stderr.write("Executing Blast for the %s...\n" %(prefix_2))
	file2= exec_blast(args.infile2, args.params, prefix_output_2)
	sys.stderr.write("\tBlast finished correctly.\n")

	multifasta1, multifasta2 = get_sequences(args.infile1, file1, args.outfile, args.params, blast_xml_2 = file2, input2 = args.infile2)
	sys.stderr.write("Running ClustalW for the %s...\n" %(prefix))
	align = clustalW(multifasta1, args.params, prefix_output+".aln")
	sys.stderr.write("\tClustalW finished correctly.\n")
	sys.stderr.write("Running ClustalW for the %s...\n" %(prefix_2))
	align2 = clustalW(multifasta2, args.params, prefix_output_2+".aln")
	sys.stderr.write("\tClustalW finished correctly.\n")
	transposed = read_clustalw(prefix_output+".aln")
	transposed_2 = read_clustalw(prefix_output_2+".aln")
	sys.stderr.write("Generating Mutual Information table... You could see it in a few seconds in the next file: %s\n" %(prefix_output+"_mi.tsv"))
	mi = mutual_information(transposed= transposed, transposed_2 = transposed_2)
	write_mi_output(mi, prefix_output + "_mi.tsv")
	sys.stderr.write("Plotting results...\n")
	plot_heatmap(mi, prefix_output+".png")

else:
	root = parse_config(args.params, "root")
	prefix_output = root.rstrip() + args.outfile
	prefix = args.outfile

	sys.stderr.write("You have only ONE protein sequence into the input with this next prefix:\n\t%s\n" %(prefix))

	sys.stderr.write("Executing Blast...\n")
	file1= exec_blast(args.infile1, args.params, prefix_output)
	sys.stderr.write("Blast finished correctly.\n")
	multifasta1 = get_sequences(args.infile1, file1, prefix_output,args.params)

	sys.stderr.write("Running ClustalW...\n")
	clustalW(prefix_output+".mfa", args.params, prefix_output+".aln")
	sys.stderr.write("ClustalW finished correctly.\n")
	module= read_clustalw(prefix_output+".aln")
	sys.stderr.write("Generating Mutual Information table. You could see it in a few seconds in the next file: %s\n" %(prefix_output+"_mi.tsv"))
	mi = mutual_information(module)
	write_mi_output(mi, prefix_output + "_mi.tsv")
	sys.stderr.write("Plotting results...\n")
	if args.heatmap == "plotly":
		plotly_heatmap(mi, prefix, args.params)
	elif args.heatmap == "png":
		plot_heatmap(mi,prefix_output+'.png')
	else:
		sys.stderr.write("This format output is not defined in this program. Try png or plotly")

	sys.stderr.write("All the files generated in this program are saved on:\n%s\n" %(prefix_output))

sys.stderr.write("The program is done. See you!\n")
