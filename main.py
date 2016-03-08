import argparse
import sys
import os
import modules.modules as modules

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
			required=False,
			help='Prefix for output files')

args = parser.parse_args()

modules.exec_blast(args.infile1, args.params, "file1")

if args.infile2:
	modules.exec_blast(args.infile2, args.params, "file2")
#modules.clustalW(args.infile1, args.params)
