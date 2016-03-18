#!/usr/local/bin/python3

"""
This is the documentation for correlated_mutation_project. Both functions below are
used in the execution of the different paths of main python. The main dependency
of this module is only os packages.
This module contain documented functions used in all the other modules that this
program contain and also in the main script.
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

def parse_config (config_file, option):
	"""
	This method is used to read the configuration file to extract the necesary information for each part of the program
	"""
	op_config = open(config_file, "r")
	if option == "blast":
		for line in op_config:
			if line.startswith("blast"):
				line = line.split("\t")
				db = line[1].strip()
				evalue = line[2].strip()
				return(db, evalue)

	elif option == "clustalw":
		for line in op_config:
			if line.startswith ("clustalw"):
				line = line.split("\t")
				clustal_path = line[1].strip()
				return (clustal_path)

	elif option == "plotly":
		for line in op_config:
			if line.startswith("plotly"):
				line = line.split("\t")
				username = line[1].strip()
				api_key = line[2].strip()
				return (username, api_key)

	elif option == "email":
		for line in op_config:
			if line.startswith("Entrez_email"):
				line = line.split("\t")
				mail = line[1].strip()
				return (mail)
	elif option == "root":
		for line in op_config:
			if line.startswith("root"):
				line = line.split("\t")
				root = line[1].strip()
				return (root)

def input_name(inputname):
	"""
	This function is used to extract the name of the input file, to be used if an output name is not given
	"""
	path= (os.path.abspath(inputname))
	return (path.split("/")[-1].split(".")[0])
