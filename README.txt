                    #####################################
                    #    correlated_mutation_project    #
                    #####################################

                        ########################
                        #        Authors:      #
                        #     Andreu Bofill    #
                        #    Marina Reixachs   #
                        ########################


This python program studies correlated mutations in a single protein sequence or between two different protein sequences.

Given a single or two fasta files with the input sequences it runs Blast in order to find homologs to perform a MSA.
Homologs are filtered according to Blast evalue in order that in the final MSA only one sequence for each species will
be represented, the one with lower evalue. If two protein sequences are given the species are matched in order that in 
the two MSA contain the same species.
With the MSA information a score of mutual information is calculated as a measure of correlation between positions.
The scores for each position pair are provided in a tsv output file and also plotted as a heatmap that can be saved
as a png file or visualised online in plotly.

DEPENDENCIES:

The program runs with python3. Please check that it is installed:

	python3 --version

If it's not, you can install it the following way:

	sudo apt-get install python3

It also requires Biopython, Numpy and Matplotlib or Plotly depending on the output format of the input.



PARAMETERS CONFIGURATION:

A parameters configuration file is provided. The file has the following format:

#######afegir format i explicacio del config

It can be modified or you can provide another config file with the same format (see USAGE section).




USAGE:

main.py [-h] -i1 INFILE1 [-i2 INFILE2] [-p PARAMS] -o OUTFILE
               [-hm HEATMAP]

optional arguments:
  -h, --help            show this help message and exit
  -i1 INFILE1, --input1 INFILE1
                        Input file name
  -i2 INFILE2, --input2 INFILE2
                        Input file name
  -p PARAMS, --param PARAMS
                        Parameters configuration file
  -o OUTFILE, --output OUTFILE
                        Prefix for output files
  -hm HEATMAP, --heatmap HEATMAP
                        Heatmap program options: png or plotly



