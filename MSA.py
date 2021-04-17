# Import Clustal Omega wrapper
from Bio.Align.Applications import ClustalOmegaCommandline
import os

# open the basic and related data sets to create a combined.fa for MSA.
basic_file = open("basic_data_set.fa", 'r')
related_file = open("related_data_set.fa", 'r')
basic = basic_file.read()
related = related_file.read().split("\n")
remove_query = []

# remove the Homo sapiens FMR1 from the related file for the combined alignment file
for line in related:
    if "Homo " in line:
        break
    else:
        remove_query.append(line)
r = "\n".join(remove_query)

# Create the combined.fa file
with open("combined.fa", 'w') as f:
    f.write(basic)
    f.write(r)

basic_file.close()
related_file.close()


def multiple_sequence_alignment(in_file, out_file):
    """
       :brief Function performs a multiple sequence alignment on an input file and outputs the results to a fasta file.
       :param in_file: input file in the form of an FASTA file.
       :param out_file: output file in the form of an aligned FASTA file.
    """
    # Get the command for Clustal Omega
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, threads=8, verbose=True,
                                                 seqtype="protein")

    # adjust the executable command to run with the binary file in the directory
    command_list = str(clustalomega_cline).split()
    for n, i in enumerate(command_list):
        if i == "clustalo":
            command_list[n] = "./clustal-omega-1.2.3-macosx"

    # create a string of the command to input into the command line
    command = " ".join(command_list)

    # call the command
    os.system(command)


in_files = ["combined.fa", "basic_data_set.fa", "related_data_set.fa"]
out_files = ["aligned.fa", "basic_aligned.fa", "related_aligned.fa"]

# Perform MSA on the combined data sets, the basic data set, and the related data set.
for i in range(len(in_files)):
    multiple_sequence_alignment(in_files[i], out_files[i])