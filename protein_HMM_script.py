# import the required modules
import os

# Call the different Python files for each step of the program
blast = "python BLAST.py"
curate_BLAST_Results = "python curate_BLAST_results.py"
msa = "python MSA.py"
adjust_fasta = "python adjust_fasta_files.py"
phylogenetic_tree = "python phylogenetic_tree.py"
build_HMM = "python build_HMM.py"

os.system(blast)
os.system(curate_BLAST_Results)
os.system(msa)
os.system(adjust_fasta)
os.system(phylogenetic_tree)
os.system(build_HMM)