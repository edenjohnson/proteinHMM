# import the required modules
import os

# Call the different Python files for each step of the program
blast = "python3 BLAST.py"
curate_BLAST_Results = "python3 curate_BLAST_results.py"
msa = "python3 MSA.py"
adjust_fasta = "python3 adjust_fasta_files.py"
phylogenetic_tree = "python3 phylogenetic_tree.py"
build_HMM = "python3 build_HMM.py"

os.system(blast)
os.system(curate_BLAST_Results)
os.system(msa)
os.system(adjust_fasta)
os.system(phylogenetic_tree)
os.system(build_HMM)
