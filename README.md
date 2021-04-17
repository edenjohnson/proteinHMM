# proteinHMM

## Introduction
This program uses BioPython perform a BLASTp search over the internet, clustal Omega to perform Multiple Sequence Alignment (MSA), and BioPython to construct a phylogenetic tree for analysis of the FMR1 gene.

The data sets generated from the BLASTp search are used to train protein HMMs that perform log-Viterbi score calculations when running the full pipeline from the following file: 
```sh
protein_HMM_script.py
```

NOTE: This file refers to a java .jar file that MUST be present in the same directory as the project files in order to run the program from the script file. Otherwise, to perform all of the steps besides protein HMM training, run the following files in the following order:

```sh
1. BLAST.py
2. curate_BLAST_results.py
3. MSA.py
4. adjust_fasta_files.py
5. phylogenetic_tree.py

```

### Prerequisites/Assumptions

Dependencies:
* Python 3+ 
* BioPython
* MatPlotLib
* Clusta Omega .exe
* Java

## Author

**Eden Johnson** 
