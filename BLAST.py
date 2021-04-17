# Import modules
from Bio.Blast import NCBIWWW
from Bio import SeqIO

def perform_BLASTp_search(entrez_query, output_file, query_fasta):
    """
    :brief Function takes in an entrez query, an output file, and a query fasta file as input and performs
           a BLASTp search over the internet with the function parameters serving as query parameters.
    :param entrez_query: entrez style query to filter to include only certain organisms in the BLASTp search.
    :param output_file: XML file containing the BLASTp hit results
    :param query_fasta: query sequence in the form of a FASTA file.
    """

    print("\nBeginning query....")
    BLAST_results = NCBIWWW.qblast(program="blastp",
                                   database="nr",
                                   sequence=query_fasta.format("fasta"),
                                   hitlist_size=200,
                                   entrez_query=entrez_query,
                                   expect=0.05)

    with open(output_file, "w") as out_handle:
        out_handle.write(BLAST_results.read())

    print("\nWriting results to XML file...")
    BLAST_results.close()


entrez_queries = ["txid40674[ORGN]", "txid8782[ORGN]"] #mammals, birds
file_names = ["basic_data_set_results.xml", "related_data_set_results.xml"]
fmr1 = SeqIO.read("FMR1.fasta", format="fasta")

# Loop to create the separate basic and related data sets
for i in range(len(entrez_queries)):
    perform_BLASTp_search(entrez_queries[i], file_names[i], fmr1)