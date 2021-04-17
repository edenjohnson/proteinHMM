# Import modules
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def filter_BLAST_results(input_file, output_file):
    """
    :brief Functions parses through a XML BLASTp results file and selects 25 unique sequences and created a new data
           set fasta file that will serve as input for Multiple Sequence Alignment (MSA).
    :param input_file: input file in the form of an XML file.
    :param output_file: output file in the form of a FASTA file.
    """
    print("\nOpening the XML file to parse the BLAST results...\n")

    Entrez.email = "eden.johnson@sjsu.edu" # For Entrez query

    handle = open(input_file) # output in input file and get them in blast record format
    blast_records = NCBIXML.read(handle)
    data_set = []
    unique_species = []
    count = 0  # count the current iteration to ensure 25 sequences for each data set.

    for item in blast_records.alignments:
        start = item.title.find("[") + 1  # obtain the species name
        end = item.title.find("]")  # obtain the species name
        species = item.title[start:end]
        if species in unique_species:  # only collect sequences from unique species.
            pass
        else:
            if count == 25:  # only collect 25 unique species sequences
                break
            unique_species.append(species)
            accession = item.accession
            entrez = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
            fasta = entrez.read()
            header = fasta.split("\n")[0]  # get the
            sequence = Seq("".join(fasta.split("\n")[1:]))  #
            id_start = header.find(">") + 1  # Get the ID for the Storage Record
            id_end = header.find(" ")
            id = header[id_start:id_end]
            description = header[id_end + 1:]
            record = SeqRecord(seq=sequence, id=id, description=description)
            data_set.append(record)
            count += 1


    print("\nWriting data set to fasta file...")
    SeqIO.write(data_set, output_file, "fasta")
    handle.close()

    # Add the query sequence to the related results for MSA
    if "related" in input_file:
        query_file = open("FMR1.fasta", 'r')
        query = query_file.read()
        with open(output_file, 'a') as f:
            f.write(query)
        query_file.close()


input_files = ["basic_data_set_results.xml", "related_data_set_results.xml"]
output_files = ["basic_data_set.fa", "related_data_set.fa"]

# Curate the basic and related data sets
for i in range(len(input_files)):
    filter_BLAST_results(input_files[i], output_files[i])