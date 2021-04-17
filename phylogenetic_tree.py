from Bio import Phylo
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import *
import matplotlib
# for font size of the phylogenetic tree


def construct_tree(alignment_file, tree_output_file):
    """
    :brief Function constructs a phylogenetic tree based upon an input alignment file. It will output the
           tree in .phy XML format with labeled branches and their corresponding length.
           Function also displays the tree in a plot and in text using the draw_ascii() function
           utilized from the BioPython Phylo module.
    :param alignment_file: MSA alignment file in the form of a FASTA file.
    :param tree_output_file: Phylogenetic tree output file in the form of a XML file.
    """

    aln = AlignIO.read(alignment_file, 'fasta')
    calculator = DistanceCalculator("blosum62") # Distance matrix
    dm = calculator.get_distance(aln)
    constructor = DistanceTreeConstructor()
    starting_tree = constructor.nj(dm)  # neighbor-joining
    phy = starting_tree.as_phyloxml()
    Phylo.write(phy, tree_output_file, 'phyloxml')

    tree = Phylo.read(tree_output_file, "phyloxml")
    matplotlib.rc('font', size=7) # change the font size to allow for full image display
    tree.ladderize()
    Phylo.draw(tree)
    Phylo.draw_ascii(tree)


in_files = ["aligned_fixed.fa", "basic_set_fixed.fa", "related_set_fixed.fa"]
tree_files = ['combined_tree.xml', 'basic_tree.xml', 'related_tree.xml']

for i in range(len(in_files)):
    construct_tree(in_files[i],tree_files[i])