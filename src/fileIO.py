from Bio import SeqIO
from Bio import Phylo
from Bio import AlignIO


# type: fasta
# return: list of Biopython record
def read_from(path: str, input_type: str):
    result = []
    with open(path) as input_handle:
        for record in SeqIO.parse(input_handle, input_type):
            result.append(record)
    return result


# type:newic
# return: Biopython tree
def load_tree(path: str, type: str):
    g_tree = Phylo.read(path, type)
    return g_tree


def write_fasta(aln, path):
    AlignIO.write(aln, path, 'fasta')


def convert_phy_fasta(phy_path, fasta_path):
    count = SeqIO.convert(phy_path, "phylip", fasta_path, "fasta")
    print("Converted %i records" % count)

# test
# read_from_phy("../data/Aligned sequences/dna.phy")
# convert_phy_fasta("../data/Aligned sequences/dna.phy", "./data/dna.fasta")
