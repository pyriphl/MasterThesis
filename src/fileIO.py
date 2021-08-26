import os
import shutil

from Bio import SeqIO
from Bio import Phylo
from Bio import AlignIO

ALN_PATH = 'data/AlignedSequences/'
SIMPHY_PATH = 'data/Simphy/HGT/1/'


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
def load_tree(path: str, input_type: str):
    g_tree = Phylo.read(path, input_type)
    return g_tree


# MUSCLE default output is fasta
# return: Biopython alinment
def load_aln(path: str, input_type: str):
    aln = AlignIO.read(path, input_type)
    return aln


def write_fasta(aln, path):
    AlignIO.write(aln, path, 'fasta')


def convert_phy_fasta(phy_path, fasta_path):
    count = SeqIO.convert(phy_path, "phylip", fasta_path, "fasta")
    print("Converted %i records" % count)


def delete_aln_file(file_name):
    os.remove(ALN_PATH + file_name)
# test
# read_from_phy("../data/AlignedSequences/dna.phy")
# convert_phy_fasta("../data/AlignedSequences/dna.phy", "./data/dna.fasta")
