from Bio import SeqIO


def read_from(path: str, input_type: str):
    result = []
    with open(path) as input_handle:
        for record in SeqIO.parse(input_handle, input_type):
            print(record)
            result.append(record)
    return result


def convert_phy_fasta(phy_path, fasta_path):
    count = SeqIO.convert(phy_path, "phylip", fasta_path, "fasta")
    print("Converted %i records" % count)

# test
# read_from_phy("../data/Aligned sequences/dna.phy")
# convert_phy_fasta("../data/Aligned sequences/dna.phy", "./data/dna.fasta")
