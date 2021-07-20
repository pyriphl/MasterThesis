from Bio import SeqIO


# with open("cor6_6.gb", "rU") as input_handle:
# with open("cor6_6.fasta", "w") as output_handle:
# sequences = SeqIO.parse(input_handle, "phylip")
# count = SeqIO.write(sequences, output_handle, "fasta")

def read_from_phy(path):
    with open(path) as input_handle:
        for record in SeqIO.parse(input_handle, "phylip"):
            print(record)


def convert_phy_fasta(phy_path, fasta_path):
    count = SeqIO.convert(phy_path, "phylip", fasta_path, "fasta")
    print("Converted %i records" % count)


read_from_phy("./data/dna.phy")
convert_phy_fasta("./data/dna.phy", "./data/dna.fasta")
