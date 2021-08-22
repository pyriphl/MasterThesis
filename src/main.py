from cogent3 import available_distances

from fileIO import read_from, load_tree, write_fasta
from seq_operations import align_seq, calculate_distance_aligned_seq
from tree_operation import lookup_by_names

if __name__ == '__main__':
    records = read_from('data/Simphy/HGT/1/data_1.fasta', 'fasta')
    tree = load_tree('data/Simphy/HGT/1/g_trees1.trees', 'newick')

    taxa = tree.find_clades(terminal=True)
    tree_by_name = lookup_by_names(tree)

    seq1 = records[0].seq
    seq2 = records[1].seq
    aligns = align_seq(str(seq1), str(seq2))#TODO global MSA
    write_fasta(aligns,'data/Temp/alignments/test_1.fasta')
    # print(available_distances())
    distance = calculate_distance_aligned_seq('data/temp/alignments/test_1.fasta', 'gtr')

