from cogent3 import available_distances

from fileIO import read_from, load_tree
from seq_operations import calculate_distance_aligned_seq, align_mult_seq
from tree_operation import lookup_by_names

if __name__ == '__main__':
    records = read_from('data/Simphy/HGT/1/data_1.fasta', 'fasta')
    tree = load_tree('data/Simphy/HGT/1/g_trees1.trees', 'newick')

    taxa = tree.find_clades(terminal=True)
    tree_by_name = lookup_by_names(tree)

    seqs = align_mult_seq('data/Simphy/HGT/1/data_1.fasta', 'data/Aligned sequences/data_1.fasta')
    # distance = calculate_distance_aligned_seq('data/Aligned sequences/data_1.fasta', 'GTR')

