from cogent3 import available_distances

from fileIO import read_from, load_tree, delete_aln_file, SIMPHY_PATH, ALN_PATH, load_aln
from seq_operations import calculate_distance_aligned_seq, align_mult_seq
from tree_operation import lookup_by_names

if __name__ == '__main__':
    records = read_from(SIMPHY_PATH + 'data_1.fasta', 'fasta')
    tree = load_tree(SIMPHY_PATH + 'g_trees1.trees', 'newick')

    taxa = tree.find_clades(terminal=True)
    tree_by_name = lookup_by_names(tree)

    align_mult_seq(SIMPHY_PATH + 'data_1.fasta', ALN_PATH + 'data_1.fasta')
    # delete_aln_file('data_1.fasta')
    seq = load_aln(ALN_PATH + 'data_1.fasta', 'fasta')
    pair_seq = {seq[0].name: seq[0].seq, seq[1].name: seq[1].seq}
    gtr_distance = calculate_distance_aligned_seq(pair_seq, 'GTR')
    print(gtr_distance)

