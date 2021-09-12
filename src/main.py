import numpy
from Bio import Phylo

import tests
from src.fileIO import delete_folder, load_tree, load_aln, SIMPHY_PATH, ALN_PATH
from src.show import plot_points
from src.sample_generation import sequence_generation_indelible, tree_generation_simphy
from src.seq_operations import get_by_name, calculate_distance_aligned_seq


def hgt_edit_tree_distance():
    tree = load_tree(SIMPHY_PATH + 'g_trees1.trees', 'newick')
    taxa = tree.get_terminals()
    seqs = load_aln(SIMPHY_PATH + 'data_1_TRUE.phy', 'phylip')
    input_seqs = {}
    for t in taxa:
        input_seqs[t.name] = get_by_name(seqs, t.name)
    gtr_distance = calculate_distance_aligned_seq(input_seqs, 'GTR')
    # print(gtr_distance)
    tree_distance = numpy.zeros((len(taxa), len(taxa)))
    for i in range (0, len(taxa)):
        for j in range (0, len(taxa)):
            tree_distance[i][j] = tree.distance(taxa[i], taxa[j])
    # print(tree_distance)
    # Phylo.draw(tree)
    return gtr_distance.to_array(), tree_distance

if __name__ == '__main__':
    # tests.test_load_tree()
    # tests.test_load_sequence()
    # tests.test_sample_generation()
    # tests.test_MSA_and_distance_models()
    # tests.test_tree_distance()

    # align_mult_seq(SIMPHY_PATH + 'data_1.fasta', ALN_PATH + 'data_1.fasta')
    # tree_generation_simphy(5, 10, 0.25, "data/Simphy/test/")
    # delete_folder("data/Simphy/test")
    # sequence_generation_indelible("data/Simphy/HGT/", "SimPhy_1.0.2/configuration_files/INDELible_simple.txt")
    gtr_distance, tree_distance = hgt_edit_tree_distance()
    print(gtr_distance)
    print(tree_distance)
    plot_points(gtr_distance, tree_distance)