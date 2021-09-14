import numpy
from Bio import Phylo

import tests
from src.fileIO import delete_folder, load_tree, load_aln, SIMPHY_PATH, ALN_PATH
from src.show import plot_points, plot_graph, show_table
from src.sample_generation import sequence_generation_indelible, tree_generation_simphy
from src.seq_operations import get_by_name, calculate_distance_aligned_seq
from src.tree_operation import lookup_by_names


def hgt_edit_tree_distance():
    tree = load_tree(SIMPHY_PATH + 'g_trees1.trees', 'newick')
    name_tree = lookup_by_names(tree)
    taxa = tree.get_terminals()
    seqs = load_aln(SIMPHY_PATH + 'data_1_TRUE.phy', 'phylip')
    input_seqs = {}
    for t in taxa:
        input_seqs[t.name] = get_by_name(seqs, t.name)
    gtr_distance = calculate_distance_aligned_seq(input_seqs, 'GTR')
    names = gtr_distance.names
    # print(gtr_distance)
    tree_distance = numpy.zeros((len(names), len(names)))
    for i in range (0, len(names)):
        for j in range (0, len(names)):
                tree_distance[i][j] = tree.distance(name_tree.get(names[i]), name_tree.get(names[j]))
    # print(tree_distance)
    # Phylo.draw(tree)
    return gtr_distance, tree_distance, names

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
    gtr_distance, tree_distance, names = hgt_edit_tree_distance()
    print(tree_distance)
    show_table(gtr_distance, tree_distance, names)
    # plot_graph(gtr_distance.to_array(), 'cubic')
    plot_points(gtr_distance.to_array(), tree_distance)