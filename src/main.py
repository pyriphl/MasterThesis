import numpy
from cogent3.evolve import models
from Bio import Phylo

import tests
from src.fileIO import delete_folder, load_tree, load_aln, SIMPHY_PATH, ALN_PATH
from src.show import plot_points_scatter, show_table, plot_surface, plot_histogram_3d, plot_histogram_2d, \
    plot_histogram_2d_onplanes, plot_histogram_2d_group, plot_histogram_2d_compare
from src.sample_generation import sequence_generation_indelible, tree_generation_simphy
from src.seq_operations import get_by_name, calculate_distance_aligned_seq, prep_input_seq, dist_window_average
from src.tree_operation import lookup_by_names, pairwaise_terminal_dist


def calculate_distances(tree,seqs):
    names, model_distance = dist_window_average(seqs, tree, 'JC69', len(seqs[0]))
    # input_seqs = prep_input_seq(seqs, tree, start, end)
    # names, model_distance = calculate_distance_aligned_seq(input_seqs, 'JC69')
    tree_distance = pairwaise_terminal_dist(names, tree)
    return model_distance, tree_distance, names


if __name__ == '__main__':
    # tests.test_load_tree()
    # tests.test_load_sequence()
    # tests.test_sample_generation()
    # tests.test_MSA_and_distance_models()
    # tests.test_tree_distance()

    tree_generation_simphy(5, 1.5, 0.1, "data/Simphy/test/")
    # tree = load_tree('data/Simphy/test/1/' + 'g_trees01.trees', 'newick')
    # Phylo.draw(tree)
    # delete_folder("data/Simphy/test")
    # sequence_generation_indelible("data/Simphy/HGT/", "SimPhy_1.0.2/configuration_files/INDELible_simple.txt")


    # s_tree = load_tree(SIMPHY_PATH + 's_tree.trees', 'newick')
    # g_tree = load_tree(SIMPHY_PATH + 'g_trees10.trees', 'newick')
    # seqs = load_aln(SIMPHY_PATH + 'data_10_TRUE.phy', 'phylip')
    # Phylo.draw(s_tree)
    # model_distance, tree_distance, names = calculate_distances(g_tree, seqs)
    # show_table(model_distance, names, 'JC69')
    # show_table(tree_distance, names, 'Tree')
    # Phylo.draw(g_tree)
    # plot_histogram_2d_compare(model_distance, tree_distance, names, 'data_01')

    # print(models.models)
