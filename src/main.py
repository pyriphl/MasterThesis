import matplotlib.pyplot as plt
import numpy
from cogent3.evolve import models
from Bio import Phylo

import tests
from src.data_processing import normalize_minmax
from src.distance_calculation import sliding_window_dist
from src.fileIO import delete_folder, load_tree, load_aln, write_data, load_distance, load_tags, create_dir, \
    load_partitions_as_aln
from src.path import ALN_PATH, SIMPHY_PATH, PICTURE_PATH, TXT_PATH
from src.show import plot_points_scatter, show_table, plot_surface, plot_histogram_3d, plot_histogram_2d, \
    plot_histogram_2d_onplanes, plot_histogram_2d_group, plot_histogram_2d_compare, plot_sliding_window
from src.sample_generation import sequence_generation_indelible, tree_generation_simphy
from src.seq_operations import get_by_name, calculate_distance_aligned_seq, dist_window_average, \
    WINDOW_SIZE, SLIDING_STEP
from src.tree_operation import lookup_by_names, pairwise_terminal_dist, pairwise_node_dist

if __name__ == '__main__':

    # tree_generation_simphy(5, 1.5, 0.1, "data/Simphy/test/")
    # tree = load_tree('data/Simphy/test/01/' + 'g_trees10.trees', 'newick')
    # Phylo.draw(tree)
    # delete_folder("data/Simphy/test")
    # sequence_generation_indelible("data/Simphy/test/", "SimPhy_1.0.2/configuration_files/INDELible_complex.txt")

    for i in range(7, 8):
        num = f'{i:02d}'
        s_tree = load_tree(SIMPHY_PATH + num + '/s_tree.trees', 'newick')
        # create_dir(PICTURE_PATH)
        # create_dir(PICTURE_PATH + num + '/')
        # Phylo.draw(s_tree, do_show=False)
        # plt.savefig(PICTURE_PATH + num + '/s_tree.png')
        seqs = load_partitions_as_aln(SIMPHY_PATH + num + '/', 2, 'phylip')
        # print(seqs)
        g_tree1 = load_tree(SIMPHY_PATH + num + '/g_trees1.trees', 'newick')
        g_tree2 = load_tree(SIMPHY_PATH + num + '/g_trees2.trees', 'newick')
        # Phylo.draw(g_tree1, do_show=False)
        # plt.savefig(PICTURE_PATH + num + '/g_tree1.png')
        # Phylo.draw(g_tree2, do_show=False)
        # plt.savefig(PICTURE_PATH + num + '/g_tree2.png')
        # names, model_distance, model_distance_list = sliding_window_dist(seqs, g_tree1, WINDOW_SIZE)
        names, model_distance, model_distance_list = dist_window_average(seqs, g_tree1, 'JC69', WINDOW_SIZE)
        s_tree_names = [n[0] for n in names]
        # tree_distance = pairwise_terminal_dist(s_tree_names, s_tree)
        show_table(model_distance, names, 'JC69')
        # create_dir(PICTURE_PATH + num + '/sliding_window_JC/')
        # plot_sliding_window(model_distance_list, len(seqs[0])-WINDOW_SIZE, 1, names, PICTURE_PATH + num + '/sliding_window_JC/')
        # show_table(tree_distance, s_tree_names, 'Tree')
        node_distance = pairwise_node_dist(s_tree_names, s_tree)
        show_table(node_distance, s_tree_names, 'Node')
        node_dist_norm = normalize_minmax(0, 3, node_distance)
        show_table(node_dist_norm, s_tree_names, 'Node normalized')
        # plot_histogram_2d_compare(model_distance, tree_distance, s_tree_names, 'data_07')

    # print(models.models)
