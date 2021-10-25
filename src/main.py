import matplotlib.pyplot as plt
import numpy
from cogent3.evolve import models
from Bio import Phylo

import tests
from src.fileIO import delete_folder, load_tree, load_aln, SIMPHY_PATH, ALN_PATH, PICTURE_PATH
from src.show import plot_points_scatter, show_table, plot_surface, plot_histogram_3d, plot_histogram_2d, \
    plot_histogram_2d_onplanes, plot_histogram_2d_group, plot_histogram_2d_compare, plot_sliding_window
from src.sample_generation import sequence_generation_indelible, tree_generation_simphy
from src.seq_operations import get_by_name, calculate_distance_aligned_seq, prep_input_seq, dist_window_average, \
    WINDOW_SIZE, SLIDING_STEP
from src.tree_operation import lookup_by_names, pairwaise_terminal_dist


def calculate_distances(tree,seqs):
    names, avg_model_distance, model_distances = dist_window_average(seqs, tree, 'JC69', WINDOW_SIZE)
    # input_seqs = prep_input_seq(seqs, tree, start, end)
    # names, model_distance = calculate_distance_aligned_seq(input_seqs, 'JC69')
    tree_distance = pairwaise_terminal_dist(names, tree)
    return avg_model_distance, tree_distance, names, model_distances


if __name__ == '__main__':
    # tests.test_load_tree()
    # tests.test_load_sequence()
    # tests.test_sample_generation()
    # tests.test_MSA_and_distance_models()
    # tests.test_tree_distance()

    # tree_generation_simphy(5, 1.5, 0.1, "data/Simphy/test/")
    # tree = load_tree('data/Simphy/test/01/' + 'g_trees10.trees', 'newick')
    # Phylo.draw(tree)
    # delete_folder("data/Simphy/test")
    # sequence_generation_indelible("data/Simphy/HGT/", "SimPhy_1.0.2/configuration_files/INDELible_simple.txt")

    # s_tree = load_tree(SIMPHY_PATH + '2/s_tree.trees', 'newick')
    # Phylo.draw(s_tree,do_show=False)
    # plt.savefig(PICTURE_PATH + '2/s_tree.png')
    for i in range(10,11):
        num = f'{i:02d}'
        g_tree = load_tree(SIMPHY_PATH + '2/g_trees'+num+'v2.trees', 'newick')
        seqs = load_aln(SIMPHY_PATH + '2/data_'+num+'v2_TRUE.phy', 'phylip')
        print(len(seqs[0]))
        model_distance, tree_distance, names, model_distance_list = calculate_distances(g_tree, seqs)
        show_table(model_distance, names, 'JC69')
        show_table(tree_distance, names, 'Tree')
        # Phylo.draw(g_tree,do_show=False)
        # plt.savefig(PICTURE_PATH +'2/g_tree'+num+'.png')
        # plot_histogram_2d_compare(model_distance, tree_distance, names, '2/data_'+num)
        plot_sliding_window(model_distance_list,(len(seqs[0])-WINDOW_SIZE),SLIDING_STEP,names,PICTURE_PATH + '2/sliding_window/')


    # print(models.models)
