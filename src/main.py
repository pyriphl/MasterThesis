import matplotlib.pyplot as plt
import numpy
from cogent3.evolve import models
from Bio import Phylo

import tests
from src.distance_calculation import sliding_window_dist
from src.fileIO import delete_folder, load_tree, load_aln, write_data, load_distance, load_tags, create_dir
from src.path import ALN_PATH, SIMPHY_PATH, PICTURE_PATH, TXT_PATH
from src.show import plot_points_scatter, show_table, plot_surface, plot_histogram_3d, plot_histogram_2d, \
    plot_histogram_2d_onplanes, plot_histogram_2d_group, plot_histogram_2d_compare, plot_sliding_window
from src.sample_generation import sequence_generation_indelible, tree_generation_simphy
from src.seq_operations import get_by_name, calculate_distance_aligned_seq, dist_window_average, \
    WINDOW_SIZE, SLIDING_STEP
from src.tree_operation import lookup_by_names, pairwise_terminal_dist, pairwise_node_dist


def calculate_distances(tree, seqs):
    names, avg_model_distance, model_distances = dist_window_average(seqs, tree, 'JC69', WINDOW_SIZE)
    # input_seqs = prep_input_seq(seqs, tree, start, end)
    # names, model_distance = calculate_distance_aligned_seq(input_seqs, 'JC69')
    tree_distance = pairwise_terminal_dist(names, tree)
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
    # sequence_generation_indelible("data/Simphy/test/", "SimPhy_1.0.2/configuration_files/INDELible_complex.txt")

    s_tree = load_tree(SIMPHY_PATH + '01/s_tree.trees', 'newick')
    create_dir(PICTURE_PATH+'01/')
    Phylo.draw(s_tree,do_show=False)
    plt.savefig(PICTURE_PATH + '01/s_tree.png')
    # for i in range(10,11):
    #     num = f'{i:02d}'
    #     g_tree = load_tree(SIMPHY_PATH + '2/g_trees'+num+'v2.trees', 'newick')
    #     seqs = load_aln(SIMPHY_PATH + '2/data_'+num+'v2_TRUE.phy', 'phylip')
    #     names, model_distance, model_distance_list = sliding_window_dist(seqs,g_tree,WINDOW_SIZE)
    #     root_distance = pairwise_node_dist(names, g_tree)
    #     show_table(root_distance, names, 'Tree')
        # model_distance, tree_distance, names, model_distance_list = calculate_distances(g_tree, seqs)
        # show_table(model_distance, names, 'JC69')
        # show_table(tree_distance, names, 'Tree')
        # write_data(model_distance,tree_distance,names,TXT_PATH)
        # load_distance(TXT_PATH+'tree_dist.txt')
        # load_tags(TXT_PATH+'tags.txt')
        # Phylo.draw(g_tree)
        # plt.savefig(PICTURE_PATH +'2/g_tree'+num+'.png')
        # plot_histogram_2d_compare(model_distance, tree_distance, names, '2/data_'+num)
        # plot_sliding_window(model_distance_list,(len(seqs[0])-WINDOW_SIZE),SLIDING_STEP,names,PICTURE_PATH + '2/sliding_window/')


    # print(models.models)
