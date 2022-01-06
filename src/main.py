import matplotlib.pyplot as plt
import numpy
from cogent3 import load_aligned_seqs
from cogent3.evolve import models
from Bio import Phylo

import data
from src.data_processing import correlation, reduce_to_list, logistic_regression
from src.distance_compare import compare_sw_results_correlation, compare_sw_results_linreg
from src.distance_calculation import sliding_window_dist
from src.fileIO import delete_folder, load_tree, load_aln, write_data, load_distance, load_tags, create_dir, \
    load_partitions_as_aln, write_distances, load_distances
from src.path import ALN_PATH, SIMPHY_PATH, PICTURE_PATH, FILE_PATH, WINDOW_SIZE
from src.show import plot_points_scatter, show_table, plot_surface, plot_histogram_3d, plot_histogram_2d, \
    plot_histogram_2d_onplanes, plot_histogram_2d_group, plot_histogram_2d_compare, plot_sliding_window, \
    plot_correlation_sw, show_compare_sw_results_linreg, plot_boxplot
from src.sample_generation import sequence_generation_indelible, tree_generation_simphy
from src.seq_operations import get_by_name, calculate_distance_aligned_seq, dist_window_average, \
    SLIDING_STEP
from src.tree_operation import lookup_by_names, pairwise_terminal_dist, pairwise_node_dist

if __name__ == '__main__':
    # tree_generation_simphy(5, 1.5, 0.1, SIMPHY_PATH)
    # tree_generation_simphy(5, 1.5, 0.1, "data/Simphy/test/")
    # tree_generation_simphy(5, 1.5, 0.1, 'data/Simphy/predefined_tree/')
    # tree1 = load_tree('data/Simphy/predefined_tree/1/' + 'g_trees1.trees', 'newick')
    # tree2 = load_tree('data/Simphy/predefined_tree/1/' + 'g_trees2.trees', 'newick')
    # stree = load_tree('data/Simphy/predefined_tree/1/' + 's_tree.trees', 'newick')
    # Phylo.draw(tree1)
    # Phylo.draw(tree2)
    # Phylo.draw(stree)
    # delete_folder("data/Simphy/test")
    # sequence_generation_indelible(SIMPHY_PATH, "SimPhy_1.0.2/configuration_files/INDELible_complex.txt")
    # sequence_generation_indelible("data/Simphy/temp/", "SimPhy_1.0.2/configuration_files/INDELible_complex.txt")

    Xs_list = []
    data_list = []
    for i in range(1, 4):
        model = 'JC69'
        num = f'{i:02d}'
        s_tree = load_tree(SIMPHY_PATH + num + '/s_tree.trees', 'newick')
        seqs = load_partitions_as_aln(SIMPHY_PATH + num + '/', 'phylip')
        g_tree1 = load_tree(SIMPHY_PATH + num + '/g_trees1.trees', 'newick')
        g_tree2 = load_tree(SIMPHY_PATH + num + '/g_trees2.trees', 'newick')
        g_trees = [g_tree1, g_tree2]
        seq_length = len(list(seqs.values())[0])
        # print(seqs.keys())
        # print(seqs.values())
        # print(seq_length)
        # d = data.dataset(num, seqs, g_trees, s_tree, model)
        # d.save_trees()
        # d.calculate_distances_names()
        # create_dir(FILE_PATH)
        # d.save_dataset(FILE_PATH + num + '/')
        # d.load_dataset(FILE_PATH + num + '/')
        # data_list.append(d)
    # data_frame = data.dataframe(data_list)
    # data_frame.save_trees()
    # data_frame.calculate_pd_dataframe()
    # data_frame.write_pd_dataframe()
    #     ####################################################################
    #     model_distance = load_distance(TXT_PATH + num + '/' + 'model_dist.txt')
    #     tree_distance = load_distance(TXT_PATH + num + '/' + 'tree_dist.txt')
    #     model_distance_list = load_distances(FILE_PATH + num + '/')
    #     names = load_tags(FILE_PATH + num + '/' + 'tags.txt')
        ####################################################################
        # create_dir(PICTURE_PATH)
        # create_dir(PICTURE_PATH + num + '/')
        # Phylo.draw(s_tree, do_show=False)
        # plt.savefig(PICTURE_PATH + num + '/s_tree.png')
        # Phylo.draw(g_tree1, do_show=False)
        # plt.savefig(PICTURE_PATH + num + '/g_tree1.png')
        # Phylo.draw(g_tree2, do_show=False)
        # plt.savefig(PICTURE_PATH + num + '/g_tree2.png')
        # Phylo.draw(g_tree3, do_show=False)
        # plt.savefig(PICTURE_PATH + num + '/g_tree3.png')
        # names, model_distance, model_distance_list = dist_window_average(seqs, model, WINDOW_SIZE)
        # s_tree_names = [n[0] for n in names]
        # tree_distance = pairwise_terminal_dist(s_tree_names, s_tree)
        # show_table(model_distance, names, model)
        # create_dir(PICTURE_PATH + num + '/sliding_window_'+model+'/')
        # plot_sliding_window(model_distance_list, seq_length - WINDOW_SIZE, 1, names, PICTURE_PATH + num + '/sliding_window_'+model+'/')
        # show_table(tree_distance, s_tree_names, 'Tree')
        #############################load data########################################################################
        # create_dir(TXT_PATH)
        # create_dir(TXT_PATH + num + '/')
        # write_distances(model_distance_list, TXT_PATH + num + '/')
        # write_data(model_distance, tree_distance, names, TXT_PATH + num + '/')
        #####################################################################################################
        # node_distance = pairwise_node_dist(s_tree_names, s_tree)
        # show_table(node_distance, s_tree_names, 'Node')
        # node_dist_norm = normalize_minmax(0, 4, node_distance)
        # show_table(node_dist_norm, s_tree_names, 'Node normalized')
        # plot_correlation_sw(node_dist_norm, model_distance_list, 0, 1.1)
        # plot_correlation_sw(tree_distance, model_distance_list, 0, 4)
        # compare_sw_results_correlation(model_distance_list, len(seqs[0])-WINDOW_SIZE, 1, names)
        # create_dir(PICTURE_PATH + num + '/linear_regression/')
        # plot_boxplot(model_distance_list, names)
        # y_preds, errors, slopes, name_pairs = compare_sw_results_linreg(model_distance_list, len(seqs[0])-WINDOW_SIZE, 1, names, PICTURE_PATH + num + '/linear_regression/')
        ######################prepare Xs###################
        # y_preds, errors, slopes, name_pairs = compare_sw_results_linreg(model_distance_list, len(seqs[0]) - WINDOW_SIZE,
        #                                                                 1, names,
        #                                                                 PICTURE_PATH + num + '/linear_regression/')
        # Xs_list.append(numpy.array((errors, slopes)))
        ###########################################
    #############################logistic regression#################################
    # Xs = Xs_list[0]
    # for j in range(1, len(Xs_list)):
    #     Xs = numpy.concatenate((Xs,Xs_list[j]), axis=1)
    # ys = [1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0]
    # ys = [1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1]
    # print(len(ys))
    # print(Xs.shape)
    # result = logistic_regression(Xs.T,ys, Xs.T)
    # print(result)
        # print(name_pairs)
        # y_test = [1,1,0,0,1,0,0,1,1,0]
        # result = logistic_regression(numpy.array(errors).reshape(-1,1),y_test,numpy.array(errors).reshape(-1,1))
        # print(result)
        # result = logistic_regression(numpy.array((errors,slopes)).T, y_test, numpy.array((errors,slopes)).T)
        # print(result)
    # print(models.models)
