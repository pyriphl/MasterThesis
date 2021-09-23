import numpy
from cogent3.evolve import models
from Bio import Phylo

import tests
from src.fileIO import delete_folder, load_tree, load_aln, SIMPHY_PATH, ALN_PATH
from src.show import plot_points, show_table
from src.sample_generation import sequence_generation_indelible, tree_generation_simphy
from src.seq_operations import get_by_name, calculate_distance_aligned_seq, prep_input_seq
from src.tree_operation import lookup_by_names, pairwaise_terminal_dist


def calculate_distances():
    tree = load_tree(SIMPHY_PATH + 'g_trees1.trees', 'newick')
    seqs = load_aln(SIMPHY_PATH + 'data_1_TRUE.phy', 'phylip')
    input_seqs = prep_input_seq(seqs, tree)
    names, model_distance = calculate_distance_aligned_seq(input_seqs, 'JC69')
    tree_distance = pairwaise_terminal_dist(names, tree)
    return model_distance.to_array(), tree_distance, names

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

    model_distance, tree_distance, names = calculate_distances()
    show_table(model_distance, names, 'JC69')
    show_table(tree_distance, names, 'Tree')
    # plot_graph(tree_distance, 'cubic')
    # plot_points(gtr_distance, tree_distance, names)

    # print(models.models)