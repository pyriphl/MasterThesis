from Bio import Phylo

from src.fileIO import read_from, SIMPHY_PATH, load_tree, ALN_PATH, load_aln
from src.sample_generation import tree_generation_simphy, sequence_generation_indelible
from src.seq_operations import calculate_distance_aligned_seq, get_by_name
from src.tree_operation import lookup_by_names


# test load newick tree from file
def test_load_tree():
    tree = load_tree(SIMPHY_PATH + 'g_trees1.trees', 'newick')
    taxa = tree.get_terminals()
    for taxon in taxa:
        print(taxon.name)
        print(taxon.branch_length)
    Phylo.draw(tree)


# calculate the tree distance and gtr distance for aligned pair of sequeces
def test_tree_distance():
    tree = load_tree(SIMPHY_PATH + 'g_trees1.trees', 'newick')
    taxa = tree.get_terminals()
    seqs = load_aln(ALN_PATH + 'data_1.fasta', 'fasta')
    input_seqs = {}
    for t in taxa:
        input_seqs[t.name] = get_by_name(seqs, t.name)
    gtr_distance = calculate_distance_aligned_seq(input_seqs, 'GTR')
    print("tree distance: " + str(tree.distance(taxa[0], taxa[1])))
    print("gtr distance: ")
    print(gtr_distance)
    # print(str(gtr_distance.array[0, 1]))


# test load sequences from fasta file
def test_load_sequence():
    records = read_from(SIMPHY_PATH + 'data_1.fasta', 'fasta')
    print(records)


# test use MUSCLE and load MSA from file and calculate distance
# def test_MSA_and_distance_models():
    # align_mult_seq(SIMPHY_PATH + 'data_1.fasta', ALN_PATH + 'data_1.fasta')
    # delete_aln_file('data_1.fasta')
    # seqs = load_aln(ALN_PATH + 'data_1.fasta', 'fasta')
    # pair_seq = {seqs[0].name: seqs[0].seq, seqs[1].name: seqs[1].seq}
    # gtr_distance = calculate_distance_aligned_seq(pair_seq, 'GTR')
    # print(gtr_distance)


# test tree and sequence generation
def test_sample_generation():
    # tree_generation_simphy(5, 10, 0.25, "data/Simphy/test/")
    sequence_generation_indelible("data/Simphy/HGT/", "SimPhy_1.0.2/configuration_files/INDELible_simple.txt")
