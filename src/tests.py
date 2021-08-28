from Bio import Phylo

from src.fileIO import read_from, SIMPHY_PATH, load_tree, ALN_PATH, load_aln
from src.sample_generation import tree_generation_simphy, sequence_generation_indelible
from src.seq_operations import align_mult_seq, calculate_distance_aligned_seq
from src.tree_operation import lookup_by_names


# test load newick tree from file
def test_load_tree():
    tree = load_tree(SIMPHY_PATH + 'g_trees1.trees', 'newick')
    taxa = tree.get_terminals()
    for taxon in taxa:
        print(taxon.name)
        print(taxon.branch_length)
    print(tree.distance(taxa[0], taxa[1]))
    Phylo.draw(tree)


# test load sequences from fasta file
def test_load_sequence():
    records = read_from(SIMPHY_PATH + 'data_1.fasta', 'fasta')
    print(records)


# test use MUSCLE and load MSA from file and calculate distance
def test_MSA_and_distance_models():
    align_mult_seq(SIMPHY_PATH + 'data_1.fasta', ALN_PATH + 'data_1.fasta')
    # delete_aln_file('data_1.fasta')
    seq = load_aln(ALN_PATH + 'data_1.fasta', 'fasta')
    pair_seq = {seq[0].name: seq[0].seq, seq[1].name: seq[1].seq}
    gtr_distance = calculate_distance_aligned_seq(pair_seq, 'GTR')
    print(gtr_distance)


# test tree and sequence generation
def test_sample_generation():
    # tree_generation_simphy(5, 10, 0.25, "data/Simphy/test/")
    sequence_generation_indelible("data/Simphy/test/", "SimPhy_1.0.2/configuration_files/INDELible_simple.txt")
