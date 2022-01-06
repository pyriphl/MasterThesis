import os
import shutil

import numpy as np
from Bio import SeqIO
from Bio import Phylo
from Bio import AlignIO
from cogent3 import load_aligned_seqs

from src.path import ALN_PATH


# type: fasta
# return: list of Biopython record
def load_from(path: str, input_type: str):
    result = []
    with open(path) as input_handle:
        for record in SeqIO.parse(input_handle, input_type):
            result.append(record)
    return result


# type:newic
# return: Biopython tree
def load_tree(path: str, input_type: str):
    g_tree = Phylo.read(path, input_type)
    return g_tree


# input phylip
# return: Biopython alinment
def load_aln(path: str, input_type: str):
    aln = AlignIO.read(path, input_type)
    return aln


# input true aln phylip
# num must be bigger than 1
# return: Biopython alinment
def load_partitions_as_aln(path: str, input_type: str):
    s_tree = load_tree(path + 's_tree.trees', 'newick')
    s_taxa = s_tree.get_terminals()
    # result = AlignIO.read(path + 'data_1_TRUE.phy', input_type)
    result = {}
    alns = []
    # num_seqs = len(result)*(num-1) # this is the number of pair of sequence partition
    for n in range(1, 3):
        aln = load_aligned_seqs(path + 'data_' + str(n) + '_TRUE.phy', moltype="dna", format=input_type)
        alns.append(aln)
    for i in range(0, len(s_taxa)):
        # result[i*num+m].id = alns[m][i].id
        name = str(alns[0].names[i])
        seq = alns[0].seqs[i] + alns[1].seqs[i]
        result[name] = seq
        # result.names[i*num+m] = str(alns[m].names[i]) + name_code

    return result


def write_fasta(aln, path):
    AlignIO.write(aln, path, 'fasta')


def convert_phy_fasta(phy_path, fasta_path):
    count = SeqIO.convert(phy_path, "phylip", fasta_path, "fasta")
    print("Converted %i records" % count)


def delete_aln_file(file_name):
    os.remove(ALN_PATH + file_name)


def delete_folder(path: str):
    shutil.rmtree(path)


# TEST only
def write_data(model_distance, tree_distance, tags, out_path):
    np.savetxt(out_path + 'model_dist.txt', model_distance, fmt='%.6f')
    np.savetxt(out_path + 'tree_dist.txt', tree_distance, fmt='%.6f')
    with open(out_path + 'tags.txt', 'w') as f:
        for item in tags:
            f.write('%s\n' % item)


def write_distance(distance, type: str, path):
    np.savetxt(path + type + '_dist.txt', distance, fmt='%.6f')


def write_tags(tags, path):
    with open(path + 'tags.txt', 'w') as f:
        for item in tags:
            f.write('%s\n' % item)


def write_distances(distances, out_path):
    create_dir(out_path)
    np.savez(out_path + 'window_distances.npz', *distances)


def load_distances(in_path):
    dict_dists = np.load(in_path + 'window_distances.npz')
    data = [dict_dists[key] for key in dict_dists]
    return data


def load_distance(in_path):
    distance = np.loadtxt(in_path)
    return distance


def load_tags(in_path):
    with open(in_path, 'r') as f:
        names = f.readlines()
    results = []
    for n in names:
        results.append(n.replace('\n', ''))
    return results


def create_dir(path):
    try:
        os.mkdir(path)
    except OSError:
        print("Creation of the directory %s failed" % path)
    else:
        print("Successfully created the directory %s " % path)
# test
# read_from_phy("../data/AlignedSequences/dna.phy")
# convert_phy_fasta("../data/AlignedSequences/dna.phy", "./data/dna.fasta")
