import math

import numpy

from src.seq_operations import get_by_name
from src.tree_operation import get_partitions


def sliding_window_dist(seqs, tree, window_size):
    dists = []
    names = []
    if window_size > len(seqs[0]):
        print('length out of boundary')
        return
    if window_size == 0:
        print('window size 0')
        return
    if window_size == len(seqs[0]):
        input_records = prepare_input_seqs(seqs, tree, 0, len(seqs[0]) - 1)
        names = [x[0] for x in input_records]
        input_seqs = [x[1].seq for x in input_records]
        dist = jukes_cantor(input_seqs)
        dists.append(dist)
    else:
        for index in range(0, len(seqs[0]) - window_size):
            start = index
            end = index + window_size
            input_records = prepare_input_seqs(seqs, tree, start, end)
            names = [x[0] for x in input_records]
            input_seqs = [x[1].seq for x in input_records]
            dist = jukes_cantor(input_seqs)
            dists.append(dist)
    return names, numpy.mean(dists, axis=0), dists


def jukes_cantor(seqs):
    # Jukes Cantor distance formula: (-3/4)ln[1-p*(4/3)]
    n = len(seqs)
    ps = percent_difference_of_nucleotides(seqs)
    results = numpy.zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            if i == j:
                d = 0
            else:
                p = ps[i][j]
                d = -0.75 * math.log(1 - (p * 4 / 3))
            results[i][j] = d
    return results


def prepare_input_seqs(seqs, tree, s_position, e_position):
    input_seqs = {}
    taxa = tree.get_terminals()
    for t in taxa:
        seq = get_by_name(seqs, t.name)
        input_seqs[t.name] = seq[s_position:e_position]

    return sorted(input_seqs.items())


def percent_difference_of_nucleotides(seqs, nucleobases=set('ACGT')):
    # percentage of nucleotide difference in two sequences
    n = len(seqs)
    results = numpy.zeros((n, n))
    for i in range(0, n):
        seq1 = seqs[i]
        for j in range(0, n):
            seq2 = seqs[j]
            diff_count = 0  # number of nucleotide differences
            valid_nucleotides_count = 0.0  # number of valid nucleotides (value is float for computing percentage)
            for a, b in zip(seq1, seq2):
                if a in nucleobases and b in nucleobases:
                    valid_nucleotides_count += 1
                    if a != b:
                        diff_count += 1
            results[i][j] = diff_count / valid_nucleotides_count if valid_nucleotides_count else 0
    return results


def get_all_partition(tree1, tree2):
    terminals1 = tree1.get_terminals()
    terminals2 = tree2.get_terminals()
    n1 = len(terminals1)
    n2 = len(terminals2)
    names1 = [n.name for n in terminals1]
    names2 = [n.name for n in terminals2]
    partition_all = []
    if n1 > n2:
        n = n1
        names = names1
    else:
        n = n2
        names = names2
    for i in range(0, n):
        partition_all.append((names[i], 'none'))
        partition_all.append(('none', names[i]))
        for j in range(i + 1, n):
            partition_all.append((names[j], names[i]))
            partition_all.append((names[i], names[j]))
    return partition_all


def Kuhner_Felsenstein_dist(tree1, tree2):
    partition_all = get_all_partition(tree1, tree2)
    partition1 = get_partitions(tree1, partition_all)
    partition2 = get_partitions(tree2, partition_all)
    branch_lengths1 = []
    branch_lengths2 = []
    for pair in partition_all:
        branch_lengths1.append(partition1[pair])
        branch_lengths2.append(partition2[pair])
    diff = numpy.subtract(branch_lengths2, branch_lengths1) / 2  # because every branch appeared twice
    return numpy.sum(diff ** 2)
