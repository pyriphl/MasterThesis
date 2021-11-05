import math

import numpy

from src.seq_operations import get_by_name

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
        input_records = prepare_input_seqs(seqs, tree, 0, len(seqs[0])-1)
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
