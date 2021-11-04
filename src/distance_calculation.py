import math

import numpy


def jukes_cantor(seqs):
    # Jukes Cantor distance formula: (-3/4)ln[1-p*(4/3)]
    n = len(seqs)
    ps = percent_difference_of_nucleotides(seqs)
    results = numpy.zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            p = ps[i][j]
            d = -0.75 * math.log(1 - (p * 4 / 3)) if p else 0
            results[i][j] = d
    return results


def percent_difference_of_nucleotides(seqs, nucleobases=set('ACGT')):
    # percentage of nucleotide difference in two sequences
    n = len(seqs)
    diff_count = 0  # number of nucleotide differences
    valid_nucleotides_count = 0.0  # number of valid nucleotides (value is float for computing percentage)
    results = numpy.zeros((len(seqs), len(seqs)))
    for i in range(0,n):
        seq1 = seqs[i]
        for j in range(0,n):
            seq2 = seqs[j]
            for a, b in zip(seq1, seq2):
                if a in nucleobases and b in nucleobases:
                    valid_nucleotides_count += 1
                    if a != b: diff_count += 1
            results[i][j] = diff_count / valid_nucleotides_count if valid_nucleotides_count else 0
    return results
