import math
import random

import numpy
from Bio.Seq import Seq

# for test:
# calculation of genetic distance between two homologous sequences,
# that is the number of of substitutions that have accumulated between them since they diverged from their common ancestor.
def test():
    N = 1000  # length of the seq
    M = 2000  # random substitution freq 1/M
    R = 10  # number of experiments
    result = [0] * M
    mean_t = [0] * R
    for t in range(0, R - 1):
        seq = gen_random_seq(N)
        seq1 = seq
        print(seq)
        for i in range(0, M - 1):
            n = random.randint(0, N - 1)
            temp = list(str(seq1))
            temp[n] = random.choice("CGTA")
            old = str(seq)
            new = str(seq1)
            seq1 = Seq(''.join(temp))
            difference = sum(1 for a, b in zip(old, new) if a != b)
            result[i] = difference
            i += 1
        mean_t[t] = numpy.mean(result)
        t += 1
    mean = numpy.mean(result)
    print(JC_distance(mean / N, N))
    # print(JC_distance(90 / 948))


def JC_distance(p, N) -> (float, float):  # p = observed proportion
    if p >= 3 / 4:
        print("p should be smaller than 3/4")
        return -1, -1
    E = math.sqrt(p * (1 - p) / (N * ((1 - (4 / 3) * p) ** 2)))  # standard deviation
    d = -3 / 4 * math.log(1 - p * 4 / 3)  # distance
    return d, E


def gen_random_seq(length):
    letters = "CGTA"
    DNA = "".join(random.choices(letters, k=length))
    seq = Seq(DNA)
    return seq
