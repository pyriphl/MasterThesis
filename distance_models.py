import math
import random

import numpy
from Bio.Seq import Seq

N = 1000  # length of the seq
M = 2000  # random substitution freq 1/M
R = 10  # number of experiments


def test():
    result = [[]]
    for t in range(1, R):
        seq = gen_random_seq(N)
        seq1 = seq
        t += 1
        for i in range(1, M):
            n = random.randint(0, N - 1)
            seq1[n] = random.choice("CGTA")
            result[i][t] = seq - seq1
            mean = numpy.mean(result)

def JC_distance(p):  # p = observed proportion
    E = math.sqrt(p * (1 - p) / (N * ((1 - (4 / 3) * p) ** 2)))  # standard deviation
    d = -3 / 4 * math.log(1 - p * 4 / 3)  # distance
    print(E)
    print(d)

def gen_random_seq(length):
    letters = "CGTA"
    DNA = "".join(random.choices(letters, k=length))
    seq = Seq(DNA)
    return seq
JC_distance(90/948)