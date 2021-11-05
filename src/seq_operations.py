import subprocess

import numpy
from cogent3 import load_aligned_seqs, get_model, make_aligned_seqs
from cogent3.evolve import distance
from Bio import Align


WINDOW_SIZE = 300
SLIDING_STEP = 1
# Models:
# 'JC69', 'K80', 'F81', 'HKY85', 'TN93', 'GTR', 'ssGN', 'GN', 'BH', 'DT', 'CNFGTR', 'CNFHKY', 'MG94HKY', 'MG94GTR',
# 'GY94', 'Y98', 'H04G', 'H04GK', 'H04GGK', 'GNC', 'DSO78', 'AH96', 'AH96_mtmammals', 'JTT92', 'WG01'
def calculate_distance_aligned_seq(input_seqs, model: str):
    aln = make_aligned_seqs(input_seqs, moltype="dna")
    d = distance.EstimateDistances(aln, submodel=get_model(model))
    d.run(show_progress=False)
    dists = d.get_pairwise_distances()
    names = dists.names
    return names, dists


def align_pair_seq(seq1, seq2):
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    return alignments


# def align_mult_seq(in_path: str, out_path: str):
# '>/dev/null 2>&1' is used to supress execution log on the command line
# subprocess.check_output('muscle -in ' + in_path + ' -out ' + out_path + '>/dev/null 2>&1', shell=True)
# os.system('muscle -in ' + in_path + ' -out ' + out_path)


def get_by_name(self, name: str):
    for seq in self:
        if seq.name == name:
            return seq


def slice_seq(seq, start, end):
    return seq[start:end]


def prep_input_seq(seqs, tree, s_position, e_position):
    input_seqs = {}
    taxa = tree.get_terminals()
    for t in taxa:
        seq = get_by_name(seqs, t.name)
        input_seqs[t.name] = seq[s_position:e_position]
    return input_seqs


def dist_window_average(seqs, tree, model: str, window_size: int):
    dists = []
    names = []
    if window_size > len(seqs[0]):
        print('length out of boundary')
        return
    if window_size == 0:
        print('window size 0')
        return
    if window_size == len(seqs[0]):
        input_seqs = prep_input_seq(seqs, tree, 0, len(seqs[0]))
        names, dist = calculate_distance_aligned_seq(input_seqs, model)
        dists.append(dist.to_array())
    else:
        for index in range(0, len(seqs[0]) - window_size):
            start = index
            end = index + window_size
            input_seqs = prep_input_seq(seqs, tree, start, end)
            names, dist = calculate_distance_aligned_seq(input_seqs, model)
            dists.append(dist.to_array())
    return names, numpy.mean(dists, axis=0), dists
