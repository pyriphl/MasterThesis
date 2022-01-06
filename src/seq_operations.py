from multiprocessing import Pool

import numpy
from cogent3 import get_model, make_aligned_seqs
from cogent3.evolve import distance
from Bio import Align

SLIDING_STEP = 1


# Models:
# 'JC69', 'K80', 'F81', 'HKY85', 'TN93', 'GTR', 'ssGN', 'GN', 'BH', 'DT', 'CNFGTR', 'CNFHKY', 'MG94HKY', 'MG94GTR',
# 'GY94', 'Y98', 'H04G', 'H04GK', 'H04GGK', 'GNC', 'DSO78', 'AH96', 'AH96_mtmammals', 'JTT92', 'WG01'
def calculate_distance_aligned_seq(input_seqs, model: str):
    aln = make_aligned_seqs(input_seqs, moltype="dna")
    d = distance.EstimateDistances(aln, submodel=get_model(model))
    d.run(show_progress=False)
    dists = d.get_pairwise_distances()
    # names = dists.names
    return dists


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


# def prep_input_seq(seqs, tree, s_position, e_position):
#     input_seqs = {}
#     taxa = tree.get_terminals()
#     for t in taxa:
#         seq = get_by_name(seqs, t.name)
#         input_seqs[t.name] = seq[s_position:e_position]
#     return input_seqs
def prep_input_seq(seqs, s_position, e_position):
    input_seqs = {}
    for name in seqs.keys():
        seq = seqs[name]
        input_seqs[name] = seq[s_position:e_position]
    return input_seqs

def dist_window_average(seqs, model: str, window_size: int):
    dists = []
    names = []
    seq_length = len(list(seqs.values())[0])
    if window_size > seq_length:
        print('length out of boundary')
        return
    if window_size == 0:
        print('window size 0')
        return
    if window_size == seq_length:
        input_seqs = prep_input_seq(seqs, 0, seq_length)
        names, dist = calculate_distance_aligned_seq(input_seqs, model)
        dists.append(dist.to_array())
    else:
        with Pool(processes=4) as pool:
            function_inputs = []
            for index in range(0, seq_length - window_size):
                start = index
                end = index + window_size
                input_seqs = prep_input_seq(seqs, start, end)
                function_inputs.append((input_seqs, model))
            dists = pool.starmap(calculate_distance_aligned_seq, function_inputs)
            array = [d.to_array() for d in dists]
            names = dists[0].names
            return names, numpy.mean(array, axis=0), array
            # for index in range(0, len(seqs[0]) - window_size):
            #     print(len(seqs[0].seq))
            #     print(len(seqs[0])-window_size)
            #     print(index)
            #     start = index
            #     end = index + window_size
            #     input_seqs = prep_input_seq(seqs, tree, start, end)
            #     names, dist = calculate_distance_aligned_seq(input_seqs, model)
            #     dists.append(dist.to_array())
    # return names, numpy.mean(dists, axis=0), dists
