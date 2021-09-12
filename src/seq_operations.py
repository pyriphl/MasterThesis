import subprocess
from cogent3 import load_aligned_seqs, get_model, make_aligned_seqs
from cogent3.evolve import distance
from Bio import Align


def calculate_distance_aligned_seq(input_seqs, model: str):
    aln = make_aligned_seqs(input_seqs, moltype="dna")
    d = distance.EstimateDistances(aln, submodel=get_model(model))
    d.run(show_progress=False)
    dists = d.get_pairwise_distances()
    return dists


def align_pair_seq(seq1, seq2):
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    return alignments


# def align_mult_seq(in_path: str, out_path: str):
    # '>/dev/null 2>&1' is used to supress execution log on the command line
    # subprocess.check_output('muscle -in ' + in_path + ' -out ' + out_path + '>/dev/null 2>&1', shell=True)
    # os.system('muscle -in ' + in_path + ' -out ' + out_path)


def check_input(input_seqs):
    if len(input_seqs) != 2:
        return False
    return True


def get_by_name(self, name: str):
    for seq in self:
        if seq.name == name:
            return seq
