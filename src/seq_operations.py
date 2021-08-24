import subprocess
from cogent3 import load_aligned_seqs, get_model
from cogent3.evolve import distance
from Bio import Align


def calculate_distance_aligned_seq(input_seqs: str, model: str):
    if not check_input(input_seqs):
        return -1
    aln = load_aligned_seqs(input_seqs, moltype="dna", format="fasta")
    d = distance.EstimateDistances(aln, submodel=get_model(model))
    d.run(show_progress=False)
    dists = d.get_pairwise_distances()
    return dists


def align_pair_seq(seq1, seq2):
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    return alignments


def align_mult_seq(in_path: str, out_path: str):
    subprocess.check_output('muscle -in ' + in_path + ' -out ' + out_path + '>/dev/null 2>&1', shell=True)
    # os.system('muscle -in ' + in_path + ' -out ' + out_path)
    # cline = MuscleCommandline(input=in_path, out=out_path)
    # print(cline)


def check_input(input_seq:str) -> bool:
    return True
