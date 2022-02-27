from src import util
from skbio.alignment import StripedSmithWaterman
from skbio.alignment import local_pairwise_align_ssw, global_pairwise_align_nucleotide
from skbio.sequence import DNA
from src.exceptions import FrameshiftException, InternalStopCodonException

import pyopa

"""
Identifies correct gene usage for inputted sequences
"""


def new_align(seq: str, ref: dict):
    """
    Align sequences using a optimised striped Smith Waterman algorithm
    Args:
        seq: user inputted sequence to align against reference
        ref: dictionary of sequences

    Returns:
        Allele with the best alignment
    """

    data = {'gap_open': -20, 'gap_ext': -3, 'pam_distance': 150, }

    best = ''
    score = 0

    seq = pyopa.Sequence(seq)
    for k, v in ref.items():
        alignment = global_pairwise_align_nucleotide(DNA(seq.replace('N', '')), DNA(v))

    return best


def define_cdr3(seq: str, v: str, j: str):
    """
    Args:
        seq: user sequence for alignment
        v: variable cdr3 sequence
        j: joining segment cdr3 sequence

    Returns: cdr3 sequence of user sequence

    """

    alignment, score, start_end = local_pairwise_align_ssw(DNA(seq), DNA(v), gap_open_penalty=1)

    start = start_end[0][0]

    alignment, score, start_end = local_pairwise_align_ssw(DNA(seq), DNA(j))

    end = start_end[0][1]

    # Translate full CDR3 for stop codons
    # minus start of v with start of j, the results should be divisible with 3

    cdr3 = seq[start - 3:end + 1]

    if '_' in util.translate(cdr3):
        raise InternalStopCodonException

    if len(cdr3) % 3 != 0:
        raise FrameshiftException
    return cdr3
