from Bio.pairwise2 import format_alignment

from src import util
from skbio.alignment import StripedSmithWaterman
from skbio.alignment import local_pairwise_align_ssw, global_pairwise_align_nucleotide
from skbio.sequence import DNA
from src.exceptions import FrameshiftException, InternalStopCodonException
from Bio import pairwise2

"""
Identifies correct gene usage for inputted sequences
"""


def align(query: str, ref: dict, debug: bool=False):
    best_gene = ''
    best_score = 0
    for gene, seq in ref.items():
        # Severe penalty for gap opens and mismatches
        alignment_score = pairwise2.align.localms(query, seq, match=2, mismatch=-10, open=-20, extend=-0.1,
                                                  one_alignment_only=True, score_only=True)

        if alignment_score > best_score:
            best_gene = gene
            best_score = alignment_score

            if debug:
                for a in  pairwise2.align.localms(query, seq, match=2, mismatch=-10, open=-20, extend=-0.1,
                                                      one_alignment_only=True):
                    print(gene)
                    print(format_alignment(*a))

    return best_gene


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
