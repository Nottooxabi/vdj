from Bio.pairwise2 import format_alignment

from src import util
from skbio.alignment import StripedSmithWaterman
from skbio.alignment import local_pairwise_align_ssw, global_pairwise_align_nucleotide
from skbio.sequence import DNA
from src.exceptions import FrameshiftException, InternalStopCodonException, Cys104NotFoundException
from Bio import pairwise2

"""
Identifies correct gene usage for inputted sequences
"""


def align(query: str, ref: dict, debug: bool = False):
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
                for a in pairwise2.align.localms(query, seq, match=2, mismatch=-10, open=-20, extend=-0.1,
                                                 one_alignment_only=True):
                    print(gene)
                    print(format_alignment(*a))

    return best_gene


def show_alignment(seq, seq_v, seq_j):
    """
    For showing alignment of query sequence with another
    Args:
        seq:
        seq_v:
        seq_j:

    Returns:
        formatted alignment
    """

    v_alignment = pairwise2.align.localms(seq, seq_v, match=2, mismatch=-10,
                                          open=-20, extend=-0.1, one_alignment_only=True)

    j_alignment = pairwise2.align.localms(seq, seq_j, match=2, mismatch=-10,
                                          open=-20, extend=-0.1, one_alignment_only=True)

    v = v_alignment[0].seqB
    j = j_alignment[0].seqB

    temp = []
    for i in range(len(v)):

        if (v[i] == '-') and (j[i] == '-'):
            temp.append('-')
        elif (v[i] == '-') and (j[i] != '-'):
            temp.append(j[i])
        elif (v[i] != '-') and (j[i] == '-'):
            temp.append(v[i])
        else:
            temp.append('.')

    return seq, ''.join(temp)

def define_cdr3(seq: str, v: str, j: str):
    """
    Args:
        seq: user sequence for alignment
        v: variable cdr3 sequence
        j: joining segment cdr3 sequence

    Returns: cdr3 sequence of user sequence

    """
    seq = DNA(seq.upper())
    v = DNA(v)
    j = DNA(j)
    try:
        alignment, score, v_start_end = local_pairwise_align_ssw(seq, v, gap_open_penalty=1)
        start = v_start_end[0][0]

        alignment, score, j_start_end = local_pairwise_align_ssw(seq, j, gap_open_penalty=10)
        end = j_start_end[0][1]

        # Translate full CDR3 for stop codons
        # minus start of v with start of j, the results should be divisible with 3

        cdr3 = str(seq[start - 3: end + 1])

        if len(cdr3) == 0:
            return '', 'No Identifiable cdr3'

        # Internal Stop codon means unproductive tcr
        if '_' in util.translate(cdr3):
            return str(cdr3), 'Stop codon in CDR3'

        # Not not divisible by three suggests frameshift
        if len(cdr3) % 3 != 0:
            return str(cdr3), 'Frameshift'

        # Must begin with a C
        if util.translate(cdr3)[0] != 'C':
            return str(cdr3), 'No Cys at 104'

        return str(cdr3), ''

    except ValueError:
        return '', ''
