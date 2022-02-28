from Bio import pairwise2
from Bio.pairwise2 import format_alignment

from src import ab1_util, imgt_reference, management, alignment
import os

if __name__ == '__main__':

    # This needs to be done once
    # imgt_reference.get_imgt_references('references/imgt_ref.fasta')

    ref = management.ReferenceManager('references/imgt_ref.fasta')
    ref.set_current(species='Homo sapiens', segment=['V', 'J'])
    x = ref.to_align(chain=['B'], segment=['J'], regions='cdr3')

    test = 'TCCTGCAGCCAGAAGACTCGGCCCTGTATCTCTGCGCCAGCAGCCAAGAAGGGTCGGACCTCGTACTTCGCGAGCAGTACTTCGGGCCGGGCACCAGGCTCACGGTCACAGAGGACNNNNAANNNNNNNNNNNNNN'

    v = alignment.align(test, x)
    pass