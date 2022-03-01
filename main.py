from Bio import pairwise2
from Bio.pairwise2 import format_alignment

from src import ab1_util, imgt_reference, management, alignment, util
import os

if __name__ == '__main__':
    # This needs to be done once
    # imgt_reference.get_imgt_references('references/imgt_ref.fasta')

    ref = management.ReferenceManager('references/imgt_ref.fasta')
    ref.set_current(species='Homo sapiens', segment=['V', 'J'])

    x = ref.get_chain_sequences(['TRAV1-2', 'TRAJ33'])

    #test = 'TCCTGCAGCCAGAAGACTCGGCCCTGTATCTCTGCGCCAGCAGCCAAGAAGGGTCGGACCTCGTACTTCGCGAGCAGTACTTCGGGCCGGGCACCAGGCTCACGGTCACAGAGGACNNNNAANNNNNNNNNNNNNN'
    #x = ref.to_align(chain=['B'], segment=['V'], regions='fr3,cdr3')
    #v = alignment.align(test, x)

    #x = ref.to_align(chain=['B'], segment=['J'], regions='cdr3,fr4')
    #j = alignment.align(test, x)


    #get_cdr3
    #x = ref.to_align(chain=['B'], segment=['V', 'J'], regions='cdr3')

    #z = alignment.define_cdr3(test, x[v], x[j])
    #print(util.translate(z))
    #pass
