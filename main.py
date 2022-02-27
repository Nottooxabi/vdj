from src import ab1_util, imgt_reference, management, alignment
import os

if __name__ == '__main__':

    # This needs to be done once
    #imgt_reference.get_imgt_references('references/imgt_ref.fasta')

    ref = management.ReferenceManager('references/imgt_ref.fasta')
    ref.set_current(species='Homo sapiens', segment=['V'])
    x = ref.to_align(chain=['B'], segment=['V'], regions=['fr3'])
    test = 'TCCTGCAGCCAGAAGACTCGGCCCTGTATCTCTGCGCCAGCAGCCAAGAAGGGTCGGACCTCGTACTTCGCGAGCAGTACTTCGGGCCGGGCACCAGGCTCACGGTCACAGAGGACNNNNAANNNNNNNNNNNNNN'

    for k, v in x.items():
        temp_best = alignment.new_align(test, v)

        print(k + ':' + temp_best )

    pass

