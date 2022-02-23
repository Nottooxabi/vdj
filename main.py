from src import ab1_util, imgt_reference, management, alignment
import os

if __name__ == '__main__':

    # This needs to be done once
    #imgt_reference.get_imgt_references('references/imgt_ref.fasta')

    ref = management.ReferenceManager('references/imgt_ref.fasta')
    ref.set_current(species='Mus musculus', segment=['V'])



    pass

