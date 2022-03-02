from Bio import pairwise2
from Bio.pairwise2 import format_alignment

from src import ab1_util, imgt_reference, management, alignment, util, user
import os
import pandas as pd

if __name__ == '__main__':
    # This needs to be done once
    # imgt_reference.get_imgt_references('references/imgt_ref.fasta')

    ref = management.ReferenceManager('references/imgt_ref.fasta')
    ref.set_current(species='Mus musculus', segment=['V', 'J'], chain=['G', 'D'])

    data = user.UserSequences(ref, species='Mus musculus')
    data.set_files('data/')
    meta = {'data/RS_54': {'chain': 'G', 'isReverse': True},
            'data/RS_55': {'chain': 'D', 'isReverse': True}}

    data.set_meta_data(meta)

    data.read_files()

    df = data._pair_chains('data/RS_54','data/RS_55', ('_g', '_d'))
    clones = df.groupby(['cdr3_g', 'cdr3_d']).size().reset_index()

    df.to_csv('temp.csv')

    pass
