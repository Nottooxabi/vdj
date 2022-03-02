"""
This script manages vj reference sequences, the ReferenceManagement class stores all
information as a pandas dataframe as well as a dictionary of gene objects
"""
import pandas as pd
import pickle
import os
from src.util import translate
import src.imgt_reference as ref
from src.exceptions import InvalidEntryException


def load_data():
    reference_storage = open('reference', 'rb')
    management_object = pickle.load(reference_storage)
    reference_storage.close()

    return management_object


class ReferenceManager:
    """
    Object maintains all reference vdj sequences as a pandas dataframe
    functions allow splicing of dataframe based on species, chain and segment etc which is
    stroed in the self.current attribute
    """

    def __init__(self, reference: str = " "):
        if os.path.isfile(reference):
            self.references_df = pd.DataFrame.from_records(
                ref.create_genes(ref.parse_fasta(reference)))

            # This attribute stores a splice of the reference genome
            # which will be referred to for alignments
            self.current = None
        else:
            raise FileExistsError("Inputted file does not exist")

    def set_current(self,
                    species: str,
                    isTCR: bool = True,
                    variety: list = None,
                    segment: list = None,
                    chain: list = None,
                    allele: list = ['01']):
        """
        Sets the current dataframe which will be used for sequence alignments, requires both a valid
        species entry to work
        """

        self.current = self.references_df[(self.references_df['species'] == species)]

        self._check_current('species')

        if isTCR:
            self.current = self.current[self.current['allele'].str.contains('TR')]
        else:
            self.current = self.current[self.current['allele'].str.contains('IG')]

        self._check_current('allele')

        if variety is not None:
            self.current = self.current[self.current['variety'].isin(variety)]
            self._check_current('variety')

        if segment is not None:
            self.current = self.current[self.current['segment'].isin(segment)]
            self._check_current('segment')

        if chain is not None:
            self.current = self.current[self.current['chain'].isin(chain)]
            self._check_current('chain')

        if allele is not None:
            self.current[self.current['sub'].isin(allele)]
            self._check_current('sub')

        self.current['gene'] = self.current['allele'].str.split('*').str[0]
        # clean DV
        self.current['temp'] = self.current['gene'].apply(
            lambda x: x.split('/')[1] if len(x.split('/')) > 1 else "")
        self.current['temp'] = self.current['temp'].str.replace('D', '')
        self.current['temp'] = self.current['temp'].apply(lambda x: 'D' + str(x) if x != '' else x)

        self.current['gene'] = self.current['gene'].str.split('/').str[0]
        self.current['gene'] = self.current['gene'].str.replace('D', '', 1).str.replace('N', '', 1)

    def _check_current(self, column):

        """
        Used to check if the self.current database contains any entries, if the length is 0
        self.current will be reset to zero and a error is raised
        A column value is asked so that correct inputs can be printed
        Args:
            column: column name or error handling
        """

        if len(self.current) == 0:
            self.current = None

            possible = list(self.references_df[column].unique())

            raise ReferenceError("Inputted " +
                                 str(column) +
                                 " value does not exist \n Only the following inputs are valid: \n" +
                                 str(possible))

    def to_align(self, chain: list, segment: list, regions: str):

        """
        Pull all sequences to align against and which regions to use for alignment
        Args:
            chain: list of A, B, G, D
            segment: list containing V, D or J
            regions: str which indicates which regions to align against

        Returns: dict of a dict
        """

        # Set contains all possible iterations of fr and cdr sections which can be used, combinations of sections are
        # formed by csv format

        accepted = {'fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3', 'fr4', 'fr1,cdr1', 'cdr1,fr2', 'fr2,cdr2', 'cdr2,fr3',
                    'fr3,cdr3', 'cdr3,fr4', 'fr1,cdr1,fr2', 'cdr1,fr2,cdr2', 'fr2,cdr2,fr3', 'cdr2,fr3,cdr3',
                    'fr3,cdr3,fr4', 'fr1,cdr1,fr2,cdr2', 'cdr1,fr2,cdr2,fr3', 'fr2,cdr2,fr3,cdr3', 'cdr2,fr3,cdr3,fr4',
                    'fr1,cdr1,fr2,cdr2,fr3', 'cdr1,fr2,cdr2,fr3,cdr3', 'fr2,cdr2,fr3,cdr3,fr4',
                    'fr1,cdr1,fr2,cdr2,fr3,cdr3', 'cdr1,fr2,cdr2,fr3,cdr3,fr4',
                    'fr1,cdr1,fr2,cdr2,fr3,cdr3,fr4'}

        if regions not in accepted:
            raise InvalidEntryException

        temp_df = self.current[(self.current['chain'].isin(chain)) &
                               (self.current['segment'].isin(segment))]

        regions = regions.split(',')
        to_align = {}
        for index, row in temp_df.iterrows():
            comb = ''
            gene = row.allele

            for region in regions:
                comb = comb + row[region].replace('.', '')
            to_align[gene] = comb

        return to_align

    def get_chain_sequences(self, usage):
        """
        Get cdr3 sequence from inputted V and J genes
        Args:
            usage: list or tuple of two elements the first being the V gene and second being a J gene

        Returns: tuple of (v, j) sequences

        """
        usage = list(usage)
        v = usage[0]
        j = usage[1]

        cv = self.current[self.current['allele'] == v].cdr3.iloc[0]
        cj = self.current[self.current['allele'] == j].cdr3.iloc[0]

        return cv, cj

    def save_data(self):
        try:
            os.remove('reference')
        except OSError:
            pass

        reference_storage = open('reference', 'ab')
        pickle.dump(self, reference_storage)
        reference_storage.close()
        return
