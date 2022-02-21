"""
This script manages vj reference sequences, the ReferenceManagement class stores all
information as a pandas dataframe as well as a dictionary of gene objects
"""
import pandas as pd
import pickle
import os
from src.util import translate
import src.imgt_reference as ref


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

            # This attribute stores a splice of the reference genome which will be referred to for alignments
            self.current = None
        else:
            raise FileExistsError("Inputted file does not exist")

    def get_reference_sequences(self, chain: list = None, segment: list = ['V'], func: list = ['F'], cdr: list = None):
        """
        Get all sequences from refererence based on conditions
        Args:
            chain: tcr chain
            segment:
            func:
            cdr:

        Returns:

        """

        if chain is None:
            data = self.current[(self.current['segment'] == segment)
                                & (self.current['functional'] == func)]

            return pd.Series(data['seq'].values, index=data['allele']).to_dict()

        data = self.current[(self.current['segment'] == segment)
                            & (self.current['functional'] == func)
                            & (self.current['chain'] == chain)]
        if cdr is None:
            return pd.Series(data['seq'].values, index=data['allele']).to_dict()
        return pd.Series(data[cdr].values, index=data['allele']).to_dict()

    def set_species(self, species, isTCR=True, variety=None, segment=None, chain=None,
                    allele='\*01'):
        """
        Sets the current dataframe which will be used for sequence alignments, requires both a valid
        species entry to work
        """

        self.current = self.references_df[(self.references_df['species'] == species)]

        self.check_current('species')

        if isTCR:
            self.current = self.current[self.current['allele'].str.contains('TR')]
        else:
            self.current = self.current[self.current['allele'].str.contains('IG')]

        self.check_current('allele')

        if variety is not None:
            self.current = self.current[self.current['variety'] == variety]
            self.check_current('variety')

        if segment is not None:
            self.current = self.current[self.current['segment'] == segment]
            self.check_current('segment')

        if chain is not None:
            self.current = self.current[self.current['chain'] == chain]
            self.check_current('chain')

        if allele is not None:
            self.current[self.current['allele'].str.contains(allele)]
            self.check_current('allele')

        self.current['gene'] = self.current['allele'].str.split('*').str[0]
        # clean DV
        self.current['temp'] = self.current['gene'].apply(
            lambda x: x.split('/')[1] if len(x.split('/')) > 1 else "")
        self.current['temp'] = self.current['temp'].str.replace('D', '')
        self.current['temp'] = self.current['temp'].apply(lambda x: 'D' + str(x) if x != '' else x)

        self.current['gene'] = self.current['gene'].str.split('/').str[0]
        self.current['gene'] = self.current['gene'].str.replace('D', '', 1).str.replace('N', '', 1)

        # self.current['temp'] = self.current['allele'].str.split('*').str[0].str.split['/']
        # self.current['type'] = self.current['allele'].str.split('*').str[1]

    def check_current(self, column):
        if len(self.current) == 0:
            self.current = None

            possible = list(self.references_df[column].unique())

            raise ReferenceError("Inputted " +
                                 str(column) +
                                 " value does not exist \n Only the following inputs are valid: \n" +
                                 str(possible))

    def get_chain_sequences(self, usage):
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
