"""
The user module is designed to read and manipulate sequences which are inputted by a user
These functions primarily involve inputted paths to specific files

Currently, the following formats can be parsed: ab1, fasta and mixcr

"""
import os
import pandas as pd

from src import ab1_util, management, util
from vj_finder import exceptions



def get_files(path=None):
    """
    Pulls relevant ab1 files from selected path

    Args:
        path: Optional, path to input files

    Returns:
        Dictionary with ab1 path string in the following format:
        {'Directory': ['Directory/file.ab1']}

    """

    files_in_data = os.listdir('Data')
    all_files = {}
    for file in files_in_data:
        path = 'Data/' + file
        try:
            files_in_path = os.listdir(path)
        except NotADirectoryError:
            pass
        accepted_files = []
        for i in files_in_path:
            try:
                extension = i.split('.')
                if extension[1] == 'ab1':
                    accepted_files.append(i)
            except IndexError:
                pass

        all_files[path] = accepted_files

    return all_files


def read_file(file_type: str, file_name: str):
    """
    Reads specified type of sequence file and appends file into read_sequences list
    Args:
        file_type: file extension type for example (fasta)
        file_name: path to file
    Returns:
    Sequence from file
    """

    if file_type == 'ab1':
        try:
            sequence = reader.read_ab1_file(file_name, handle_n=True, reverse_transcribe=True)
            return sequence
        except OSError:
            return "Needs new solution"
    elif file_name == 'fasta':
        return "Not implemented"
    else:
        raise ValueError("Invalid file type")

class UserSequences:
    """
    User inputted ab1 files are parsed and stored in this object
    Is paired with a reference class which stores reference sequecnes

    Attributes:
    __read_sequences -> list of sequences which have read and aligned against reference
    references -> reference object
    isTCR -> boolean indicating is supplied files contain TCR sequences or not
    species -> which species
    """

    def __init__(self,
                 references: management.ReferenceManager = None,
                 is_tcr: bool = True,
                 species: str = '',
                 variety: str = ''):
        self.__read_sequences = []
        self.references = references
        self.is_tcr = is_tcr
        self.species = species
        self.variety = variety

    def set_references(self):
        """
        Changes attributes in reference instance

        Args:
            references: Reference Object
        """
        self.references.set_species(species=self.species, variety=self.variety, isTCR=self.is_tcr)

    def get_targets(self, chain: str, segment: str, locations: list):

        """
        Returns a list of targets for which each user sequence will be aligned against can be more than one section
        TODO: allow for FR1-4 and combinations of FR and CDR
        Args:
            chain:
            segment:
            locations:

        Returns:

        """

        allowed = set(['cdr1', 'cdr2', 'cdr3'])

        references = []

        for location in locations:
            if location in allowed:
                references.append(self.references.get_reference_sequences(chain=chain,
                                                                          segment=segment,
                                                                          cdr=location))
            else:
                print(allowed)
                continue

        return references

    def align_sequence(self, seq, chain, cdr=['cdr1', 'cdr2', 'cdr3']):

        """
        Aligns a tcr chain against references and finds best matching sequence and cdr3
        Args:
            cdr: which regions to align against
            chain: TCR chain

        Returns: list containing best vj pairing

        """
        error = False

        v_ref = self.get_targets(chain=chain, segment='V', locations=cdr)
        best_v = set()
        if len(v_ref) == 0:
            return 0

        for ref in v_ref:
            best_v.add(vj.find_best_allele(seq, ref))

        # TODO: Solve issue if more than one v is found when multiple locations are used for
        #  determination of allele
        best_v = best_v.pop()

        j_ref = self.get_targets(chain=chain, segment='J', locations=['cdr3'])
        best_j = vj.find_best_allele(seq, j_ref[0])

        try:
            v = self.references.current[self.references.current['allele'] == best_v].iloc[0]['cdr3']
            j = self.references.current[self.references.current['allele'] == best_j].iloc[0]['cdr3']
            cdr3 = vj.define_cdr3(seq, v, j)
        except exceptions.FrameshiftException:
            cdr3 = 'Frameshift'
            error = True
        except exceptions.InternalStopCodonException:
            cdr3 = 'Stop Codon'
            error = True
        except ValueError:
            cdr3 = '?'
            error = True

        return {'Variable': best_v,
                'Joining': best_j,
                'nuc_cdr3': cdr3,
                'prot_cdr3': util.translate(cdr3),
                'errors': error}

    def set_all_usage(self, chain):
        """
        Runs set_usage_information method for each sequence

        Args:
            chain: sequence chain

        Returns: Nothing

        """

        for i in self.__read_sequences:
            print(i)
            print(self.align_sequence(i['sequence'], chain))

    # TODO improve this function to allow custom directories
    def read_all_files(self):
        """
        Real each file based on files in directory

        Returns: Nothing

        """
        directory = get_files()
        root = os.getcwd()

        for path, names in directory.items():
            try:
                os.chdir(path)
                for j in names:
                    sequence = read_file(file_type='ab1', file_name=j)
                    if sequence is not None:
                        sequence['directory'] = path.split('/')[-1]
                        sequence['coordinate'] = sequence['seq_id'][0:2]
                        self.__read_sequences.append(sequence)
                os.chdir(root)

            except NotADirectoryError:
                pass

    def get_seq(self):
        """

        Returns: list of sequences stored in instance of object

        """
        return self.__read_sequences

    def get_clones(self, chain_1: list, chain_2: list):

        """
        Pairs two chains from the same cell together
        Requires both lists to be of the same length

        Args:
            chain_1: list containing directories for one chain
            chain_2: list containing directories to corresponding chain
        """

        if len(chain_1) != len(chain_2):
            raise IndexError('equal number of inputs required')

        # Current implementation has no error handling
        # Poor inputs will result is incorrect pairings
        dataframe = pd.DataFrame.from_records(self.__read_sequences)

        for i in range(len(chain_1)):
            c_1 = chain_1[i]
            c_2 = chain_2[i]

            df_c1 = dataframe[dataframe['directory'] == c_1]
            df_c2 = dataframe[dataframe['directory'] == c_2]

            df_c1 = df_c1[
                'seq_id', 'coordinate', 'chain', 'directory',
                'v_gene', 'j_gene', 'cdr3_seq', 'cdr_aa', 'errors']

            df_c2 = df_c2[
                'seq_id', 'coordinate', 'chain', 'directory',
                'v_gene', 'j_gene', 'cdr3_seq', 'cdr_aa', 'errors']

        return dataframe.groupby(['cdr3_seq', 'v_gene', 'j_gene'])

    def convert_to_csv(self):
        """
        Exports data into csv file, user inputs file name
        Returns: Nothing

        """
        print("Enter file name: ", end=" ")
        file_name = input()
        if '.csv' not in file_name:
            file_name = file_name.split('.')[0]
            file_name = file_name + '.csv'

        path = os.getcwd() + "/Output/" + file_name
        try:
            pd.DataFrame.from_records(self.__read_sequences).to_csv(path)
        except PermissionError:
            pd.DataFrame.from_records(self.__read_sequences).to_csv(file_name)

    def debug(self, seq_id):
        error = []
        for i in self.__read_sequences:
            if i['seq_id'] == seq_id:
                error.append(i)
                break
        self.__read_sequences = error
