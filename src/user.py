"""
The user module is designed to read and manipulate sequences which are inputted by a user
These functions primarily involve inputted paths to specific files

Currently, the following formats can be parsed: ab1, fasta and mixcr

"""
import os
import pandas as pd

from src import ab1_util, management, util, alignment
from src import exceptions


class UserSequences:
    """
    User inputted ab1 files are parsed and stored in this object
    Is paired with a reference class which stores reference sequecnes

    Attributes:
    __read_sequences -> list of sequences which have read and aligned against reference -> reference object
    isTCR -> boolean indicating is supplied files contain TCR sequences or not
    species -> which species
    """

    def __init__(self,
                 references: management.ReferenceManager = None,
                 is_tcr: bool = True,
                 species: str = '',
                 variety: str = ''):

        self.references = references
        self.is_tcr = is_tcr
        self.species = species
        self.variety = variety
        self.__processed_files = {'plate': [],'sample': [], 'sequence': [], 'V': [], 'J': [], 'cdr3': [], 'cdr3_protein': [],
                                  'errors': []}
        self.__read_sequences = {}
        self.__meta = None

    def set_references(self):
        """
        Changes attributes in reference instance

        Args:
            references: Reference Object
        """
        self.references.set_species(species=self.species, variety=self.variety, isTCR=self.is_tcr)

    def set_files(self, path):
        """
        Pulls relevant ab1 files from selected path

        Args:
            path: path to input files

        Returns:
            Dictionary with ab1 path string in the following format:
            {'Directory': ['Directory/file.ab1']}

        """

        directories = os.listdir(path)

        invalid_files = []
        valid = {}
        for directory in directories:

            temp_path = path + directory
            try:
                valid[temp_path] = os.listdir(temp_path)
            except NotADirectoryError:
                if '.ab1' in directory:
                    invalid_files.append(directory)

        if len(invalid_files) != 0:
            print('Following files were not read, please place into directory: ' + '. '.join(invalid_files))

        self.__read_sequences = valid

        print(self.__read_sequences.keys())

    def set_meta_data(self, meta: dict):
        """
        Specify meta data for each set of sequences such as chain etc. as a list
        some meta info can be used to find
        Args:
            meta: key is path to file, needs to be the same as keys from self.__read_sequences, values are a dictionary
            each dictionary must contain the key (chain) specifying the tcr chain of the files
            in addition isReverse needed for reverse sequenmces

        Returns:

        """
        self.__meta = None

        if self.__read_sequences is None:
            print('no sequence files found')
            return
        if meta.keys() != self.__read_sequences.keys():
            print('meta keys dont match read sequence keys')
            return

        for k, v in meta.items():
            if 'chain' not in v:
                print('all meta entries requires {\'chain\': value}')
                return

            if 'isReverse' not in v:
                meta['isReverse'] = False

        self.__meta = meta

    def read_files(self, v_loc='fr3,cdr3', j_loc='cdr3,fr4'):

        if self.__meta is None:
            raise ValueError('No meta data found')

        for key, value in self.__meta.items():
            v_ref = self._get_targets(chain=value['chain'], segment=['V'], locations=v_loc)
            j_ref = self._get_targets(chain=value['chain'], segment=['J'], locations=j_loc)

            for entry in self.__read_sequences[key]:
                self.__processed_files['plate'].append(key)
                self.__processed_files['sample'].append(entry)
                path = key + '/' + entry
                try:
                    seq = ab1_util.read_ab1_file(path, handle_n=True, reverse_transcribe=value['isReverse'])
                    self.__processed_files['sequence'].append(seq['sequence'])
                    temp = self._get_usages(seq['sequence'], v_ref, j_ref)
                    for k, v in temp.items():
                        self.__processed_files[k].append(v)

                except FileNotFoundError:
                    print('file not found: ' + path)

    # TODO: Needed?
    def _get_targets(self, chain: list, segment: str, locations: str):

        """
        Returns a list of targets for which each user sequence will be aligned against can be more than one section
        Args:
            chain:
            segment:

        Returns:

        """

        try:
            references = self.references.to_align([chain], segment, locations)
        except exceptions.InvalidEntryException:
            return {}
        return references

    def _get_usages(self, seq, for_v, for_j):

        """
        Aligns a tcr chain against references and finds best alleles for each sequence
        Args:
            seq: query sequence which references will be aligned against
            chain: TCR chain A, B ,G ,D
            for_v: which v segments will be used for references
            for_j:
        Returns: list containing best vj pairing

        """

        data = {'V': alignment.align(seq, for_v), 'J': alignment.align(seq, for_j)}
        usage = self.references.get_chain_sequences(data.values())
        cdr3, errors = alignment.define_cdr3(seq, usage[0], usage[1])

        data['cdr3'] = cdr3
        data['cdr3_protein'] = util.translate(cdr3)
        data['errors'] = errors

        return data

    def get_seq(self):
        """

        Returns: list of sequences stored in instance of object

        """
        return self.__read_sequences

    def _to_table(self):
        """
        Converts data into a pandas dataframe
        Returns:

        """

        return pd.DataFrame(self.__processed_files)

    def _pair_chains(self, plate_1: str, plate_2: str, plate_suffices: tuple):

        """
        Pairs two chains together with
        This requires two 'plates' with paired TCRs having the same sample name

        Args:
            plate_1: list containing directories for one chain
            plate_2: list containing directories to corresponding chain
            plate_suffices: suffices used for column names

        """

        if len(self.__processed_files['plate']) == 0:
            raise ValueError("no plates have beeen processed yet")

        df = self._to_table()

        samples_1 = df[df['plate'] == plate_1]
        samples_2 = df[df['plate'] == plate_2]

        return pd.merge(samples_1, samples_2, on='sample', suffixes=plate_suffices)



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
            sequence = ab1_util.read_ab1_file(file_name, handle_n=True, reverse_transcribe=True)
            return sequence
        except OSError:
            return "Needs new solution"
    elif file_name == 'fasta':
        return "Not implemented"
    else:
        raise ValueError("Invalid file type")


