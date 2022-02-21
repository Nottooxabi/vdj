import re

from Bio import SeqIO
from src import util
import pandas as pd


class SequenceReader:

    def __init__(self):
        self.processed_files = {}





# TODO remake function
def read_fasta_file(file_name, chain: str):
    """

    Args:
        file_name: str, path to fasta file
        chain: str, can be A, B, D, G

    Returns: list, [{'name': file_name, 'seq': AT, chain: 'A'}, ... ]
    """

    return


class MixcrData:

    def __init__(self, file, sep):
        self.data = pd.read_csv(file, sep=sep)

    def parse_row(self, sequence, v_alignments, j_alignments):
        """
        :param sequence: sequences of inputted gene
        :param v_alignments: string containing all variable segments the sequence was matched to
        :param j_alignments: string containing all joining segments the sequences was matched to

        :return: a dictionary with the following structure: {'sequence': 'ATGC', 'v_alignments': ['TRAV1': 12], 'j_alignments': ['TRAJ12']: 11}
        """

        temp = {'sequence': sequence, 'v_alignments': self.get_alignments(v_alignments),
                'j_alignments': self.get_alignments(j_alignments)}

        return temp

    @staticmethod
    def get_alignments(matches: str):
        """
        Args:
            matches:

        Returns:

        """

        alignments = []
        matches = matches.split(',')
        for i in matches:
            temp = {}
            gene = i.split('(')[0]
            score = float(
                re.search(r'\(([^)]+)\)', i).group(0)[1:-2])  # Get score for each alignment
            temp[gene[:-3]] = score
            alignments.append(temp)
        return alignments

    def get_mixcr_align_scores(self):
        """
        Returns:
            list of parsed scores
        """

        return [self.parse_row(x, y, z) for x, y, z in
                zip(self.data['targetSequences'],
                    self.data['allVHitsWithScore'],
                    self.data['allJHitsWithScore'])]
