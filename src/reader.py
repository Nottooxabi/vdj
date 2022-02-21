from Bio import SeqIO
from vj_finder import util
import pandas as pd


# TODO: Encapsulate functions into class, storing read sequences into dictionary



def read_ab1_file(file_name: str, handle_n: bool = False, naming: str = None,
                  chain: str = None, reverse_transcribe: bool = False):
    """
    Reads a .ab1 file returns either basecalls or basecalls and trace information.

    Args:
        file_name: Input filename of ab1 file to be read
        handle_n: bool, If true, will attempt to call unknown nucleotides, requires ab1 file
        naming: Naming nomenclature for the file name.
                If no naming convention is indicated will use whole name as name
        chain: TCR chain of the sequence, can be alpha, beta, gamma or delta,
        If no chain is indicated, will pass as NaN
        reverse_transcribe: bool, reverse transcribe sequence, Default is false

    Returns:
        Either a string sequence or dictionary for different attributes
    """

    if naming == 'agrf':
        # AGRF naming adds index to end of file name separated with an underscore
        name = file_name.split('_')[0]
    else:
        name = file_name
    read = SeqIO.read(file_name, 'abi')

    read_data = read.annotations["abif_raw"]
    call_index = read_data["PLOC2"]  # Index location of each called base
    primary_basecall = read_data['PBAS1'].decode('utf-8')  # Predicted sequence of file
    secondary_basecall = read_data['P2BA1'].decode('utf-8')

    # These channels contain trace data for each nucleotide
    data9 = read_data['DATA9']
    data10 = read_data['DATA10']
    data11 = read_data['DATA11']
    data12 = read_data['DATA12']

    key = read_data['FWO_1'].decode('utf-8')  # Maps correct nucleotide to each channel

    if handle_n:
        data = {'seq_id': name, 'sequence': primary_basecall, 'chain': chain,
                'secondary': secondary_basecall, 'call_index': call_index,
                key[0]: data9, key[1]: data10, key[2]: data11, key[3]: data12}

        primary_basecall = call_n_ab1(data)

    if reverse_transcribe:
        primary_basecall = util.reverse_complement(primary_basecall)

    return {'seq_id': name, 'sequence': primary_basecall, 'chain': chain}


# TODO: Find more test files
def call_n_ab1(ab1_data: dict):
    """
    Forces a nucleotide call for 'N' labeled nucleotides in ab1 sequence data.
    Trace signal from the initial 20 (arbitrary) or so nucleotides are highly inaccurate
    (outside range of detection) and therefore are not processed.
    Heterogeneous traces are not detected, this function may result in invalid sequences as a result.
    Double check files which contain a high number of poor quality base calls

    Args:
        ab1_data: dict, trace and sequence information, must be from read_ab1_file function

    Returns:
        Nucleotide string with unknown nucleotides resolved
    """

    # Run through primary call sequence stop at any 'N' characters encountered
    sequence = list(ab1_data['sequence'])
    peak_index = ab1_data['call_index']
    for i in range(len(sequence)):
        if (sequence[i] == 'N') and (i > 20):
            index_for_n = peak_index[i]
            values_at_peak = {'A': ab1_data['A'][index_for_n], 'T': ab1_data['T'][index_for_n + 1],
                              'G': ab1_data['G'][index_for_n + 1],
                              'C': ab1_data['C'][index_for_n + 1]}
            # sort based on which one have the highest trace value,
            highest = sorted(values_at_peak, key=values_at_peak.get, reverse=True)
            # For each nucleotide get index for adjacent peaks, get middle value
            # Compared if middle values are smaller or greater than peak value
            # If both are smaller, replace N with that value
            try:
                left_peak_index = peak_index[i - 1]
                right_peak_index = peak_index[i + 1]
            except IndexError:
                continue

            left_difference = round(index_for_n - ((index_for_n - left_peak_index) / 3))
            right_difference = round(index_for_n + ((right_peak_index - index_for_n) / 3))

            #
            for nucleotide in highest:
                centre = ab1_data[nucleotide][index_for_n]
                left = ab1_data[nucleotide][left_difference]
                right = ab1_data[nucleotide][right_difference]

                if centre < 50:  # Arbitrary value, anything below is assumed to be background
                    break  # TODO: Could also call nucleotide based on highest value if above a certain threshold
                elif centre > left and centre > right:
                    sequence[i] = nucleotide.lower()
                    break
                elif (left < centre < right) and centre > 100:
                    sequence[i] = nucleotide.lower()
                    break
                elif (left > centre > right) and centre > 100:
                    sequence[i] = nucleotide.lower()
                    break

    return ''.join(sequence)


def read_fasta_file(file_name, chain: str):
    """

    Args:
        file_name: str, path to fasta file
        chain: str, can be A, B, D, G

    Returns: list, [{'name': file_name, 'seq': AT, chain: 'A'}, ... ]
    """

    accepted_chain = ['A', 'B', 'G', 'D']

    if chain.upper() not in accepted_chain:
        raise ValueError("Invalid chain type")

    with open(file_name, 'r') as file:
        raw_text = file.read()
        raw_text = raw_text.split('>')
        if len(raw_text[0]) == 0:
            raw_text = raw_text[1:]

    all_genes = []

    if file_type == 'clc':
        for i in raw_text:
            sequence = fasta_export(i, chain)
            all_genes.append(sequence)

    return all_genes


def fasta_export(fasta_section, chain: str):
    """

    Args:
        fasta_section:
        chain:

    Returns:

    """
    fasta_section = fasta_section.split('\n')
    attributes = fasta_section[0]
    seq = ''.join(fasta_section[1:])

    attributes = attributes.split(';')
    name = attributes[0].split(" ")[0]

    ## Shit naming convention Michael
    if chain == 'A':
        name = name.replace('a', '')
    elif chain == 'B':
        name = name.replace('b', '')
    length = len(seq)
    ambiguous = seq.count('N')

    return {'name': name, 'chain': chain, 'seq': seq}


# TODO script for reading seq file


def mixcr_export(file, sep='\t', index=None):
    """
    Parses alignment file from mixcr exportAlignment output
    By default assumes that the inputted file in in a tsv format
    """

    return pd.read_csv(file, sep=sep, index_col=index)


def parse_row(sequence, v_alignments, j_alignments):
    """
    :param sequence: sequences of inputted gene
    :param v_alignments: string containing all variable segments the sequence was matched to
    :param j_alignments: string containing all joining segments the sequences was matched to

    :return: a dictionary with the following structure: {'sequence': 'ATGC', 'v_alignments': ['TRAV1': 12], 'j_alignments': ['TRAJ12']: 11}
    """

    temp = {}

    temp['sequence'] = sequence
    temp['v_alignments'] = get_alignments(v_alignments)
    temp['j_alignments'] = get_alignments(j_alignments)

    return temp


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
        score = float(re.search(r'\(([^)]+)\)', i).group(0)[1:-2])  # Get score for each alignment
        temp[gene[:-3]] = score
        alignments.append(temp)
    return alignments


def get_mixcr_align_scores(df):
    """
    Args:
        df: pandas Dataframe containing mixcr alignment export data

    Returns:
        list of parsed scores
    """

    return [parse_row(x, y, z) for x, y, z in
            zip(df['targetSequences'], df['allVHitsWithScore'], df['allJHitsWithScore'])]
