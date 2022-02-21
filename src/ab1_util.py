import util
from Bio import SeqIO


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

            for nucleotide in highest:
                centre = ab1_data[nucleotide][index_for_n]
                left = ab1_data[nucleotide][left_difference]
                right = ab1_data[nucleotide][right_difference]

                if centre < 50:  # Arbitrary value, anything below is assumed to be background
                    break
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
