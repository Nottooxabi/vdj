"""Collection of functions which allow for manipulation of dna sequences.

    -> translate: translates inputted DNA sequence into amino acids based on frame
    -> reverse complement: creates reverse complement of a given sequence

"""


def reverse_complement(seq: str):
    """
    Reverse complements inputted DNA sequence
    Args:
        seq: str, DNA sequence

    Returns: str, Reverse complement of inputted sequence
    """

    sequence = seq.upper()
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    complement_sequence = []
    for i in sequence:
        try:
            complement_sequence.append(complement[i])
        except KeyError as error:
            print('Unexpected value encountered, replaced ' + str(error) + ' with ' + '_')
            complement_sequence.append('_')

    return ''.join(complement_sequence[::-1])


def translate(seq: str, frame: int = 1, debug = False):
    """
    Translation of a DNA sequence into a Amino Acid Sequence

    Args:
        seq: str, DNA sequence
        frame: int, Reading frame. Default is 1

    Returns:
            str, Translated protein sequence
            Unknown values will be treated as X, stop codons terminate translation

    """

    # prevents errors associated by frame value being either negative or a float
    frame = round(abs(frame))

    # prevents errors caused by incorrect case and whitespace
    seq = seq.upper().replace(' ', '')
    protein_sequence = []
    translation_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    triplets = [seq[i:i + 3] for i in range((frame - 1), len(seq), 3)]

    for i in triplets:
        try:
            if len(i) != 3:  # Non-triplets which should only appear at the end are ignored
                if debug:
                    protein_sequence.append(str(len(i)))
                    continue
                else:
                    continue
            protein_sequence.append(translation_table[i])
        except KeyError:  # Used to handle Unknown triplets
            protein_sequence.append('X')

    return ''.join(protein_sequence)


def get_groups(df):
    alleles = df.current.groupby('gene')
    all_alleles = {}
    for allele, values in alleles:
        all_alleles[allele] = []
        for index, row in values.iterrows():
            temp = {row['allele']: row['seq']}
            all_alleles[allele].append(temp)
    return all_alleles


def export_fasta(genes: dict):
    for key, values in genes.items():
        total = ''
        file_name = 'reference_fasta/' + str(key) + '.fasta'
        for value in values:
            for allele, seq in value.items():
                temp_string = '>' + str(allele) + '\n' + str(seq) + '\n'
                total += temp_string

        file = open(file_name, 'w')
        file.write(total)
        file.close()
