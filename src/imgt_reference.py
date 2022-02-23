import urllib.request
from src.util import translate


def get_imgt_references(file_name='imgt_ref.fasta',
                        url='http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences'
                            '.fasta-nt-WithoutGaps-F+ORF+inframeP'):
    """
    pulls all records from imgt databse
    Args:
        file_name: filename of a textfile to store data
        url: imgt database url

    Returns: filename

    """

    file = open(file_name, 'w')
    file.write(urllib.request.urlopen(url).read().decode('utf-8'))
    file.close()

    return file_name


def parse_fasta(file_name):
    """
    :param file_name: location of file in fasta format
    :return: list of gene entries
    """
    with open(file_name, 'r') as file:
        raw_text = file.read()
        raw_text = raw_text.split('>')
        if len(raw_text[0]) == 0:
            raw_text = raw_text[1:]  # First entry contains nothing
    return raw_text


def get_imgt_coordinate(value: str):
    """
    gets coordinate from inputted string
    can be either int..int or str,int..int
    """

    if ',' in value:
        for i in value.split(','):
            if '..' in i:
                value = i
                break
    coor = value.split('..')
    return [int(coor[0]), int(coor[1])]


def create_genes(raw_file):
    """
    :param raw_file:
    FASTA file imported from imtg website
    :return:
    Dictionary: containing information for each TCR gene
    """
    all_genes = []
    for i in raw_file:
        gene_seq = i.split('\n')
        attributes = gene_seq[0].split('|')
        seq = ''.join(gene_seq[1:]).upper()

        acc_no = attributes[0]
        temp = attributes[1].split('*')
        allele = temp[0]
        sub = temp[1]


        species_full = attributes[2].split('_')
        species = species_full[0]
        try:
            variety = species_full[1]
        except IndexError:
            variety = ''

        functional = attributes[3]
        exon_name = attributes[4]

        # Handling coordinates attribute, need to isolate start and end location values of gene
        coordinates = attributes[5]
        try:
            coordinates = get_imgt_coordinate(coordinates)
        except ValueError:
            coordinates = [0, 1]
        start = coordinates[0]
        end = coordinates[1]

        # remove 'nt' from length
        length = attributes[6]
        length_type = type(length)
        length = int(length_type().join(filter(length_type.isdigit, length)))

        start_codon = int(attributes[7])
        v1 = attributes[8]
        v2 = attributes[9]
        v3 = attributes[10]
        aa_length = attributes[11]

        # Convert into int
        seq_length = attributes[12]
        seq_length = seq_length.split('=')
        seq_length = int(seq_length[1])

        partial = attributes[13]
        rev_comp = attributes[14]

        # Make constant segment allele values unique

        if ('TRAC' in allele) or ('TRBC' in allele) or ('TRDC' in allele) or ('TRGC' in allele):
            comb = [allele, exon_name]
            allele = '_'.join(comb)

        chain = allele[2]
        segment = allele[3]
        protein_seq = translate(seq, start_codon)
        cdr3 = ''
        # combine into a dictionary
        attribute_list = {'acc_no': acc_no, 'allele': allele, 'sub': sub,'species': species,
                          'variety': variety,
                          'functional': functional, 'exon_name': exon_name, 'start': start,
                          'end': end, 'length': length, 'start_codon': start_codon,
                          'v1': v1, 'v2': v2, 'v3': v3, 'aa_length': aa_length,
                          'seq_length': seq_length, 'partial': partial, 'rev_comp': rev_comp,
                          'seq': seq,
                          'protein_seq': protein_seq, 'segment': segment, 'chain': chain,
                          'cdr3': cdr3}

        # Create a tuple which contains list of attributes for each gene and a string for gene sequence
        gene = attribute_list
        all_genes.append(gene)
    return all_genes
