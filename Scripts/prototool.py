#!/usr/bin/env python3

"""
prototool
"""

from typing import List, Optional, Tuple, Union


AMINOACIDS_DICT = {
    'Ala': {'one_letter_code': 'A',
            'coding_RNAs': {'GCU', 'GCC', 'GCA', 'GCG'},
            'pKa_aminoacids': 0.0,
            'molecular_weights': 89},
    'Arg': {'one_letter_code': 'R',
            'coding_RNAs': {'CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'},
            'pKa_aminoacids': 12.4,
            'molecular_weights': 174},
    'Asn': {'one_letter_code': 'N',
            'coding_RNAs': {'AAU', 'AAC'},
            'pKa_aminoacids': 3.2,
            'molecular_weights': 132},
    'Asp': {'one_letter_code': 'D',
            'coding_RNAs': {'GAU', 'GAC'},
            'pKa_aminoacids': 3.9,
            'molecular_weights': 133},
    'Cys': {'one_letter_code': 'C',
            'coding_RNAs': {'UGU', 'UGC'},
            'pKa_aminoacids': 8.3,
            'molecular_weights': 121},
    'Glu': {'one_letter_code': 'E',
            'coding_RNAs': {'GAA', 'GAG'},
            'pKa_aminoacids': 4.3,
            'molecular_weights': 147},
    'Gln': {'one_letter_code': 'Q',
            'coding_RNAs': {'CAA', 'CAG'},
            'pKa_aminoacids': 3.2,
            'molecular_weights': 146},
    'Gly': {'one_letter_code': 'G',
            'coding_RNAs': {'GGU', 'GGC', 'GGA', 'GGG'},
            'pKa_aminoacids': 0.0,
            'molecular_weights': 75},
    'His': {'one_letter_code': 'H',
            'coding_RNAs': {'CAU', 'CAC'},
            'pKa_aminoacids': 6.0,
            'molecular_weights': 155},
    'Ile': {'one_letter_code': 'I',
            'coding_RNAs': {'AUU', 'AUC', 'AUA'},
            'pKa_aminoacids': 0.0,
            'molecular_weights': 131},
    'Leu': {'one_letter_code': 'L',
            'coding_RNAs': {'CUU', 'CUC', 'CUA', 'CUG'},
            'pKa_aminoacids': 0.0,
            'molecular_weights': 131},
    'Lys': {'one_letter_code': 'K',
            'coding_RNAs': {'AAA', 'AAG'},
            'pKa_aminoacids': 10.5,
            'molecular_weights': 146},
    'Met': {'one_letter_code': 'M',
            'coding_RNAs': {'AUG'},
            'pKa_aminoacids': 0.0,
            'molecular_weights': 149},
    'Phe': {'one_letter_code': 'F',
            'coding_RNAs': {'UUU', 'UUC'},
            'pKa_aminoacids': 0.0,
            'molecular_weights': 165},
    'Pro': {'one_letter_code': 'P',
            'coding_RNAs': {'CCU', 'CCC', 'CCA', 'CCG'},
            'pKa_aminoacids': 0.0,
            'molecular_weights': 115},
    'Ser': {'one_letter_code': 'S',
            'coding_RNAs': {'UCU', 'UCC', 'UCA', 'UCG'},
            'pKa_aminoacids': 2.2,
            'molecular_weights': 105},
    'Thr': {'one_letter_code': 'T',
            'coding_RNAs': {'ACU', 'ACC', 'ACA', 'ACG'},
            'pKa_aminoacids': 2.6,
            'molecular_weights': 119},
    'Tyr': {'one_letter_code': 'W',
            'coding_RNAs': {'UAU', 'UAC'},
            'pKa_aminoacids': 11.0,
            'molecular_weights': 181},
    'Trp': {'one_letter_code': 'Y',
            'coding_RNAs': {'UGG'},
            'pKa_aminoacids': 10.1,
            'molecular_weights': 204},
    'Val': {'one_letter_code': 'V',
            'coding_RNAs': {'GUU', 'GUC', 'GUA', 'GUG'},
            'pKa_aminoacids': 0.0,
            'molecular_weights': 117},
}


# A dictionary where keys are 1-letter and values are 3-letters codes
TO_3_DICT = {nested_dict['one_letter_code']: key for key,
             nested_dict in AMINOACIDS_DICT.items()}


def is_one_letter(seq: str) -> bool:
    """
    Defines whether the sequence is 1 coded.

    Args:
    - seq - sequence to check

    Returns:
    - bool
    """
    return all(aa.isalpha() and aa.isupper() for aa in seq)


def recode(seq: str) -> str:
    """
    Translate 1-letter to 3-letter encoding if 1-letter
    encoded sequence is given and vice versa.

    Args:
    - seq - sequence to recode

    Returns:
    - either a 3- or a 1-letter recoded sequence
    """

    if is_one_letter(seq):
        # Translate 1-letter to 3-letter coded sequence
        three_letter_sequence = ""
        for aa in seq:
            three_letter_code = TO_3_DICT.get(aa, aa)
            three_letter_sequence += three_letter_code
        return three_letter_sequence
    # Translate 3-letter to 1-letter coded sequence
    one_letter_sequence = ""
    for aa in range(0, len(seq), 3):
        amino_acid = seq[aa:aa+3]
        one_letter_sequence += AMINOACIDS_DICT[amino_acid]['one_letter_code']
    return one_letter_sequence


def prettify_alignment(aligned_seq_on: str, aligned_seq2: str) -> None:
    """
    Prettifies alignment output by printing out two
    sequences on top of each other

    Finds the start of aligned sequence in the longer of sequences.\\
    Prints the longer sequence as an upper one and aligned sequence
    is bellow separated via vertical lines

    Args:
    - aligned_seq_on, aligned_seq2 - sequences
    from the local_alignment()

    Returns:
    None \\
    Prints out the prettified view in stdout
    """

    print(aligned_seq_on)
    print('|' * len(aligned_seq2))
    print(aligned_seq2)


def local_alignment(seq_on: str,
                    seq2: str,
                    match=2,
                    mismatch=-1,
                    gap=-1,
                    prettify: bool = True) -> dict:
    """
    Perform a local alignment of 2 given sequences

    Args:
    - seq_on - the sequence to align onto
    - seq2 - sequences to align
    - match, mismatch, gap - alignment scoring and penalty values
    defaulted to 2, -1, -1
    - prettify - if True (default) prints out the prettified version
    of sequences aligned on top of each other

    Returns:
    - a a dictionary with alignment resluts
    """

    len_seq_on, len_seq2 = len(seq_on), len(seq2)

    # Initialize the score matrix and traceback matrix
    score_matrix = [[0] * (len_seq2 + 1) for _ in range(len_seq_on + 1)]
    traceback_matrix = [[None] * (len_seq2 + 1) for _ in range(len_seq_on + 1)]

    alignment_score = 0  # To keep track of the maximum score in the matrix
    max_i, max_j = 0, 0  # To store the position of the maximum score

    # Fill in the score matrix
    for i in range(1, len_seq_on + 1):
        for j in range(1, len_seq2 + 1):
            if seq_on[i - 1] == seq2[j - 1]:
                match_score = score_matrix[i - 1][j - 1] + match
            else:
                match_score = score_matrix[i - 1][j - 1] + mismatch

            delete_score = score_matrix[i - 1][j] + gap
            insert_score = score_matrix[i][j - 1] + gap

            # Calculate the maximum score for the current cell
            score = max(0, match_score, delete_score, insert_score)

            # Update the score matrix and traceback matrix
            score_matrix[i][j] = score

            if score > alignment_score:
                alignment_score = score
                max_i, max_j = i, j

            if score == match_score:
                traceback_matrix[i][j] = "match"
            elif score == delete_score:
                traceback_matrix[i][j] = "delete"
            elif score == insert_score:
                traceback_matrix[i][j] = "insert"
            else:
                traceback_matrix[i][j] = "none"

    # Traceback to find the aligned sequences
    aligned_seq_on = []
    aligned_seq2 = []

    counter_identity: int = 0
    counter_gaps: int = 0

    i, j = max_i, max_j

    while i >= 0 and j >= 0:
        if i == 0 or j == 0 or traceback_matrix[i][j] == "none":
            break
        if traceback_matrix[i][j] == "match":
            aligned_seq_on.append(seq_on[i - 1])
            aligned_seq2.append(seq2[j - 1])
            counter_identity += 1
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == "delete":
            aligned_seq_on.append(seq_on[i - 1])
            aligned_seq2.append("-")
            counter_gaps += 1
            i -= 1
        elif traceback_matrix[i][j] == "insert":
            aligned_seq_on.append("-")
            aligned_seq2.append(seq2[j - 1])
            counter_gaps += 1
            j -= 1

    # Reverse the aligned sequences
    aligned_seq_on = "".join(aligned_seq_on[::-1])
    aligned_seq2 = "".join(aligned_seq2[::-1])

    alignment_length = (len(aligned_seq_on)
                        if len(aligned_seq_on) < len(aligned_seq2)
                        else len(aligned_seq2))

    try:
        identity = round(counter_identity/alignment_length, 4)
    except ZeroDivisionError:
        print('It appears the alignment length of your sequences is 0.')
        identity = 'NaN'
        prettify = False

    # Prettify an alignment output
    if prettify is True:
        prettify_alignment(aligned_seq_on, aligned_seq2)
    else:
        pass

    return {'aligned_seq_on': aligned_seq_on,
            'aligned_seq': aligned_seq2,
            'length': alignment_length,
            'score': alignment_score,
            'identity': identity,
            'gaps': counter_gaps}


def calc_protein_molecular_weight(seq: str) -> float:
    """
    Args:
    - seq - a sequence to calculate molecular weight of

    Returns:
    protein_weight - molecular weight of a given protein sequence
    """

    protein_weight = 0
    aminoacids = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    for i, aminoacid in enumerate(aminoacids):
        if aminoacid in AMINOACIDS_DICT:
            aminoacid_weight = (AMINOACIDS_DICT[aminoacid]
                                ['molecular_weights'])
            protein_weight += aminoacid_weight
    return protein_weight


def from_proteins_seqs_to_rna(seq: str) -> set:
    """
    List all possible RNAs for a given protein sequence
    Args:
    - seq - a protein sequence (needs 3-letter coded)

    Returns:
    - rna_sequences - a set with all possible RNA sequences for a given
    protein sequence
    """

    rna_sequences = ['']

    for amino_acid in [seq[i:i + 3] for i in range(0, len(seq), 3)]:
        possible_codons = AMINOACIDS_DICT[amino_acid]['coding_RNAs']
        rna_sequences = [seq + codon for seq in rna_sequences
                         for codon in possible_codons]
    return set(rna_sequences)


def calculate_pI(seq: str) -> float:
    """
    Calculates an isoelectric point (pI) for a given sequence

    Args:
    - seq - a sequence of a protein to calculate pI for

    Returns:
    - pI - a pH when the net charge of a given sequence equals 0
    """
    pKa_NH2 = 8.0
    pKa_COOH = 3.1
    step = 0.01
    pH = 0
    charge = 0.0

    # Iterate through a range of pH values to find the isoelectric point
    while pH <= 14:
        # Calculate the charge contributed by ionizable groups at this pH
        for amino_acid in seq:
            if amino_acid == 'D' or amino_acid == 'E':
                charge += 1 / (1 + 10 ** (pKa_side_chains[amino_acid] - pH))
            elif amino_acid == 'H':
                charge += 1 / (1 + 10 ** (pH - pKa_side_chains[amino_acid]))
            elif amino_acid == 'K' or amino_acid == 'R':
                charge += 1 / (1 + 10 ** (pH - pKa_side_chains[amino_acid]))

        # Calculate the charge contributed by the N-terminus and C-terminus
        charge_NH2 = 1 / (1 + 10 ** (pKa_NH2 - pH))
        charge_COOH = -1 / (1 + 10 ** (pH - pKa_COOH))

        # Calculate the net charge at this pH
        net_charge = charge + charge_NH2 + charge_COOH

        # Check if the net charge is close to zero (within a tolerance)
        if abs(net_charge) < 0.01:
            return pH

        # Increment pH and reset the charge
        pH += step
        charge = 0

    # If no pI is found, return None
    return None


def is_protein(seq: str) -> bool:
    """
    Checks if all aminoacids are present in the dictionary

    Args:
    - seq - a sequence to check

    Returns:
    - bool
    """
    if not is_one_letter(seq):
        return all(aa in AMINOACIDS_DICT for aa
                   in [seq[i:i + 3] for i in range(0, len(seq), 3)])
    else:
        return all(aa in TO_3_DICT for aa in seq)


def check_input(*args: Union[List[str], str], method: str) -> \
                                    Tuple[List[str], Optional[str]]:
    """
    Function to check the validity of the input.

    Args:
    - *args - are supposed to be all sequences to process
    - method - the method to process with method

    Returns:
    - seqs_list - list of sequences
    - seq_on (optional) - in case of local_alignment method
    """

    if len(args) == 0:
        # Handle the case where there are no arguments
        raise ValueError('No input defined.')
    else:
        # Form a list with sequences from the input
        seqs_list = list(args)
        valid_seqs_list = []
        seq_on = None
        if method in ['from_proteins_seqs_to_rna',
                      'calculate_pI',
                      'calc_protein_molecular_weight']:
            # Needs protein 3-letter coded
            for seq in seqs_list:
                if not is_protein(seq):
                    print(seq, ' is not a valid protein.')
                    continue
                if is_one_letter(seq):
                    print(f'Warning! Function {method}() needs '
                          '3-letter encoded sequences. Your sequence '
                          'will be mutated to a 3-letter encoding.')
                    valid_seqs_list.append(recode(seq))
                    print(seq, ' sequence has been mutated into: ',
                          valid_seqs_list[-1])
                else:
                    valid_seqs_list.append(seq)
            return valid_seqs_list, seq_on

        elif method == 'local_alignment':
            # Needs at least 2 proteins 1-letter coded
            for seq in seqs_list:
                if not is_protein(seq):
                    print(seq, ' is not a valid protein.')
                    continue
                if not is_one_letter(seq):
                    print('Warning! Function local_alignment() needs '
                          '1-letter encoded sequences. Your sequence '
                          'will be mutated to a 1-letter encoding.')
                    print(seq, ' sequence has been mutated into: ',
                          valid_seqs_list[-1])
                    valid_seqs_list.append(recode(seq))
                else:
                    valid_seqs_list.append(seq)
            if len(valid_seqs_list) < 2:
                # Handle the case where there s only 1 sequence for
                # a local_alignment()
                raise IndexError('Need at least two sequences to align.')
            seq_on = valid_seqs_list.pop(0)
            return valid_seqs_list, seq_on

        elif method == 'recode':
            # Needs protein
            for seq in seqs_list:
                if not is_protein(seq):
                    print(seq, ' is not a valid protein.')
                    continue
                valid_seqs_list.append(seq)
            return valid_seqs_list, seq_on

        else:
            raise ValueError(method, ' is not a valid method.')


def prototools(*args: Tuple[Union[List[str], str]],
               method: Optional[str] = None) -> dict:
    """
    This function provides the access to the following methods:

    1. Translate 1 letter to 3 letter encoding and vice versa - the last
    argument: 'recode'
        - needs at least 1 sequence 1- or 3- letter encoded.
        - returns a recoded sequence

    2. Local Alignment of two sequences - the last argument: 'local_alignment'
       - needs 2 protein sequences 1-letter encoded.
       - performs an alignment using Smith-Waterman algorithm

    3. Find possible RNA sequences for defined protein sequence - the
    last argument: 'from_proteins_seqs_to_rna'
        - needs an input containing at least 1 3-letter coded aminoacid
        sequence
        - returns a set of RNA codones, which encode this protein

    4. Determinate isoelectric point - the last argument:
    'calculate_pI'
        - needs an input containing at least 1 3-letter coded aminoacid
        sequence
        - returns a value is an isoelectric point of a given protein

    5. Calculate protein molecular weight - the last argument:
    'calc_protein_molecular_weight'
        - needs an input containing at least 1 3-letter coded aminoacid
        sequence
        - returns a calculated molecular weight for a given sequence

    Args:
    - *args - are supposed to be all sequences to process
    - method is a kwarg - the method to process with.

    Returns:
    function_result - a dictionary with the result of a chosen function
    """

    seqs_list, seq_on = check_input(*args, method=method)
    print(f'Your sequences are: {seqs_list}',
          f'The method is: {method}', sep='\n')

    match method:
        case 'recode':
            recode_dict: dict = {}
            for seq in seqs_list:
                recode_dict[seq] = recode(seq)
            return recode_dict

        case 'local_alignment':
            alignments_dict: dict = {}
            print(f'Your sequence align on is: {seq_on}')
            for seq in seqs_list:
                alignments_dict[seq] = local_alignment(seq_on=seq_on,
                                                       seq2=seq,
                                                       prettify=True)
            return alignments_dict

        case 'from_proteins_seqs_to_rna':
            possible_rnas_dict: dict = {}
            for seq in seqs_list:
                possible_rnas_dict[seq] = from_proteins_seqs_to_rna(seq)
            return possible_rnas_dict

        case 'calc_protein_molecular_weight':
            molecular_weights_dict: dict = {}
            for seq in seqs_list:
                molecular_weights_dict[seq] = \
                    calc_protein_molecular_weight(seq)
            return molecular_weights_dict

        case 'calculate_pI':
            iet_dict: dict = {}
            for seq in seqs_list:
                iet_dict[seq] = calculate_pI(seq)
            return iet_dict


if __name__ == '__main__':
    pass




test1 = prototools('AlaTyr', 'A', 'ljadlksg', method='calculate_pI')
print(test1)
