#!/usr/bin/env python3

"""
prototool
"""

from typing import List, Optional, Tuple, Union


AMINOACIDS_DICT = {
    'Ala': {'one_letter_code': 'A',
            'coding_RNAs': {'GCU', 'GCC', 'GCA', 'GCG'},
            'pKa': 6.01,
            'molecular_weights': 89},
    'Arg': {'one_letter_code': 'R',
            'coding_RNAs': {'CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'},
            'pKa': 10.76,
            'molecular_weights': 174},
    'Asn': {'one_letter_code': 'N',
            'coding_RNAs': {'AAU', 'AAC'},
            'pKa': 5.41,
            'molecular_weights': 132},
    'Asp': {'one_letter_code': 'D',
            'coding_RNAs': {'GAU', 'GAC'},
            'pKa': 2.85,
            'molecular_weights': 133},
    'Cys': {'one_letter_code': 'C',
            'coding_RNAs': {'UGU', 'UGC'},
            'pKa': 5.05,
            'molecular_weights': 121},
    'Glu': {'one_letter_code': 'E',
            'coding_RNAs': {'GAA', 'GAG'},
            'pKa': 3.15,
            'molecular_weights': 147},
    'Gln': {'one_letter_code': 'Q',
            'coding_RNAs': {'CAA', 'CAG'},
            'pKa': 5.65,
            'molecular_weights': 146},
    'Gly': {'one_letter_code': 'G',
            'coding_RNAs': {'GGU', 'GGC', 'GGA', 'GGG'},
            'pKa': 6.06,
            'molecular_weights': 75},
    'His': {'one_letter_code': 'H',
            'coding_RNAs': {'CAU', 'CAC'},
            'pKa': 7.60,
            'molecular_weights': 155},
    'Ile': {'one_letter_code': 'I',
            'coding_RNAs': {'AUU', 'AUC', 'AUA'},
            'pKa': 6.05,
            'molecular_weights': 131},
    'Leu': {'one_letter_code': 'L',
            'coding_RNAs': {'CUU', 'CUC', 'CUA', 'CUG'},
            'pKa': 6.01,
            'molecular_weights': 131},
    'Lys': {'one_letter_code': 'K',
            'coding_RNAs': {'AAA', 'AAG'},
            'pKa': 9.60,
            'molecular_weights': 146},
    'Met': {'one_letter_code': 'M',
            'coding_RNAs': {'AUG'},
            'pKa': 5.74,
            'molecular_weights': 149},
    'Phe': {'one_letter_code': 'F',
            'coding_RNAs': {'UUU', 'UUC'},
            'pKa': 5.49,
            'molecular_weights': 165},
    'Pro': {'one_letter_code': 'P',
            'coding_RNAs': {'CCU', 'CCC', 'CCA', 'CCG'},
            'pKa': 6.30,
            'molecular_weights': 115},
    'Ser': {'one_letter_code': 'S',
            'coding_RNAs': {'UCU', 'UCC', 'UCA', 'UCG'},
            'pKa': 5.68,
            'molecular_weights': 105},
    'Thr': {'one_letter_code': 'T',
            'coding_RNAs': {'ACU', 'ACC', 'ACA', 'ACG'},
            'pKa': 5.60,
            'molecular_weights': 119},
    'Tyr': {'one_letter_code': 'W',
            'coding_RNAs': {'UAU', 'UAC'},
            'pKa': 5.64,
            'molecular_weights': 181},
    'Trp': {'one_letter_code': 'Y',
            'coding_RNAs': {'UGG'},
            'pKa': 5.89,
            'molecular_weights': 204},
    'Val': {'one_letter_code': 'V',
            'coding_RNAs': {'GUU', 'GUC', 'GUA', 'GUG'},
            'pKa': 6.0,
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


def calc_mw(seq: str) -> float:
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


def prot_to_rna(seq: str) -> set:
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


def calculate_pI(seq: str,
                 pH_range: tuple = (0, 14),
                 pH: float = 6.5,
                 pH_step: float = 0.01) -> float:
    """
    Calculates an isoelectric point (pI) for a given sequence

    Args:
    - seq - a sequence of a protein to calculate pI for
    - pH_range - defaulted to (0, 14), lower and upper boundaries for pH
    to look in
    - pH - supossed starting pH value
    - pH_step - an incremental step for the search

    Returns:
    - pI - a pH when the net charge of a given sequence equals 0
    """

    pH_lower, pH_upper = pH_range[0], pH_range[1]

    asp_count: int = 0
    glu_count: int = 0
    cys_count: int = 0
    tyr_count: int = 0
    his_count: int = 0
    lys_count: int = 0
    arg_count: int = 0

    for aminoacid in [seq[i:i + 3] for i in range(0, len(seq), 3)]:
        match aminoacid:
            case 'Asp':
                asp_count += 1
            case 'Glu':
                glu_count += 1
            case 'Cys':
                cys_count += 1
            case 'Tyr':
                tyr_count += 1
            case 'His':
                his_count += 1
            case 'Lys':
                lys_count += 1
            case 'Arg':
                arg_count += 1

    while not ((pH - pH_lower < pH_step) &
               (pH_upper - pH < pH_step)):
        # Negatives
        c_term = -1 / (1 + 10 ** (3.65 - pH))
        asp_q = -asp_count / (1 + 10 ** (AMINOACIDS_DICT['Asp']['pKa'] - pH))
        glu_q = -glu_count / (1 + 10 ** (AMINOACIDS_DICT['Glu']['pKa'] - pH))
        cys_q = -cys_count / (1 + 10 ** (AMINOACIDS_DICT['Cys']['pKa'] - pH))
        tyr_q = -tyr_count / (1 + 10 ** (AMINOACIDS_DICT['Tyr']['pKa'] - pH))
        # Positives
        n_term = 1 / (1 + 10 ** (pH - 8.2))
        his_q = his_count / (1 + 10 ** (pH - AMINOACIDS_DICT['His']['pKa']))
        lys_q = lys_count / (1 + 10 ** (pH - AMINOACIDS_DICT['Lys']['pKa']))
        arg_q = arg_count / (1 + 10 ** (pH - AMINOACIDS_DICT['Arg']['pKa']))

        # Net charge in given pH
        net_charge: float = (c_term +
                             asp_q +
                             glu_q +
                             cys_q +
                             tyr_q +
                             his_q +
                             n_term +
                             lys_q +
                             arg_q)

        # We are out of range, thus the new pH value must be smaller
        if net_charge < 0:
            temp = pH
            pH = pH - ((pH - pH_lower) / 2)
            pH_upper = temp
        # We used to small pH value, so we have to increase it
        else:
            temp = pH
            pH = pH + ((pH_upper - pH) / 2)
            pH_lower = temp
        if abs(net_charge) < 0.001:
            break
    return pH


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
        if method in ['prot_to_rna',
                      'calculate_pI',
                      'calc_mw']:
            # Needs protein 3-letter coded
            for seq in seqs_list:
                if not is_protein(seq):
                    print(seq, ' is not a valid protein.')
                    continue
                if is_one_letter(seq):
                    valid_seqs_list.append(recode(seq))
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


def run_prototools(*args: Tuple[Union[List[str], str]],
                   method: Optional[str] = None,
                   prettify: bool = True) -> dict:
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
    last argument: 'prot_to_rna'
        - needs an input containing at least 1 3-letter coded aminoacid
        sequence
        - returns a set of RNA codones, which encode this protein

    4. Determinate isoelectric point - the last argument:
    'calculate_pI'
        - needs an input containing at least 1 3-letter coded aminoacid
        sequence
        - returns a value is an isoelectric point of a given protein

    5. Calculate protein molecular weight - the last argument:
    'calc_mw'
        - needs an input containing at least 1 3-letter coded aminoacid
        sequence
        - returns a calculated molecular weight for a given sequence

    Args:
    - *args - are supposed to be all sequences to process
    - method - the method to process with
    - prettify - defaulted to True - wether to prettify the
    alignment

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
                                                       prettify=prettify)
            return alignments_dict

        case 'prot_to_rna':
            possible_rnas_dict: dict = {}
            for seq in seqs_list:
                possible_rnas_dict[seq] = prot_to_rna(seq)
            return possible_rnas_dict

        case 'calc_mw':
            molecular_weights_dict: dict = {}
            for seq in seqs_list:
                molecular_weights_dict[seq] = \
                    calc_mw(seq)
            return molecular_weights_dict

        case 'calculate_pI':
            iet_dict: dict = {}
            for seq in seqs_list:
                iet_dict[seq] = calculate_pI(seq)
            return iet_dict


if __name__ == '__main__':
    pass
