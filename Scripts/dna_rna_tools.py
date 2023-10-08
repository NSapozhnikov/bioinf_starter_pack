#!/usr/bin/env python
"""
dna_rna_tools
"""
from typing import List


def complement(sequence: str, nucleic_acid: str) -> str:
    """
    Returns a complementary sequence to the given

    Args:
    - sequence
    - nucleic_acid - a type of nucleic acid ('DNA'/'RNA')

    Returns:
    - complementary_sequence
    """

    complementary_sequence = (sequence
                              .translate(str.maketrans('ATGC', 'TACG')))
    complementary_sequence = (complementary_sequence
                              .translate(str.maketrans('atgc', 'tacg')))

    if nucleic_acid == 'RNA':
        complementary_sequence = (sequence
                                  .translate(str.maketrans('AUGC', 'UACG')))
        complementary_sequence = (complementary_sequence
                                  .translate(str.maketrans('augc', 'uacg')))
    return complementary_sequence


def transcribe(sequence: str,
               nucleic_acid: str,
               trans_rna: bool = True) -> str:
    """
    Returns a transcribed sequence to the given. Works with both DNA
    and RNA

    Args:
    - sequence - a sequence to transcribe
    - nucleic_acid - type of a sequence ('DNA'/'RNA')
    - trans_RNA - defaulted to True - wether to transcribe RNA
    sequences

    Returns:
    - transcribed_sequence
    """

    if nucleic_acid == 'DNA':
        transcribed_sequence = sequence.replace('T', 'U')
        transcribed_sequence = transcribed_sequence.replace('t', 'u')
    else:
        if trans_rna:
            transcribed_sequence = sequence.replace('U', 'T')
            transcribed_sequence = transcribed_sequence.replace('u', 't')
    return transcribed_sequence


def reverse(seq: str) -> str:
    """
    This function reverses the sequence.

    Args:
    - seq - is a DNA or RNA sequence

    Returns:
    - a string with a reversed sequence
    """
    return seq[::-1]


def check_input(*seqs: str) -> List[str]:

    """
    Function to check the validity of the input.

    Args:
    *seqs - is supposed to be all sequences to process and the method
    to process with.
    The method is supposed to be the last argument.

    Returns:
    - seqs_list - a list of all sequences
    - method - method to use
    """

    if len(seqs) < 1:
        # Handle the case where there are no arguments
        raise ValueError("No input defined.")
    # Check the last element of the input is a valid method
    method = seqs[-1]
    if method not in ['transcribe',
                      'reverse',
                      'complement',
                      'reverse_complement']:

        raise ValueError(method, " is not a valid method.")
    return list(seqs[:-1]), method


def check_seq(seq: str) -> str:
    """
    This function checks whether the sequence is a DNA or an RNA.

    Args:
    - seq - is a sequence
    It can contain characters typical for a DNA or an RNA. If any other symbols
    are present the function raise a ValueError.

    Returns:
    - a string of type 'DNA' or 'RNA'
    """

    # Define the if sequence is a DNA
    if all(char in 'ATGCatgc' for char in seq):
        return 'DNA'
    # Define the if sequence is an RNA
    elif all(char in 'AUGCaugc' for char in seq):
        return 'RNA'
    else:
        raise ValueError(seq,
                         " sequence is not a valid nucleic acid sequence.")


def run_dna_rna_tools(*seqs: List[str]) -> List[str]:
    """
    This function processes sequences with one of the available methods:
    1. 'transcribe'
    2. 'reverse'
    3. 'complement'
    4. 'reverse_complement'

    Args:
    - seqs - DNA or RNA sequences
    - method - a valid method to use

    Returns:
    - processed sequence or a list of processed sequences
    """

    seqs_list, method = check_input(*seqs)
    output_list = []
    # Process sequences according to the method
    for seq in seqs_list:
        nucleic_acid = check_seq(seq)
        match method:
            case 'transcribe':
                result_seq = transcribe(seq,
                                        nucleic_acid,
                                        trans_rna=True)
            case "reverse":
                result_seq = reverse(seq)
            case "complement":
                result_seq = complement(seq, nucleic_acid)
            case "reverse_complement":
                result_seq = reverse(seq)
                result_seq = complement(result_seq, nucleic_acid)

        output_list.append(result_seq)
    return output_list


if __name__ == '__main__':
    pass
