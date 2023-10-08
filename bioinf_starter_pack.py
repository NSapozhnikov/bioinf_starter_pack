#!/usr/bin/env python3
"""
BioInf_Starter_Pack

is a python toolkit to interact with nucleic acid and protein sequences.
"""

from typing import List, Union, Tuple, Optional
from Scripts.dna_rna_tools import run_dna_rna_tools
from Scripts.FASTQ_tools import run_FASTQ_tools
from Scripts.prototools import run_prototools


def proto_tools(*args: Tuple[Union[List[str], str]],
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
    - method is a kwarg - the method to process with.

    Returns:
    function_result - a dictionary with the result of a chosen function
    """

    return run_prototools(*args, method=method)


def filter_fastq(
        seqs: dict,
        gc_bounds: Union[tuple, float] = (0, 100),
        length_bounds: Union[tuple, float] = (0, 2**32),
        quality_threshold: int = 0
        ) -> dict:
    """
    a filter function to sort out the sequences that pass the setting.

    Args:
    - *seqs - an unlimited amount of sequences in a dictionary where key is
    the name of a sequence and the value is a tuple of strings
    (sequence, quality)
    - gc_bounds - defaulted to (0, 100) - GC ratio threshold. Can be either
    a tuple or a float, if latter it will be the upper boundary
    - length_bounds - defaulted to (0, 2**32) - seqeunces length range. Can be
    either a tuple or a float, if latter it will be the upper boundary
    - quality_threshold - defaulted to 0 - a mean quality threshold (phred33).
    All sequences with a mean quality lower than set by this parameter will
    be filtered out

    Returns:
    filtered_fastq - result of a filtering a dictionary with sequence names as
    keys and sequences as values

    """

    return run_FASTQ_tools(seqs=seqs,
                           gc_bounds=gc_bounds,
                           length_bounds=length_bounds,
                           quality_threshold=quality_threshold)


def nucleic_tools(*seqs: List[str]) -> List[str]:
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

    return run_dna_rna_tools(*seqs)
