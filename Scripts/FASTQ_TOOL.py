#!/usr/bin/env python3
"""
FASTQ_TOOL
"""

from typing import List, Union, Tuple, Optional


def calc_gc(seq: str) -> float:
    """
    calculates GC ratio of a given sequence

    Args:
    seq - a sequence to use

    Returns:
    gc_ratio - a GC-content of a given sequence
    """

    gc_ratio = round(((seq.count('G') + seq.count('C'))/len(seq)) * 100, 2)
    return gc_ratio


def calc_mean_quality(seq_quality: str) -> float:
    """
    calculates a mean quality of a sequence using a phred33

    Args:
    seq_quality - sequence quality

    Returns:
    mean_seq_quality - mean of all scores of a sequence
    """

    mean_quality = sum((ord(score) - 33)
                       for score in seq_quality) / len(seq_quality)
    return round(mean_quality, 2)


def fastq_filter(
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

    filtered_fastq = {}

    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)

    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)

    print(f'''Filtering parameters:
        - gc_bounds\t{gc_bounds}
        - length_bounds\t{length_bounds}
        - quality_threshold\t{quality_threshold}''')

    for seq_name in seqs:
        sequence, seq_quality = seqs[seq_name][0], seqs[seq_name][1]
        gc_content = calc_gc(sequence)
        mean_quality = calc_mean_quality(seq_quality)
        if (gc_bounds[0] <= gc_content <= gc_bounds[1] and
                length_bounds[0] <= len(sequence) <= length_bounds[1] and
                mean_quality >= quality_threshold):
            filtered_fastq[seq_name] = {'sequence': sequence,
                                        'sequence_length': len(sequence),
                                        'gc-content': gc_content,
                                        'mean_quality': mean_quality}

    return filtered_fastq


if __name__ == '__main__':
    pass
