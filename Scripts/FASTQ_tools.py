#!/usr/bin/env python3
"""
FASTQ_TOOL
"""

from typing import Union, Optional
import os


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


def filter_fastq(fastq_dict: dict,
                 gc_bounds: tuple,
                 length_bounds: tuple,
                 quality_threshold: float) -> dict:
    """
    applies filters, returns resulting dict with filtered entries

    Args:
    - fastq_dict - a read dictionary
    - gc_bounds - GC ratio threshold bounds
    - length_bounds - seqeunces length range
    - quality_threshold - a mean quality threshold (phred33)

    Returns:
    - filtered_fastq - a filtered dictionary with fastq
    """

    filtered_fastq = {}

    for seq_key, seq_values in fastq_dict.items():
        sequence, seq_info, seq_quality = (entry for entry in seq_values)
        gc_content = calc_gc(sequence)
        mean_quality = calc_mean_quality(seq_quality)
        if (gc_bounds[0] <= gc_content <= gc_bounds[1] and
                length_bounds[0] <= len(sequence) <= length_bounds[1] and
                mean_quality >= quality_threshold):
            filtered_fastq[seq_key] = [sequence, seq_info, seq_quality]
    return filtered_fastq


def read_fastq(input_path: str) -> dict:
    """
    reads a .fastq file and returns a dict

    Args:
    - input_path - a path to the .fastq file

    Returns:
    - fastq_dict - a dictionary where keys are the names of
    sequences and the values are tuples of strings
    """

    fastq_dict = {}

    with open(input_path, 'r', encoding='utf-8') as fastq_file:
        content = fastq_file.readlines()
        for line_num, line in enumerate(content):
            line = line.strip()
            if (line.startswith('@')) and (line_num % 4 in [0, 1]):
                key = line
                fastq_dict[key] = []
            else:
                fastq_dict[key].append(line)
    return fastq_dict


def write_fastq(filtered_fastq: dict,
                output_filename: str) -> None:
    """
    write filtered fastq file

    Args:
    - filtered_fastq - a filtered fastq dictionary
    - output_filename - path to the output file

    Returns:
    None
    """

    if not os.path.exists('fastq_filtrator_resuls'):
        os.makedirs('fastq_filtrator_resuls')
    out_path = os.path.join('fastq_filtrator_resuls', output_filename)
    with open(out_path, 'w', encoding='utf-8') as out_file:
        for entry in filtered_fastq.items():
            out_file.write(entry[0] + '\n')
            for seq_row in entry[1]:
                out_file.write(seq_row + '\n')


def run_fastq_tools(
        input_path: str,
        output_filename: Optional[str] = None,
        gc_bounds: Union[tuple, float] = (0, 100),
        length_bounds: Union[tuple, int] = (0, 2**32),
        quality_threshold: float = 0
        ) -> None:
    """
    a filter function to sort out the sequences that pass the setting.

    Args:
    - input_path - a file path to the .fastq file
    - output_filename - a file path to the output file. If none is given
    using the same name as input
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

    if output_filename is None:
        output_filename = os.path.split(input_path)[-1]

    fastq_dict = read_fastq(input_path)

    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)

    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)

    print(f'''Filtering parameters:
        - gc_bounds\t{gc_bounds}
        - length_bounds\t{length_bounds}
        - quality_threshold\t{quality_threshold}''')

    filtered_fastq = filter_fastq(fastq_dict,
                                  gc_bounds,
                                  length_bounds,
                                  quality_threshold)
    write_fastq(filtered_fastq, output_filename)


if __name__ == '__main__':
    pass
