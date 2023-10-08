# Bioinf_starter_pack

is a toolkit to process amino- and nucleic acids.

### FASTQ_tools

  A filter function to sort out the sequences that pass the setting.

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

### prototools

This function provides the access to the following methods:

  1. Translate 1 letter to 3 letter encoding and vice versa - the last
  argument: 'recode'
      - needs at least 1 sequence 1- or 3- letter encoded.
      - returns a recoded sequence

  2. Local Alignment of two sequences - the last argument: 'local_alignment'
     - needs 2 protein sequences 1-letter encoded.
     - performs an alignment using Smith-Waterman algorithm [^1]

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

### dna_rna_tools
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

### Overview
`prototool.py` includes 7 methods to treatment of polyaminoacid sequences.
`prototool.py` can be used for the following purposes:
- recoding 1-letter coded polyaminoacid seqeunces into 3-letter coded and vice versa;
- polyaminoacid sequences aligment with Smith-Waterman algorithm [^1];
- finding possinle RNA sequences for given polyaminoacid sequences;
- determining polyaminoacid isoelectric point;
- calculating polyaminoacid molecular weight;
- finding possinle DNA sequences for given polyaminoacid sequences; 
- determining GC-content of a corresponding DNA sequence to a given polyaminoacid sequence

***

### Usage
This tool can be used both standalone and as module (the latter is recommended).
- to use the tool as a module: `python import bioin_starter_pack`
- you can also use any of the secondary functions e.g.: `from Scripts.prototools import local_alignment`. *Check the wiki for all methods description.*
- the module is logically devided into 3 parts:
  - `dna_rna_tools` - to work with simple DNA and RNA sequences
  - `prototools` - to work with polyaminoacid sequences
  - `FASTQ_tools` - to work with FASTQ format


### Examples



### Troubleshooting
Please read function docstrings to pass valid inputs. Any bug reported is also appreciated

### Contributions and contacts

Feel free to report any bugs and problems encountered.
Email: nikita.sapozhnikov1@gmail.com 

### References

[^1]: T.F. Smith, M.S. Waterman, (1981). [Identification of common molecular subsequences](https://doi.org/10.1016/0022-2836(81)90087-5). Journal of Molecular Biology.
