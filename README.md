# Bioinf_starter_pack

is a toolkit to process amino- and nucleic acids.

### Overview
1. `run_FASTQ_tools` includes a function to filter by:
  - GC-content
  - sequence length
  - sequence average quality

2. `prototools` includes 7 methods to treatment of polyaminoacid sequences.
    `run_prototools` can be used for the following purposes:
  - recoding 1-letter coded polyaminoacid seqeunces into 3-letter coded and vice versa;
  - polyaminoacid sequences aligment with Smith-Waterman algorithm [^1];
  - finding possinle RNA sequences for given polyaminoacid sequences;
  - determining polyaminoacid isoelectric point;
  - calculating polyaminoacid molecular weight;
  - finding possinle DNA sequences for given polyaminoacid sequences; 
  - determining GC-content of a corresponding DNA sequence to a given polyaminoacid sequence

3. `run_dna_rna_tools` can be used to perform simple manipulations on nucleic acid sequences:
  - 'transcribe'
  - 'reverse'
  - 'complement'
  - 'reverse_complement'


### Usage
This tool can be used both standalone and as module (the latter is recommended).
1. General usecase:
   - `python import bioin_starter_pack`
   - `bioinf_starter_pack.run_FASTQ_tools(example_fastq)`
  

  ![image](https://github.com/NSapozhnikov/bioinf_starter_pack/assets/81642791/8a780e2a-6c55-4bc5-9ce1-1a4f0ba81d54)
2. You can also use any of the secondary functions e.g.: `from Scripts.prototools import local_alignment`. *Check the wiki for all methods description.*
- the module is logically devided into 3 parts:
  - `dna_rna_tools` - to work with simple DNA and RNA sequences
  - `prototools` - to work with polyaminoacid sequences
  - `FASTQ_tools` - to work with FASTQ format

![image](https://github.com/NSapozhnikov/bioinf_starter_pack/assets/81642791/d9040d8f-4adb-45b1-a9fa-bbeea89d17d2)


### Troubleshooting
Note that all checks for validity are performed in the main functions(run_dna_rna_tools(), run_prototools(), run_FASTQ_tools()). If you are using any secondary function please read its docstring and check your input for validity or the function may work incorrectly.
e.g.:

![image](https://github.com/NSapozhnikov/bioinf_starter_pack/assets/81642791/929cb304-b21d-4869-8e75-60099709a14e)

### Contributions and contacts

Feel free to report any bugs and problems encountered. Any bug reported is appreciated.
Email: nikita.sapozhnikov1@gmail.com 

### References

[^1]: T.F. Smith, M.S. Waterman, (1981). [Identification of common molecular subsequences](https://doi.org/10.1016/0022-2836(81)90087-5). Journal of Molecular Biology.
