"""
Read trimming for Accel-Amplicon Panel
(https://swiftbiosci.com/wp-content/uploads/2017/03/17-1397-Amplicon-Bioinf-Guidelines.pdf)

Requirements
Samples have to be Paired-End Illumina sequences.

Expects the global variable config of at least the following structure
...............................................................................
---

illuminaclip_file: illumina.fa

accel_panels:
  panel1:
    5p_primer_file: 5p_primers.fa
    3p_primer_file: 3p_primers.fa
  panel2:
    5p_primer_file: 5p_primers.fa
    3p_primer_file: 3p_primers.fa

...............................................................................

Expects a samples.tsv file
With the following columns, sample, panel, input_folder, sample_number and
lane. Sample must be the first column, the order of the other columns
aren't important. The panel column must contain a panel name that can be found
under accel_panels in the config.yaml file. input_folder isn't required, when
not specified it will be assumed that the fastq files are located under the
working directory.

Example of a sample.tsv, columns need to be tab separated
...............................................................................
sample     panel     fq1                fq2
sample1    panel1    sample1.R1.fastq   sample1.R2.fastq
sample2    panel2    sample2.R1.fastq   sample2.R2.fastq
...............................................................................

Output files will be:
1 - seqdata_trimmed/{sample}.R1.trimmomatic_cutadapt.fastq.gz
2 - seqdata_trimmed/{sample}.R2.trimmomatic_cutadapt.fastq.gz
3 - seqdata_trimmed/{sample}.trimmomatic_cutadapt.qc.txt
"""

__author__ = "Patrik Semds"
__license__ = "MIT"

import pandas as pd

configfile: "config.yaml"
samples = pd.read_table("samples.tsv", index_col=0)

file_endings = [".R1.trimmomatic_cutadapt.fastq.gz", ".R2.trimmomatic_cutadapt.fastq.gz", ".trimmomatic_cutadapt.qc.txt"]

def generate_file_output():
    return [os.path.join("seqdata_trimmed", str(row.Index) + ending)
        for row in samples.itertuples()
            for ending in file_endings]

rule all:
    input:
        generate_file_output()

include: "rules/accel_amplicon.smk"
