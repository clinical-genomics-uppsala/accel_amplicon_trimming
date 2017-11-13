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
