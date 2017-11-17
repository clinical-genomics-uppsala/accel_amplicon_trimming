__author__ = "Patrik Semds"
__license__ = "MIT"

import pandas as pd

configfile: "config.yaml"
samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

file_endings = [".R1.trimmomatic_cutadapt.fastq.gz", ".R2.trimmomatic_cutadapt.fastq.gz", ".trimmomatic_cutadapt.qc.txt"]

def generate_file_output():
    return [os.path.join("trimmed", str(row.Index[0]) + "-" + str(row.Index[1]) + ending)
        for row in units.itertuples()
            for ending in file_endings]

rule all:
    input:
        generate_file_output()

include: "rules/accel_amplicon.smk"
