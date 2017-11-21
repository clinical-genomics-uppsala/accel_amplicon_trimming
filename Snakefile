__author__ = "Patrik Semds"
__license__ = "MIT"

import pandas as pd

configfile: "config.yaml"
samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

files = [("trimmed",".R1.trimmomatic_cutadapt.fastq.gz"),
             ("trimmed",".R2.trimmomatic_cutadapt.fastq.gz"),
                 ("logs/trimmed",".trimmomatic_cutadapt.qc.txt")]

def generate_file_output():
    return [f[0] + "/" + str(row.Index[0]) + "-" + str(row.Index[1]) + f[1]
        for row in units.itertuples()
            for f in files]

rule all:
    input:
        generate_file_output()

include: "rules/accel_amplicon.smk"
