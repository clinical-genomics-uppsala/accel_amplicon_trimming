__author__ = "Patrik Semds"
__license__ = "MIT"

import pandas as pd

configfile: "config.yaml"
samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

rule all:
    input:
        expand(["trimmed/{unit.sample}.{unit.unit}.R{read}.trimmomatic_cutadapt.fastq.gz",
                "qc/trimmed/{unit.sample}.{unit.unit}.trimmomatic_cutadapt.qc.txt"],
               unit=units.reset_index().itertuples(),
               read=[1, 2])

include: "rules/accel_amplicon.smk"
