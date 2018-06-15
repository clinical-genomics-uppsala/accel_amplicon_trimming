# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Patrik Smeds"
__license__ = "MIT"

###############################################################################
# Step 1:
# Remove any Illummina adaptor sequence using Trimmomatic
# (http://www.usadellab.org/cms/?page=trimmomatic).
###############################################################################

def _get_fastq(wildcards,units,read_pair='fq1'):
    return units.loc[(wildcards.sample, wildcards.unit), [read_pair]].dropna()[0]

rule trimmomatic:
    input:
        r1 = lambda wildcards: _get_fastq(wildcards, units, 'fq1'),
        r2 = lambda wildcards: _get_fastq(wildcards, units, 'fq2')
    output:
        temp("logs/trimmed/{sample}.{unit}.trimmomatic.qc.txt"),
        r1 = temp("trimmed/{sample}.{unit}.R1.trimmomatic.fastq"),
        r2 = temp("trimmed/{sample}.{unit}.R2.trimmomatic.fastq"),
        r1_unpaired=temp("trimmed/{sample}.{unit}.R1.trimmomatic.up.fastq"),
        r2_unpaired=temp("trimmed/{sample}.{unit}.R2.trimmomatic.up.fastq")
    threads: 12
    log:
        "logs/trimmed/{sample}.{unit}.trimmomatic.qc.txt"
    params:
        extra="-threads 12",
        trimmer=["ILLUMINACLIP:" + config["illuminaclip_file"] + ":2:30:10", "MINLEN:30"]
    wrapper:
        "0.17.4/bio/trimmomatic/pe"


###############################################################################
# Step 2
# Anchored 5’ and 3’ trimming of primer sequences, using cutadapt with
# design specific primer files.
# (http://cutadapt.readthedocs.io/en/stable/guide.html)
################################################################################

################################################################################
# Step 1 and 2
# anchored 5’ trimming of primer sequences with 5 primer design file.
###############################################################################

rule cutadapt_step1:
    input:
        "trimmed/{sample}.{unit}.R1.trimmomatic.fastq",
         "trimmed/{sample}.{unit}.R2.trimmomatic.fastq"
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -g file:" + config["accel_panels"][samples["panel"][wildcards.sample]]["5p_primer_file"]
    output:
        fastq1=temp("trimmed/{sample}.{unit}.tmpR1.fastq"),
        fastq2=temp("trimmed/{sample}.{unit}.tmpR2.fastq"),
        qc=temp("qc/trimmed/{sample}.{unit}.cutadapt_STEP1.qc.txt")
    log:
        "logs/trimmed/{sample}.{unit}.cutadapt_STEP1.log"
    wrapper:
        "0.17.4/bio/cutadapt/pe"

rule cutadapt_step2:
    input:
        "trimmed/{sample}.{unit}.tmpR2.fastq",
         "trimmed/{sample}.{unit}.tmpR1.fastq"
    output:
        fastq1=temp("trimmed/{sample}.{unit}.5ptmpR2.fastq"),
        fastq2=temp("trimmed/{sample}.{unit}.5ptmpR1.fastq"),
        qc=temp("qc/trimmed/{sample}.{unit}.cutadapt_STEP2.qc.txt")
    log:
        "logs/trimmed/{sample}.{unit}.cutadapt_STEP2.log"
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -g file:" + config["accel_panels"][samples["panel"][wildcards.sample]]["5p_primer_file"]
    wrapper:
        "0.17.4/bio/cutadapt/pe"

################################################################################
# Step 1 and 2
# anchored 3’ trimming of primer sequences with 3 primer design file.
###############################################################################

rule cutadapt_step3:
    input:
        ["trimmed/{sample}.{unit}.5ptmpR1.fastq",
         "trimmed/{sample}.{unit}.5ptmpR2.fastq"]
    output:
        fastq1=temp("trimmed/{sample}.{unit}.tmp3R1.fastq"),
        fastq2=temp("trimmed/{sample}.{unit}.tmp3R2.fastq"),
        qc=temp("qc/trimmed/{sample}.{unit}.cutadapt_STEP3.qc.txt")
    log:
        "logs/trimmed/{sample}.{unit}.cutadapt_STEP3.log"
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -a file:" + config["accel_panels"][samples["panel"][wildcards.sample]]["3p_primer_file"]
    wrapper:
        "0.17.4/bio/cutadapt/pe"

rule cutadapt_step4:
    input:
        "trimmed/{sample}.{unit}.tmp3R2.fastq",
         "trimmed/{sample}.{unit}.tmp3R1.fastq"
    output:
        fastq1="trimmed/{sample}.{unit}.R2.trimmomatic_cutadapt.fastq.gz",
        fastq2="trimmed/{sample}.{unit}.R1.trimmomatic_cutadapt.fastq.gz",
        qc=temp("qc/trimmed/{sample}.{unit}.cutadapt_STEP4.qc.txt")
    log:
        "logs/trimmed/{sample}.{unit}.cutadapt_STEP4.log"
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -a file:" + config["accel_panels"][samples["panel"][wildcards.sample]]["3p_primer_file"]
    wrapper:
        "0.17.4/bio/cutadapt/pe"

################################################################################
# Final step
# Merge all generate log files into one,
###############################################################################

rule merge_qc:
    input:
        qc=expand("qc/trimmed/{{sample}}.{{unit}}.{steps}.qc.txt",
                    steps=["cutadapt_STEP1","cutadapt_STEP2","cutadapt_STEP3","cutadapt_STEP4"])
    output:
         qc="qc/trimmed/{sample}.{unit}.trimmomatic_cutadapt.qc.txt"
    run:
        with open(output.qc,"w") as out:
            for qc_file in input.qc:
                with open(qc_file,"r") as qc_input:
                    out.write(qc_input.read())
