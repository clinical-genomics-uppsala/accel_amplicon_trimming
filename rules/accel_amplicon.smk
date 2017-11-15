# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Patrik Smeds"
__license__ = "MIT"

###############################################################################
# Step 1:
# Remove any Illummina adaptor sequence using Trimmomatic
# (http://www.usadellab.org/cms/?page=trimmomatic).
###############################################################################
def _sample_name(wildcards):
    return os.path.split(wildcards.sample)[-1]

def _get_fastq(wildcards,units,read_pair='fq1'):
    return units.loc[(_sample_name(wildcards), wildcards.unit), [read_pair]].dropna()[0]

# Rule to perform trimmomatic operations on none compressed fastq files.
rule trimmomatic:
    input:
        r1 = lambda wildcards: _get_fastq(wildcards, units, 'fq1'),
        r2 = lambda wildcards: _get_fastq(wildcards, units, 'fq2')
    output:
        temp("{sample}-{unit}.trimmomatic.qc.txt"),
        r1 = temp("{sample}-{unit}.R1.trimmomatic.fastq"),
        r2 = temp("{sample}-{unit}.R2.trimmomatic.fastq"),
        r1_unpaired=temp("{sample}-{unit}.R1.trimmomatic.up.fastq"),
        r2_unpaired=temp("{sample}-{unit}.R2.trimmomatic.up.fastq")
    threads: 12
    log:
        "{sample}-{unit}.trimmomatic.qc.txt"
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
        ["{sample}-{unit}.R1.trimmomatic.fastq",
         "{sample}-{unit}.R2.trimmomatic.fastq"]
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -g file:" + config["accel_panels"][samples["panel"][_sample_name(wildcards)]]["5p_primer_file"]
    output:
        fastq1=temp("{sample}-{unit}.tmpR1.fastq"),
        fastq2=temp("{sample}-{unit}.tmpR2.fastq"),
        qc=temp("{sample}-{unit}.cutadapt_STEP1.qc.txt")
    wrapper:
        "0.17.4/bio/cutadapt/pe"

rule cutadapt_step2:
    input:
        ["{sample}-{unit}.tmpR2.fastq",
         "{sample}-{unit}.tmpR1.fastq"]
    output:
        fastq1=temp("{sample}-{unit}.5ptmpR2.fastq"),
        fastq2=temp("{sample}-{unit}.5ptmpR1.fastq"),
        qc = temp("{sample}-{unit}.cutadapt_STEP2.qc.txt")
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -g file:" + config["accel_panels"][samples["panel"][_sample_name(wildcards)]]["5p_primer_file"]
    wrapper:
        "0.17.4/bio/cutadapt/pe"

################################################################################
# Step 1 and 2
# anchored 3’ trimming of primer sequences with 3 primer design file.
###############################################################################

rule cutadapt_step3:
    input:
        ["{sample}-{unit}.5ptmpR1.fastq",
         "{sample}-{unit}.5ptmpR2.fastq"]
    output:
        fastq1=temp("{sample}-{unit}.tmp3R1.fastq"),
        fastq2=temp("{sample}-{unit}.tmp3R2.fastq"),
        qc = temp("{sample}-{unit}.cutadapt_STEP3.qc.txt")
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -a file:" + config["accel_panels"][samples["panel"][_sample_name(wildcards)]]["3p_primer_file"]
    wrapper:
        "0.17.4/bio/cutadapt/pe"

rule cutadapt_step4:
    input:
        ["{sample}-{unit}.tmp3R2.fastq",
         "{sample}-{unit}.tmp3R1.fastq"]
    output:
        fastq1="{sample}-{unit}.R1.trimmomatic_cutadapt.fastq.gz",
        fastq2="{sample}-{unit}.R2.trimmomatic_cutadapt.fastq.gz",
        qc = temp("{sample}-{unit}.cutadapt_STEP4.qc.txt")
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -a file:" + config["accel_panels"][samples["panel"][_sample_name(wildcards)]]["3p_primer_file"]
    wrapper:
        "0.17.4/bio/cutadapt/pe"

################################################################################
# Final step
# Merge all generate log files into one,
###############################################################################

rule merge_logs:
    input:
        qc=expand("{{sample}}-{{unit}}.{steps}.qc.txt",
                    steps=["trimmomatic","cutadapt_STEP1","cutadapt_STEP2","cutadapt_STEP3","cutadapt_STEP4"])
    output:
         qc="{sample}-{unit}.trimmomatic_cutadapt.qc.txt"
    run:
        with open(output.qc,"w") as out:
            for qc_file in input.qc:
                with open(qc_file,"r") as qc_input:
                    out.write(str(qc_input.read()))
