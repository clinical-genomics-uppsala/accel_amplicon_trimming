# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

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

Expects a samples config
With the following dictionary structure

samples = {
    'panel': {sample_id: value},
    'input_folder': {samplei_id: value},
    'sample_number' {sample_id: value},
    'lane': {sample_id: value}
}
the panel key should contain a value that can be found under
accel_panels in the config.yaml file.

Example of a sample.tsv that can be imported using pandas
...............................................................................
sample     panel     input_folder    sample_number  lane
sample1    panel1    fastq/sample1   S1             L001
sample2    panel2    fastq           S2             L001
...............................................................................
"""

__author__ = "Patrik Smeds"
__license__ = "MIT"

###############################################################################
# Step 1:
# Remove any Illummina adaptor sequence using Trimmomatic
# (http://www.usadellab.org/cms/?page=trimmomatic).
###############################################################################

#Function used to generate input
def get_fastq_file(wildcards, read_pair="R1", compressed= False):
    if samples is None or \
        samples.get('input_folder',{}).get(wildcards.sample, None) is None:
        return wildcards.sample + "_" + \
               wildcards.sample_number + "_" + \
               wildcards.lane + "_"  + \
               read_pair + "_" + \
               wildcards.counter + ".fastq" + ("" if not compressed else ".gz")
    else:
        return os.path.join(samples['input_folder'][wildcards.sample],
           wildcards.sample + "_" +
           wildcards.sample_number + "_" +
           wildcards.lane + "_"  +
           read_pair + "_" +
           wildcards.counter + ".fastq" + ("" if not compressed else ".gz"))



# If the fastq files already have been decompressed use those files, to not
# waste calculations on decompressing the files once more.
ruleorder: trimmomatic > trimmomatic_compressed

# Rule to perform trimmomatic operations on none compressed fastq files.
rule trimmomatic:
    input:
        r1 = lambda wildcards: get_fastq_file(wildcards, "R1"),
        r2 = lambda wildcards: get_fastq_file(wildcards, "R2")
    output:
        temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.trimmomatic.qc.txt"),
        r1 = temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_R1_{counter,\d+}.trimmomatic.fastq"),
        r2 = temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_R2_{counter,\d+}.trimmomatic.fastq"),
        r1_unpaired=temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_R1_{counter,\d+}.trimmomatic.up.fastq"),
        r2_unpaired=temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_R2_{counter,\d+}.trimmomatic.up.fastq")
    threads: 12
    log:
        "{sample}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.trimmomatic.qc.txt"
    params:
        extra="-threads 12",
        trimmer=["ILLUMINACLIP:" + config["illuminaclip_file"] + ":2:30:10", "MINLEN:30"]
    wrapper:
        "0.17.4/bio/trimmomatic/pe"

# Rule to perform trimmomatic operations on compressed fastq files.
rule trimmomatic_compressed:
    input:
        r1 = lambda wildcards: get_fastq_file(wildcards, "R1", compressed=True),
        r2 = lambda wildcards: get_fastq_file(wildcards, "R2", compressed=True)
    output:
        temp("temp{sample}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.trimmomatic.qc.txt"),
        r1 = temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_R1_{counter,\d+}.trimmomatic.fastq"),
        r2 = temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_R2_{counter,\d+}.trimmomatic.fastq"),
        r1_unpaired=temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_R1_{counter,\d+}.trimmomatic.up.fastq"),
        r2_unpaired=temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_R2_{counter,\d+}.trimmomatic.up.fastq")
    threads: 12
    log:
        temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.trimmomatic.qc.txt")
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
        ["{sample}_{sample_number}_{lane}_R1_{counter}.trimmomatic.fastq",
         "{sample}_{sample_number}_{lane}_R2_{counter}.trimmomatic.fastq"]
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -g file:" + config["accel_panels"][samples["panel"][wildcards.sample]]["5p_primer_file"]
    output:
        fastq1=temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.tmpR1.fastq"),
        fastq2=temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.tmpR2.fastq"),
        qc=temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.cutadapt_STEP1.qc.txt")
    wrapper:
        "0.17.4/bio/cutadapt/pe"

rule cutadapt_step2:
    input:
        ["{sample}_{sample_number}_{lane}_{counter}.tmpR2.fastq",
         "{sample}_{sample_number}_{lane}_{counter}.tmpR1.fastq"]
    output:
        fastq1=temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.5ptmpR2.fastq"),
        fastq2=temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.5ptmpR1.fastq"),
        qc = temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.cutadapt_STEP2.qc.txt")
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
        ["{sample}_{sample_number}_{lane}_{counter}.5ptmpR1.fastq",
         "{sample}_{sample_number}_{lane}_{counter}.5ptmpR2.fastq"]
    output:
        fastq1=temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.tmp3R1.fastq"),
        fastq2=temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.tmp3R2.fastq"),
        qc = temp("{sample}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.cutadapt_STEP3.qc.txt")
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -a file:" + config["accel_panels"][samples["panel"][wildcards.sample]]["3p_primer_file"]
    wrapper:
        "0.17.4/bio/cutadapt/pe"

ruleorder: cutadapt_step4_with_outdirectory > cutadapt_step4

rule cutadapt_step4:
    input:
        ["{sample}_{sample_number}_{lane}_{counter}.tmp3R2.fastq",
         "{sample}_{sample_number}_{lane}_{counter}.tmp3R1.fastq"]
    output:
        fastq1="{sample,[A-Za-z0-9-]+}_{sample_number,S\d+}_{lane,L\d+}_R2_{counter,\d+}.trimmomatic_cutadapt.fastq.gz",
        fastq2="{sample,[A-Za-z0-9-]+}_{sample_number,S\d+}_{lane,L\d+}_R1_{counter,\d+}.trimmomatic_cutadapt.fastq.gz",
        qc = temp("{sample,[A-Za-z0-9-]+}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.cutadapt_STEP4.qc.txt")
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -a file:" + config["accel_panels"][samples["panel"][wildcards.sample]]["3p_primer_file"]
    wrapper:
        "0.17.4/bio/cutadapt/pe"

rule cutadapt_step4_with_outdirectory:
    input:
        ["{sample}_{sample_number}_{lane}_{counter}.tmp3R2.fastq",
         "{sample}_{sample_number}_{lane}_{counter}.tmp3R1.fastq"]
    output:
        fastq1="{output_dir}/{sample,[A-Za-z0-9-]+}_{sample_number,S\d+}_{lane,L\d+}_R2_{counter,\d+}.trimmomatic_cutadapt.fastq.gz",
        fastq2="{output_dir}/{sample,[A-Za-z0-9-]+}_{sample_number,S\d+}_{lane,L\d+}_R1_{counter,\d+}.trimmomatic_cutadapt.fastq.gz",
        qc = temp("{output_dir}/{sample,[A-Za-z0-9-]+}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.cutadapt_STEP4.qc.txt")
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

rule merge_logs:
    input:
        qc = expand("{{sample}}_{{sample_number}}_{{lane}}_{{counter}}.{steps}.qc.txt",
                    steps=["trimmomatic","cutadapt_STEP1","cutadapt_STEP2","cutadapt_STEP3","cutadapt_STEP4"])
    output:
         qc="{sample,[A-Za-z0-9-]+}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.trimmomatic_cutadapt.qc.txt"
    run:
        with open(output.qc,"w") as out:
            for qc_file in input.qc:
                with open(qc_file,"r") as qc_input:
                    out.write(str(qc_input.read()))

def generate_input_rule_path(wildcards, steps, pre_path = None):
    def generate_path(wildcards, step):
        return "_".join([wildcards.sample, wildcards.sample_number, wildcards.lane, wildcards.counter]) + "." + step + ".qc.txt"
    if pre_path is None:
        return [ generate_path(wildcards,end) for end in steps]
    else:
        return [ os.path.join(pre_path, generate_path(wildcards,end)) for end in steps]

rule merge_logs_with_outdirectory:
    input:
        qc=lambda wildcards:
            generate_input_rule_path(wildcards, ["trimmomatic","cutadapt_STEP1","cutadapt_STEP2","cutadapt_STEP3"]) +
            generate_input_rule_path(wildcards, ["cutadapt_STEP4"], wildcards.output_dir )
    output:
         qc="{output_dir}/{sample,[A-Za-z0-9-]+}_{sample_number,S\d+}_{lane,L\d+}_{counter,\d+}.trimmomatic_cutadapt.qc.txt"
    run:
        with open(output.qc,"w") as out:
            for qc_file in input.qc:
                with open(qc_file,"r") as qc_input:
                    out.write(str(qc_input.read()))
