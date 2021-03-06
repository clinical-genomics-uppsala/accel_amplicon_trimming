# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Patrik Smeds"
__license__ = "MIT"

###############################################################################
# Step 1:
# Remove any Illummina adaptor sequence using Trimmomatic
# (http://www.usadellab.org/cms/?page=trimmomatic).
###############################################################################
from pytools.persistent_dict import PersistentDict

storage = PersistentDict("accelamplicon_storage")

def _accel_get_num_splits(config):
    return int(config.get("cgu_accel_num_fastq_split",config.get("num_fastq_split",1)))

def _accel_get_fastq(wildcards,units,read):
    return units.loc[(wildcards.sample, wildcards.unit), [read]].dropna()[0]

rule cgu_accel_extract_fastq_files:
   input:
      lambda wildcards: _accel_get_fastq(wildcards,units, 'fq1' if wildcards.read == "R1" else 'fq2')
   output:
      temp("trimmed/.temp_accel/{sample}.{unit}.{read,[R12]+}.fastq")
   run:
      if input[0].endswith("gz"):
         shell("zcat {input} > {output}")
      else:
         shell("cat {input} > {output}")

rule cgu_accel_count_lines_in_fastq:
    input:
      "trimmed/.temp_accel/{sample}.{unit}.{read}.fastq"
    output:
      temp("trimmed/.temp_accel/{sample}.{unit}.{read}.var")
    wildcard_constraints:
      sample="[A-Za-z0-9-_]+",
      unit="[A-Za-z0-9]+",
      read="[R12]+"
    run:
      import subprocess, os
      lines = int(float(subprocess.run("wc -l " + str(input[0]) + " |  awk '{print($1/4)}'", stdout=subprocess.PIPE,shell=True).stdout.decode('utf-8').rstrip("\n")))
      storage.store(wildcards.sample + "." + wildcards.unit + "." + wildcards.read + ".var",str(lines))
      shell("echo 'reads: '" + str(lines) + "'' > "  + output[0])

rule cgu_accel_split_fastq_file:
   input:
      "trimmed/.temp_accel/{sample}.{unit}.{read}.fastq",
      "trimmed/.temp_accel/{sample}.{unit}.{read}.var"
   output:
      temp(['trimmed/.temp_accel/{sample}.{unit}.%04d.{read}.fastq' % num for num in range(0,_accel_get_num_splits(config))])
   wildcard_constraints:
      sample="[A-Za-z0-9-_]+",
      unit="[A-Za-z0-9]+",
      read="[R12]+"
   params:
      output_prefix=lambda wildcards: "trimmed/.temp_accel/" + wildcards.sample + "." + wildcards.unit + ".",
      output_suffix=lambda wildcards: "." + wildcards.read + ".fastq"
   run:
     import math
     num_reads = int(storage.fetch(wildcards.sample + "." + wildcards.unit + "." + wildcards.read + ".var"))
     num_split = _accel_get_num_splits(config)
     lines_per_file = 4*math.ceil(num_reads / num_split)
     shell('cat {input[0]} | awk \'BEGIN{{ file = 0; filename = sprintf("{params.output_prefix}%.04d{params.output_suffix}", file) }}{{ print > filename}} NR % {lines_per_file} == 0 {{ close(filename); file = file + 1; filename = sprintf("{params.output_prefix}%.04d{params.output_suffix}",file)}}\'')
     num_files_generated = 4*math.floor(num_reads / lines_per_file)
     while num_files_generated < num_split:
        shell("touch {params.output_prefix}%04d{params.output_suffix}" % num_split)
        num_split -= 1

rule cgu_accel_trimmomatic:
    input:
        r1 = "trimmed/.temp_accel/{sample}.{unit}.{part}.R1.fastq",
        r2 = "trimmed/.temp_accel/{sample}.{unit}.{part}.R2.fastq"
    output:
        temp("logs/trimmed/.temp_accel/{sample}.{unit}.{part}.trimmomatic.qc.txt"),
        r1 = temp("trimmed/.temp_accel/{sample}.{unit}.{part}.R1.trimmomatic.fastq"),
        r2 = temp("trimmed/.temp_accel/{sample}.{unit}.{part}.R2.trimmomatic.fastq"),
        r1_unpaired=temp("trimmed/.temp_accel/{sample}.{unit}.{part}.R1.trimmomatic.up.fastq"),
        r2_unpaired=temp("trimmed/.temp_accel/{sample}.{unit}.{part}.R2.trimmomatic.up.fastq")
    wildcard_constraints:
        sample="[A-Za-z0-9-_]+",
        unit="[A-Za-z0-9]+",
        part="[0-9]+"
    threads: 8
    log:
        "logs/trimmed/.temp_accel/{sample}.{unit}.{part}.trimmomatic.qc.txt"
    params:
    #ToDo fix so that threads are configurable!!!
        extra=lambda wildcards: "-threads 8 " + config.get("phread_flag",""),
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

rule cgu_accel_cutadapt_step1:
    input:
        "trimmed/.temp_accel/{sample}.{unit}.{part}.R1.trimmomatic.fastq",
         "trimmed/.temp_accel/{sample}.{unit}.{part}.R2.trimmomatic.fastq"
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -g file:" + config["cgu_accel_panels"][samples["panel"][wildcards.sample]]["5p_primer_file"]
    output:
        fastq1=temp("trimmed/.temp_accel/{sample}.{unit}.{part}.tmpR1.fastq"),
        fastq2=temp("trimmed/.temp_accel/{sample}.{unit}.{part}.tmpR2.fastq"),
        qc=temp("qc/trimmed/.temp_accel/{sample}.{unit}.{part}.cutadapt_STEP1.qc.txt")
    wildcard_constraints:
        sample="[A-Za-z0-9-_]+",
        unit="[A-Za-z0-9]+",
        part="[0-9]+"
    log:
        "logs/trimmed/{sample}.{unit}.{part}.cutadapt_STEP1.log"
    wrapper:
        "0.17.4/bio/cutadapt/pe"

rule cgu_accel_cutadapt_step2:
    input:
        "trimmed/.temp_accel/{sample}.{unit}.{part}.tmpR2.fastq",
         "trimmed/.temp_accel/{sample}.{unit}.{part}.tmpR1.fastq"
    output:
        fastq1=temp("trimmed/.temp_accel/{sample}.{unit}.{part}.5ptmpR2.fastq"),
        fastq2=temp("trimmed/.temp_accel/{sample}.{unit}.{part}.5ptmpR1.fastq"),
        qc=temp("qc/trimmed/.temp_accel/{sample}.{unit}.{part}.cutadapt_STEP2.qc.txt")
    wildcard_constraints:
        sample="[A-Za-z0-9-_]+",
        unit="[A-Za-z0-9]+",
        part="[0-9]+"
    log:
        "logs/trimmed/.temp_accel/{sample}.{unit}.{part}.cutadapt_STEP2.log"
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -g file:" + config["cgu_accel_panels"][samples["panel"][wildcards.sample]]["5p_primer_file"]
    wrapper:
        "0.17.4/bio/cutadapt/pe"

################################################################################
# Step 1 and 2
# anchored 3’ trimming of primer sequences with 3 primer design file.
###############################################################################

rule cgu_accel_cutadapt_step3:
    input:
        ["trimmed/.temp_accel/{sample}.{unit}.{part}.5ptmpR1.fastq",
         "trimmed/.temp_accel/{sample}.{unit}.{part}.5ptmpR2.fastq"]
    output:
        fastq1=temp("trimmed/.temp_accel/{sample}.{unit}.{part}.tmp3R1.fastq"),
        fastq2=temp("trimmed/.temp_accel/{sample}.{unit}.{part}.tmp3R2.fastq"),
        qc=temp("qc/trimmed/.temp_accel/{sample}.{unit}.{part}.cutadapt_STEP3.qc.txt")
    wildcard_constraints:
        sample="[A-Za-z0-9-_]+",
        unit="[A-Za-z0-9]+",
        part="[0-9]+"
    log:
        "logs/trimmed/{sample}.{unit}.{part}.cutadapt_STEP3.log"
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -a file:" + config["cgu_accel_panels"][samples["panel"][wildcards.sample]]["3p_primer_file"]
    wrapper:
        "0.17.4/bio/cutadapt/pe"

rule cgu_accel_cutadapt_step4:
    input:
        "trimmed/.temp_accel/{sample}.{unit}.{part}.tmp3R2.fastq",
         "trimmed/.temp_accel/{sample}.{unit}.{part}.tmp3R1.fastq"
    output:
        fastq1=temp("trimmed/.temp_accel/{sample}.{unit}.{part}.R2.trimmomatic_cutadapt.fastq.gz"),
        fastq2=temp("trimmed/.temp_accel/{sample}.{unit}.{part}.R1.trimmomatic_cutadapt.fastq.gz"),
        qc=temp("qc/trimmed/.temp_accel/{sample}.{unit}.{part}.cutadapt_STEP4.qc.txt")
    wildcard_constraints:
        sample="[A-Za-z0-9-_]+",
        unit="[A-Za-z0-9]+",
        part="[0-9]+"
    log:
        "logs/trimmed/{sample}.{unit}.{part}.cutadapt_STEP4.log"
    params:
        " --minimum-length 40",
        " -e 0.12",
        lambda wildcards: \
            " -a file:" + config["cgu_accel_panels"][samples["panel"][wildcards.sample]]["3p_primer_file"]
    wrapper:
        "0.17.4/bio/cutadapt/pe"

################################################################################
# Final step
# Merge all generate log files into one,
###############################################################################

def get_parts(config):
  return [ "%04d"  % part for part in range(0,_accel_get_num_splits(config)) for unit in units]

rule cgu_accel_merge_qc_split:
    input:
        qc=lambda wildcards: ["qc/trimmed/.temp_accel/" + wildcards.sample + "." + wildcards.unit + "." + part + "." + wildcards.step + ".qc.txt" for part in get_parts(config)]
    output:
         qc="qc/trimmed/.temp_accel/{sample}.{unit}.{step}.qc.txt"
    wildcard_constraints:
        sample="[A-Za-z0-9-_]+",
        unit="[A-Za-z0-9]+",
        step="[A-Za-z0-9_]+",
        part="[0-9]+"
    run:
        with open(output.qc,"w") as out:
            for qc_file in input.qc:
                with open(qc_file,"r") as qc_input:
                    out.write(qc_input.read())

rule cgu_accel_merge_qc_final:
    input:
        qc=expand("qc/trimmed/.temp_accel/{{sample}}.{{unit}}.{steps}.qc.txt",
                    steps=["cutadapt_STEP1","cutadapt_STEP2","cutadapt_STEP3","cutadapt_STEP4"])
    output:
         qc="qc/trimmed/{sample}.{unit}.trimmomatic_cutadapt.qc.txt"
    wildcard_constraints:
         sample="[A-Za-z0-9-_]+",
         unit="[A-Za-z0-9]+",
         step="[A-Za-z0-9_]+",
         part="[0-9]+"
    run:
        with open(output.qc,"w") as out:
            for qc_file in input.qc:
                with open(qc_file,"r") as qc_input:
                    out.write(qc_input.read())

rule cgu_accel_merge_qc_final_split:
    input:
        qc=expand("qc/trimmed/.temp_accel/{{sample}}.{{unit}}.{{part}}.{steps}.qc.txt",
                    steps=["cutadapt_STEP1","cutadapt_STEP2","cutadapt_STEP3","cutadapt_STEP4"])
    output:
         qc="qc/trimmed/{sample}.{unit}.{part}.trimmomatic_cutadapt.qc.txt"
    wildcard_constraints:
         sample="[A-Za-z0-9-_]+",
         unit="[A-Za-z0-9]+",
         step="[A-Za-z0-9_]+",
         part="[0-9]+"
    run:
        with open(output.qc,"w") as out:
            for qc_file in input.qc:
                with open(qc_file,"r") as qc_input:
                    out.write(qc_input.read())

rule cgu_accel_move_fastq:
    input:
        "trimmed/.temp_accel/{sample}.{unit}.{part}.{read}.trimmomatic_cutadapt.fastq.gz"
    output:
        "trimmed/{sample}.{unit}.{part}.{read}.trimmomatic_cutadapt.fastq.gz"
    wildcard_constraints:
         sample="[A-Za-z0-9-_]+",
         unit="[A-Za-z0-9]+",
         step="[A-Za-z0-9_]+",
         part="[0-9]+",
         read="[R12]+"
    shell: "mv {input} {output}"

rule cgu_accel_merge_split:
    input:
        lambda wildcards: ["trimmed/.temp_accel/" + wildcards.sample + "." + wildcards.unit + "." + part + "." + wildcards.read + ".trimmomatic_cutadapt.fastq.gz" for part in get_parts(config)]
    output:
        "trimmed/{sample}.{unit}.{read}.trimmomatic_cutadapt.fastq.gz"
    wildcard_constraints:
         sample="[A-Za-z0-9-_]+",
         unit="[A-Za-z0-9]+",
         step="[A-Za-z0-9_]+",
         part="[0-9]+"
    shell: "zcat {input} | gzip > {output}"
