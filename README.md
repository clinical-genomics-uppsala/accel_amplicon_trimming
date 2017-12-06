# Snakemake workflow: accel-amplicon-trimming

[![Snakemake](https://img.shields.io/badge/snakemake-≥4.3.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/clinical-genomics-uppsala/accel_amplicon_trimming.svg?branch=master)](https://travis-ci.org/clinical-genomics-uppsala/accel_amplicon_trimming)


This workflow performs read trimming on Accel Amplicon Panel data, using the
recommended guidelines provided by [Swift](https://swiftbiosci.com/wp-content/uploads/2017/03/17-1397-Amplicon-Bioinf-Guidelines.pdf)

## Authors

Patrik Smeds (@smeds)

## Requirements

Samples have to be Paired-End Illumina sequences.

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the latest release. If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, once available, its DOI.

### Step 2: Configure workflow

Configure the workflow according to your needs by editing the config.yaml and the sample sheet samples.tsv.

#### Config.yaml
The following entries are expected in the config file
```
illuminaclip_file: /path/to/illumina.fa

accel_panels:
  panel1:
    5p_primer_file: /path/to/5p_primers.fa
    3p_primer_file: /path/to/3p_primers.fa
  panel2:
    5p_primer_file: /path/to/5p_primers.fa
    3p_primer_file: /path/to/3p_primers.fa
```

#### Sample.tsv
Example of a sample.tsv file (columns need to be
tab separated)
```
sample     panel
sample1    panel1
sample2    panel2
```
The panel column must contain a panel name that can be found
in the accel_panels entry in the config.yaml

#### unit.tsv
Example of a units.tsv file (columns need to be
tab separated)
```
sample     panel     fq1                fq2
sample1    panel1    /path/to/sample1.R1.fastq   /path/to/sample1.R2.fastq
sample2    panel2    /path/to/sample2.R1.fastq   /path/to/sample2.R2.fastq
```
sample	unit	fq1	fq2
A	lane1	A.1.fastq	A.2.fastq
```

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

```
snakemake -n
```
Execute the workflow locally via

```
snakemake --cores $N
```
using $N cores or run it in a cluster environment via

```
snakemake --cluster qsub --jobs 100
```
or

```
snakemake --drmaa --jobs 100
```
See the Snakemake documentation for further details.
