# Snakemake workflow: accel-amplicon-trimming

[![Snakemake](https://img.shields.io/badge/snakemake-≥4.3.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/clinical-genomics-uppsala/accel_amplicon_trimming.svg?branch=master)](https://travis-ci.org/clinical-genomics-uppsala/accel_amplicon_trimming)


This workflow performs read trimming on Accel Amplicon Panel data, using the
recommended guidelines provided by [Swift](https://swiftbiosci.com/wp-content/uploads/2017/03/17-1397-Amplicon-Bioinf-Guidelines.pdf)

## Authors

Patrik Smeds (@smeds)
## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the latest release. If you intend to modify and further develop this workflow, fork this reposity. Please consider providing any generally applicable modifications via a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, once available, its DOI.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file config.yaml and the sample sheet samples.tsv.

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
