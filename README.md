# Snakemake workflow: `<name>`

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.3.7-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for DNM (de novo mutation) calling.

## Introduction

This branch is to develop TriosCompass for general public use.  There are several major changes having been made: 
+ Remove the dependencies of the NIH biowulf modules.
  + Replace with either conda or docker container.
+ Break down the whole workflow into multiple sub-workflows to enhance management.
+ Move the hard-coded resource settings to two configure files:
  + workflow/profiles/biowulf/config.yaml
  + config/config.yaml (see the property ***threads***)


## Usage

Most Snakemake command-line option has also been configured at workflow/profiles/biowulf/config.yaml, so as to simplify the launch of the workflow.

### Dependencies
+ conda + mamba
+ Snakemake (Version 7.3.7+)
+ Singularity
+ python modules
  + peppy
  + peds
  + eido
  + uuid
  + itertools

### Resources
+ Homo_sapiens_assembly38.fasta
+ hg38.wgs_interval.bed
  
### Input files
+ fastq input files
+ pedigree files for each family

+ Example of the ped file
  + Special requirements:
    + The first line is for father, the second line is for mother, and the 3rd line (also the last line) is for child of the trio family.
    + Basename of the ped file is the concatenation of the family ID and the extension *.ped*.
```bash
cat ped/AJ.ped 
AJ      HG003   0       0       1       1
AJ      HG004   0       0       2       1
AJ      HG002   HG003   HG004   1       1
```

### Get started
```bash
cd /data/DCEG_Trios/ChernobylTrios/dev

git clone -b general_use git@github.com:NCI-CGR/TriosCompass_v2.git

conda activate snakemake
module load singularity

mkdir -p TMP
export TMPDIR=TMP

snakemake --configfile TriosCompass_v2/config/config.yaml --profile TriosCompass_v2/workflow/profiles/biowulf 2>&1 | tee snakemake.log 

```


### workflow/profiles/biowulf/config.yaml
:bookmark: Note that *set-threads* setting cannot work properly as expected, so that we took an alternative solution to put threads setting in [config/config.yaml](config/config.yaml).

+ workflow/profiles/biowulf/config.yaml
```yml

cluster-cancel: "scancel"
jobscript: "slurm-jobscript.sh"
cluster: "slurm-submit.py"
cluster-status: "slurm-status.py"


# Example resource configuration
# default-resources:
#   - runtime=100
#   - mem_mb=6000
#   - disk_mb=1000000
# # set-threads: map rule names to threads
# set-threads:
#   - single_core_rule=1
#   - multi_core_rule=10
# # set-resources: map rule names to resources in general
# set-resources:
#   - high_memory_rule:mem_mb=12000
#   - long_running_rule:runtime=1200


snakefile: TriosCompass_v2/workflow/Snakefile
verbose: True
skip-script-cleanup: True
latency-wait: 60
reason: True
show-failed-logs: True
keep-going: True
printshellcmds: True
rerun-incomplete: True
use-conda: True
use-envmodules: True
use-singularity: True
local-cores: 1
printshellcmds: True

restart-times: "0"

# Cluster submission
jobname: "{rule}.{jobid}"   
max-jobs-per-second: 1     
max-status-checks-per-second: 10     
jobs: 800     

# Job resources
set-resources:
  - bwa_index:mem_mb=16000
  - fastp:mem_mb=40000
  - fastqc:mem_mb=10000
  - fastq_screen:mem_mb=50000
  - fq2bam:mem_mb=40000
  - fq2bam:runtime="8h"
  - fq2bam:partition=gpu
  - fq2bam:slurm=gres=gpu:v100x:1 
  - flagstat:mem_mb=80000
  - collectwgsmetrics:runtime="1d"
  - collectwgsmetrics:mem_mb=80000
  - collectmultiplemetrics:mem_mb=80000
  - collectmultiplemetrics:partition=gpu
  - collectmultiplemetrics:slurm=gres=gpu:v100x:1 
  - gatkhc_pb:mem_mb=240000
  - gatkhc_pb:runtime="1d"
  - gatkhc_pb:partition=gpu
  - gatkhc_pb:slurm=gres=gpu:v100x:1
  - gatk_combine_gvcf:mem_mb=40000
  - gatk_combine_gvcf:runtime="6d"
  - gatk_cgp:mem_mb=40000
  - gatk_cgp:runtime="1d"
  - gatk_genotype_gvcf_pb:runtime="1d"
  - gatk_genotype_gvcf_pb:mem_mb=80000
  - gatk_genotype_gvcf_pb:partition=gpu
  - gatk_genotype_gvcf_pb:slurm=gres=gpu:v100x:1 
  - deepvariant_pb:mem_mb=180000
  - deepvariant_pb:runtime="1d"
  - deepvariant_pb:slurm=gres=gpu:v100x:1
  - deepvariant_pb:partition=gpu
  - glnexus_dv:mem_mb=100000
  - glnexus_dv:runtime="1d"
  - call_dnm_dv:mem_mb=20000
  - call_dnm_dv:runtime="1d"
  - call_dnm_hc:mem_mb=20000
  - call_dnm_hc:runtime="1d"
  - call_JIGV:mem_mb=60000


# Define the number of threads used by rules
# buggy about this setting parse (see 121M)
# set-threads:
#   - "bwa_index=16"

# For some reasons time needs quotes to be read by snakemake
default-resources:
  - mem_mb=2000
  - runtime="10:00:00"

```
