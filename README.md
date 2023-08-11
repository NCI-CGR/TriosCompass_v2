# TriosCompass_v2
This is a trios analysis workflow written in Snakemake.

## Introduction

This workflow takes fastq files and pedigree files as inupt and generate lists of de novo mutations for each tios in the end. Three variant callers are used here: DeepVariant, GATK HaplotypeCaller and Strelka.  

## Installation

The workflow has been tested under biowulf with snakemake version 7.3.7.  Snakemake is installed under conda.

## Prepare input files
###  Manifest file and its schema specification

The manifest file is a csv file with typical NCI-CGR sample sheet format, with additional two columns (for the absolute paths to the paired Fastq files): R1, R2 (see https://github.com/NCI-CGR/TriosCompass_v2/blob/main/pep/manifest_fastq.csv).  The columns "CGR_ID", "INDEX", "FLOWCELL", "R1/2" are required and all the others are optional.  Of them, "CGR_ID", "INDEX", and "FLOWCELL" are used to define the read groups in the bam files, in the format:  
```
@RG\\tPL:ILLUMINA\\tID:{FLOWCELL}_{LANE}\\tSM:{CGF_ID}\\tPU:{CGF_ID}_{FLOWCELL}\\tLB:{CGF_ID}_{INDEX}
```

The format of the manifest file is specified by schemas/cgr_manifest_schema.yaml. The manifest file will be automatically validated at the beginning of the Snakemake workflow via the new Snakemake feature of PEP (protable encapsulated project).  Users may vist [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#configuring-scientific-experiments-via-peps) to learn more about this PEP feature. 

:bookmark: Note that one sample is allowed multiple fastq files from different flow cells.  This workflow will combine those fastq files and generate one bam file for each sample.



---

## Get started 


sbatch -J snakemake -t 200:00:00 --export=ALL --mem=12g -p norm  --wrap='./run_it_gpu.sh '
