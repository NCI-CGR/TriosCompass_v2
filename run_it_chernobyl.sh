#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH -o ${PWD}/snakemake.%j.out
#SBATCH -e ${PWD}/snakemake.%j.err


# conda activate snakemake
# export TMPDIR=/data/DCEG_Chernobyl/NP0436-HE5/GIAB_DATA/TMP
mkdir -p TMP
export TMPDIR=TMP

module load singularity 

# -k to keep going even on errors
snakemake --verbose --skip-script-cleanup -k  --keep-incomplete --ignore-incomplete --profile workflow/profiles/biowulf --verbose -p --use-conda --jobs 800 --use-singularity --use-envmodules --latency-wait 120 -T 0  -s Snakefile_chernobyl

 
