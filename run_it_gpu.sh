#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH -o ${PWD}/snakemake.%j.out
#SBATCH -e ${PWD}/snakemake.%j.err


# conda activate snakemake
export TMPDIR=/data/DCEG_Chernobyl/NP0436-HE5/GIAB_DATA/TMP

snakemake --profile workflow/profile/biowulf --verbose -p --use-conda --jobs 400 --use-envmodules --latency-wait 120 -T 0 

