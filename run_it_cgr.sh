#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH -o ${PWD}/snakemake.%j.out
#SBATCH -e ${PWD}/snakemake.%j.err


# conda activate snakemake
export TMPDIR=/data/DCEG_Trios/new_cgr_data/TriosCompass_v2/TMP

snakemake --skip-script-cleanup -k  --keep-incomplete --rerun-incomplete --profile workflow/profiles/biowulf --verbose -p --use-conda --jobs 400 --use-envmodules --latency-wait 600 -T 0 -s Snakefile

