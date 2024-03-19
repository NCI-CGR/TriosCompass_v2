#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH -o ${PWD}/logs/snakemake.%j.out
#SBATCH -e ${PWD}/logs/snakemake.%j.err


# conda activate snakemake
mkdir -p /data/DCEG_Trios/ChernobylTrios/output/TMP
export TMPDIR=/data/DCEG_Trios/ChernobylTrios/output/TMP

mkdir -p logs

snakemake --skip-script-cleanup -k  --keep-incomplete --rerun-incomplete --profile biowulf --verbose -p --use-conda --jobs 400 --use-envmodules --latency-wait 600 -T 0 -s realign12bams.snakefile