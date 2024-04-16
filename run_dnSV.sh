#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH -o ${PWD}/snakemake.%j.out
#SBATCH -e ${PWD}/snakemake.%j.err


# conda activate snakemake
mkdir -p TMP
export TMPDIR=TMP

export REF_PATH=/data/DCEG_Trios/new_cgr_data/TriosCompass_v2/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
export REF_CACHE=/data/DCEG_Trios/new_cgr_data/TriosCompass_v2/cache/%2s/%2s/%s

module load singularity 

snakemake --skip-script-cleanup -k  --keep-incomplete --rerun-incomplete --profile workflow/profiles/biowulf --verbose -p --use-conda --jobs 400 --use-singularity --singularity-args " -B /vf,/spin1,/data,/fdb,/gpfs " --use-envmodules --latency-wait 600 -T 0 -s Snakefile_dnSV