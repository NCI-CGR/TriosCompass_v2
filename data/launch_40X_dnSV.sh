#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH -o ${PWD}/snakemake.%j.out
#SBATCH -e ${PWD}/snakemake.%j.err

### or make sure singularity is available in alternative way
module load singularity

mkdir -p TMP
export TMPDIR=TMP

snakemake --profile TriosCompass_v2/workflow/profiles/slurm --configfile GIAB_40X_dnSV.yaml \
          --conda-frontend mamba \
          --verbose --skip-script-cleanup -k  --keep-incomplete --ignore-incomplete --latency-wait 600 
        #   \
        #   --groups vizaln=group0 extract_parental_origin=group1 \
        #   --group-components group0=20 group1=20 