#! /bin/bash

mkdir -p output/dragen_joint/t0321c1

# export PATH=/usr/local/slurm/bin:/usr/local/apps/samtools/1.19/bin:/data/zhuw10/conda_envs/snakemake/bin:/data/zhuw10/miniconda3/condabin:/usr/local/slurm/bin:/usr/local/bin:/usr/X11R6/bin:/usr/local/jdk/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/local/mysql/bin:/home/zhuw10/.local/bin:/home/zhuw10/bin:/opt/edico/bin
source /etc/profile.d/edico.sh

ulimit -n 65535

dragen_lic -f genome &> output/dragen_joint/t0321c1/dragen_log.start.log

# load reference for your analysis, for DNA, RNA, CNV, and HLA analysis, use
# dragen -l -r /staging/human/hg38

dragen -f             -r /staging/human/hg38             --enable-joint-genotyping true             --variant output/dragen_gvcf/SC501111/SC501111.hard-filtered.gvcf.gz --variant output/dragen_gvcf/SC501096/SC501096.hard-filtered.gvcf.gz --variant output/dragen_gvcf/SC501108/SC501108.hard-filtered.gvcf.gz             --output-directory "output/dragen_joint/t0321c1"             --output-file-prefix t0321c1             --pedigree-file ped/t0321c1.ped

dragen_lic -f genome &> output/dragen_joint/t0321c1/dragen_log.end.log
