#! /bin/bash

# set up paths etc
source /etc/profile.d/edico.sh

source /etc/profile.d/edico.sh

dragen_lic -f genome &> output/dragen_gvcf/SC501108/dragen_log.start.log

#load reference for your analysis, for DNA, RNA, CNV, and HLA analysis, use
dragen -l -r /staging/human/hg38

mkdir -p output/dragen_gvcf/SC501108

dragen -f             -r /staging/human/hg38             --vc-target-bed ./ref/chr22.bed --fastq-list output/fq_list/SC501108.fastq_list.csv             --fastq-list-sample-id SC501108             --enable-variant-caller true             --enable-map-align-output true             --enable-map-align true             --enable-duplicate-marking true             --vc-emit-ref-confidence GVCF             --vc-enable-vcf-output true             --output-directory "$(dirname output/dragen_gvcf/SC501108/SC501108.gvcf.gz)"             --output-file-prefix SC501108             --output-format cram             --enable-cnv true             --cnv-enable-self-normalization true --cnv-interval-width 1000              --enable-sv true


dragen_lic -f genome &> output/dragen_gvcf/SC501108/dragen_log.end.log
    