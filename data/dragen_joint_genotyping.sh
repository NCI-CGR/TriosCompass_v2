#! /bin/bash

# Input parameters
# 3 bam files and pedigree file and the output directory, ID

BAM1=$1
BAM2=$2
BAM3=$3
PEDIGREE=$4
OUTPUT_DIR=$5
ID=$6

# set up paths etc
export PATH="/opt/dragen/4.4.4/bin:$PATH"
genome=/staging/ref/current/hg38-alt_masked.cnv.graph.hla.methyl_cg.rna-11-r5.0-1
annotation=/staging/ref/current/annotations/hg38_all

# set up paths etc
RESULTPATH=${OUTPUT_DIR}/${ID}

mkdir -p $RESULTPATH

source /etc/profile.d/edico.sh
ulimit -n 65535
ulimit -u 16384

dragen_lic -f genome &> ${RESULTPATH}/dragen_log.start.log

# dragen -l -r $genome

dragen -f \
    -r $genome \
    --enable-joint-genotyping true \
    --variant $BAM1 \
    --variant $BAM2 \
    --variant $BAM3 \
    --output-directory "$RESULTPATH" \
    --output-file-prefix ${ID} \
    --pedigree-file $PEDIGREE

dragen_lic -f genome &> ${RESULTPATH}/dragen_log.end.log
