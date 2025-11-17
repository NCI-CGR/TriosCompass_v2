#! /bin/bash

BAM=$1
ID=$2
OUTPUT_DIR=$3

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
--bam-input $BAM \
--output-directory $RESULTPATH \
--output-file-prefix ${ID} \
--enable-duplicate-marking true \
--enable-variant-caller true \
--vc-emit-ref-confidence GVCF \
--enable-cnv true \
--cnv-enable-self-normalization true \
--enable-sv true \
--enable-hla true \
--hla-enable-class-2 true \
--repeat-genotype-enable true \
--qc-detect-contamination true \
--output-format cram

dragen_lic -f genome &> ${RESULTPATH}/dragen_log.end.log


