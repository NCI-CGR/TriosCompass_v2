
--- 
### Prepare GIAB 40X
#### Download Illumina NGS data from GIAB
GIAB Trio WGS fastq data (Illumina; 150bp paired-end reads; ~300X coverage) were downloaded from https://ftp-trace.ncbi.nlm.nih.gov. A detailed manifest file, including md5 checksum, is available at https://42basepairs.com/download/s3/giab/data_indexes/AshkenazimTrio/sequence.index.AJtrio_Illumina300X_wgs_07292015.  

```bash
### Get the manifest file
wget https://raw.githubusercontent.com/genome-in-a-bottle/giab_data_indexes/master/AshkenazimTrio/sequence.index.AJtrio_Illumina300X_wgs_07292015

The above http link is dead but the same file is still available at https://42basepairs.com/download/s3/giab/data_indexes/AshkenazimTrio/sequence.index.AJtrio_Illumina300X_wgs_07292015 .

### Download the fastq files in the manifest file from ftp-trace.ncbi.nlm.nih.gov and validate with md5checksum

```
#### Downsampling to 40X
The GIAB fastq data were trimmed by *fastp* and then aligned to the human reference genome hg38 using *parabricks fq2bam*. Aligned read counts *$X* were retrieved by *samtools idxstat*, the fraction to be downsampled was calculated using the formula: 

  Fraction=3.3*10^9 * 40 /(148 * $X),

where 3.3X10^9 is the size of the human genome, 40 is the target coverage, and 148 is the actual read length after *fastp* trimming.  

*Sambamba* was utilized to downsample the resulting bam files to 40X coverage.  Further details regarding this data processing are provided within the Snakemake workflow available at https://github.com/NCI-CGR/TriosCompass_v2/blob/GIAB_Trios/Snakefile_parabricks. 


---

### Compare DeepTrio with TriosCompass
#### Run DeepTrio (V1.8.0)
+ run_it_module_wg.sh
```bash
#!/bin/bash

module load deepvariant/1.8.0-deeptrio


# Define input and output directories (adjust these as needed)
INPUT_DIR=/data/DCEG_Trios/TriosCompass_Manuscript/GIAB
OUTPUT_DIR="${PWD}/output_GIAB40X_module_WG"
mkdir -p "${OUTPUT_DIR}"

BIN_VERSION=1.8.0
N_SHARDS=47

REF_FASTA="${INPUT_DIR}/ref/Homo_sapiens_assembly38.fasta"
CHILD_BAM="${INPUT_DIR}/bam/HG002_NA24385_son_40X.bam"
PARENT1_BAM="${INPUT_DIR}/bam/HG003_NA24149_father_40X.bam"
PARENT2_BAM="${INPUT_DIR}/bam/HG004_NA24143_mother_40X.bam"

SAMPLE_CHILD="HG002"
SAMPLE_PARENT1="HG003"
SAMPLE_PARENT2="HG004"

# Set environment variables (important for DeepVariant stability)
export OPENBLAS_NUM_THREADS=1
export KEY_CPU_THREADS_NUM=1

MYTMP=TMP_MODULE
mkdir -p  $MYTMP $OUTPUT_DIR
mkdir -p "${OUTPUT_DIR}/intermediate_results_dir"


export TMPDIR=$MYTMP

# Run DeepTrio using Singularity
/usr/bin/time -v run_deeptrio \
    --model_type=WGS \
    --ref=$REF_FASTA \
    --reads_child=$CHILD_BAM \
    --reads_parent1=$PARENT1_BAM \
    --reads_parent2=$PARENT2_BAM \
    --output_vcf_child "${OUTPUT_DIR}/${SAMPLE_CHILD}.output.vcf.gz" \
    --output_vcf_parent1 "${OUTPUT_DIR}/${SAMPLE_PARENT1}.output.vcf.gz" \
    --output_vcf_parent2 "${OUTPUT_DIR}/${SAMPLE_PARENT2}.output.vcf.gz" \
    --sample_name_child "${SAMPLE_CHILD}" \
    --sample_name_parent1 "${SAMPLE_PARENT1}" \
    --sample_name_parent2 "${SAMPLE_PARENT2}" \
    --num_shards "${N_SHARDS}" \
    --intermediate_results_dir "${OUTPUT_DIR}/intermediate_results_dir" \
    --output_gvcf_child "${OUTPUT_DIR}/${SAMPLE_CHILD}.g.vcf.gz" \
    --output_gvcf_parent1 "${OUTPUT_DIR}/${SAMPLE_PARENT1}.g.vcf.gz" \
    --output_gvcf_parent2 "${OUTPUT_DIR}/${SAMPLE_PARENT2}.g.vcf.gz" 
```

```bash
cd /data/DCEG_Trios/TriosCompass_Manuscript/test_deeptrio

sbatch -t 80:00:00 --mem=80g --export=ALL  --cpus-per-task=48 --wrap='./run_it_module_wg.sh '

ls -altr output_GIAB40X_module_WG
total 1697040
drwxr-s---. 2 zhuw10 DCEG_Trios      4096 Apr 24 18:56 intermediate_results_dir
-rw-r-----. 1 zhuw10 DCEG_Trios 117731959 Apr 24 22:01 HG004.output.vcf.gz
-rw-r-----. 1 zhuw10 DCEG_Trios 427252884 Apr 24 22:01 HG004.g.vcf.gz
-rw-r-----. 1 zhuw10 DCEG_Trios   1626084 Apr 24 22:01 HG004.output.vcf.gz.tbi
-rw-r-----. 1 zhuw10 DCEG_Trios 117203520 Apr 24 22:02 HG002.output.vcf.gz
-rw-r-----. 1 zhuw10 DCEG_Trios 468798498 Apr 24 22:02 HG002.g.vcf.gz
-rw-r-----. 1 zhuw10 DCEG_Trios    635438 Apr 24 22:02 HG004.g.vcf.gz.tbi
-rw-r-----. 1 zhuw10 DCEG_Trios   1633425 Apr 24 22:02 HG002.output.vcf.gz.tbi
-rw-r-----. 1 zhuw10 DCEG_Trios    650204 Apr 24 22:02 HG002.g.vcf.gz.tbi
-rw-r-----. 1 zhuw10 DCEG_Trios 116233479 Apr 24 22:02 HG003.output.vcf.gz
-rw-r-----. 1 zhuw10 DCEG_Trios 483713157 Apr 24 22:02 HG003.g.vcf.gz
-rw-r-----. 1 zhuw10 DCEG_Trios   1632118 Apr 24 22:02 HG003.output.vcf.gz.tbi
drwxr-s---. 2 zhuw10 DCEG_Trios      4096 Apr 24 22:02 .
-rw-r-----. 1 zhuw10 DCEG_Trios    654758 Apr 24 22:02 HG003.g.vcf.gz.tbi
drwxr-s---. 2 zhuw10 DCEG_Trios      4096 Apr 24 22:02 ..

dj 54881480
jobid: 54881480
jobidarray: -
jobname: wrap
user: zhuw10
state: COMPLETED
submit_time: 2025-04-24T08:25:31
state_reason: -
state_desc: -
priority: 917349
partition: norm
qos: global
reservation: -
nodes: 1
cpus: 48
gpus: -
mem: 80 GB
timelimit: 10-00:00:00
gres: -
dependency: -
alloc_node: cn0062
command: -
work_dir: /vf/users/DCEG_Trios/TriosCompass_Manuscript/test_deeptrio
std_in: /dev/null
std_out: -
std_err: -
start_time: 2025-04-24T08:29:21
queued_time: 3:50
nodelist: cn0020
end_time: 2025-04-24T22:02:50
elapsed_time: 13:33:29
exit_code: 0
cpu_cur: 1
cpu_min: 1
cpu_max: 103
cpu_avg: 40.76
cpu_util: 84.9%
mem_cur: 258 MB
mem_min: 253 MB
mem_max: 44 GB
mem_avg: 19 GB
mem_util: 23.5%
gpu_cur: -
gpu_min: -
gpu_max: -
gpu_avg: -
gpu_util: -
D_cur: 0
D_min: 0
D_max: 20
D_avg: 1.00
eval: -
comment: -
eligible_time: 2025-04-24T08:25:31
suspended_time: -

### Get sbatch std output from the file slurm-54881480.out
User time (seconds): 1769619.88
        System time (seconds): 8695.86
        Percent of CPU this job got: 3644%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 13:33:20
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 15540728
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 24551
        Minor (reclaiming a frame) page faults: 2952877814
        Voluntary context switches: 154707075
        Involuntary context switches: 40654410
        Swaps: 0
        File system inputs: 4729575802
        File system outputs: 834337544
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
```

####  Merge VCFs using GLnexus
+ Ref: https://github.com/google/deepvariant/blob/r1.8/docs/deeptrio-wgs-case-study.md
```bash
cd /data/DCEG_Trios/TriosCompass_Manuscript/test_deeptrio

module load singularity

BIN_VERSION=1.8.0
OUTPUT_DIR="${PWD}/output_GIAB40X_module_WG"

# by default use TMPDIR=/data/zhuw10/TMP
singularity run \
  -B "${OUTPUT_DIR}":"/output" \
  -B /data:/data \
  docker://quay.io/mlin/glnexus:v1.2.7 \
  /usr/local/bin/glnexus_cli \
  --config DeepVariant_unfiltered \
  /output/HG002.g.vcf.gz \
  /output/HG003.g.vcf.gz \
  /output/HG004.g.vcf.gz \
  | singularity run docker://google/deepvariant:deeptrio-"${BIN_VERSION}" \
    bcftools view - \
  | singularity run docker://google/deepvariant:deeptrio-"${BIN_VERSION}" \
    bgzip -c > ${OUTPUT_DIR}/HG002_trio_merged.vcf.gz

ls -altr output_GIAB40X_module_WG/HG002_trio_merged.vcf.gz
-rw-r----- 1 zhuw10 DCEG_Trios 224470226 Apr 25 10:18 output_GIAB40X_module_WG/HG002_trio_merged.vcf.gz

### normalize
module load bcftools 
bcftools norm -f ../GIAB/ref/Homo_sapiens_assembly38.fasta -m - -O z -o ${OUTPUT_DIR}/HG002_trio_merged.norm.vcf.gz ${OUTPUT_DIR}/HG002_trio_merged.vcf.gz
```

#### Call DNMs 
##### Approach 1. Call DNMs using RTG
+ Ref: https://github.com/google/deepvariant/blob/r1.8/docs/deeptrio-wgs-case-study.md
```bash
OUTPUT_DIR="${PWD}/output_GIAB40X_module_WG"

singularity run \
  -B "../GIAB/ref":"/reference" \
  docker://realtimegenomics/rtg-tools format \
  -o /reference/Homo_sapiens_assembly38.sdf "/reference/Homo_sapiens_assembly38.fasta"

### run rtg mendelian on the normalized vcf data
singularity run \
    --pwd /mnt -B /data --bind $(pwd):/mnt \
    -B "../GIAB/ref":"/reference" \
    -B "${OUTPUT_DIR}":"/output" \
    docker://realtimegenomics/rtg-tools mendelian \
    -i "/output/HG002_trio_merged.norm.vcf.gz" \
    -o "/output/HG002_trio_annotated.norm.output.vcf.gz" \
    --pedigree=GIAB.ped \
    -t /reference/Homo_sapiens_assembly38.sdf \
    | tee ${OUTPUT_DIR}/deepvariant.input_rtg_norm_output.txt


### Call DNM
bcftools filter  -i "MCV='HG002:0/0+0/0->0/1'" ${OUTPUT_DIR}/HG002_trio_annotated.norm.output.vcf.gz | bcftools view -s HG002 -Ov -o HG002.dnm.vcf 

grep -c -v "^#" HG002.dnm.vcf
14053
```

#### Approach 2. Call DNMs using slivar
```bash
slivar --version
# > slivar version: 0.2.8 23df117c3809a2bf57eb0dd1fffdca0df06252a3

# Ref: https://github.com/NCI-CGR/TriosCompass_v2/blob/manuscript/config/GIAB_40X.yaml#L47-L62
slivar expr  --vcf output_GIAB40X_module_WG/HG002_trio_merged.norm.vcf.gz \
             --ped  GIAB.ped     \
             --pass-only \
             --out-vcf HG002_slivar.dnm.vcf \
             --trio "denovo:( \
              kid.GQ >= 15 && kid.het && \
              kid.AB > 0.25 && kid.AB < 1-0.25 && \
              (kid.AD[0]+kid.AD[1]) >= 15 && (kid.AD[0]+kid.AD[1]) < 125 && \
              (mom.GQ >= 20 && mom.hom_ref) && \
              (dad.GQ >= 20 && dad.hom_ref) && \
              mom.AD[1]/(mom.AD[0]+mom.AD[1]) < 0.02 && \
              dad.AD[1]/(dad.AD[0]+dad.AD[1]) < 0.02 && \
              (mom.AD[0]+mom.AD[1]) >= 15 && (mom.AD[0]+mom.AD[1]) < 125 && \
              (dad.AD[0]+dad.AD[1]) >= 15 && (dad.AD[0]+dad.AD[1]) < 125 \
            )"

### Extract the sample HG002 from HG002_slivar.dnm.vcf
bgzip -c HG002_slivar.dnm.vcf > HG002_slivar.dnm.vcf.gz
tabix HG002_slivar.dnm.vcf.gz
bcftools view -s HG002 -R ../GIAB/ref/hg38.wgs_interval.bed HG002_slivar.dnm.vcf.gz -O z -o HG002_silvar_final.dnm.vcf.gz
tabix HG002_silvar_final.dnm.vcf.gz

zgrep -c -v "^#" HG002_silvar_final.dnm.vcf.gz
#1635

# The processing for 80X is similar, except some changes in the slivar expression 
# ref: https://github.com/NCI-CGR/TriosCompass_v2/blob/manuscript/config/GIAB_80X.yaml#L47-L62
slivar expr  --vcf ${OUTPUT_DIR}/HG002_trio_merged.80X_norm.vcf.gz \
             --ped  GIAB.ped     \
             --pass-only \
             --out-vcf HG002_slivar.80X_dnm.vcf \
             --trio "denovo:( \
              kid.GQ >= 13 && kid.het && \
              kid.AB > 0.3 && kid.AB < 1-0.3 && \
              (kid.AD[0]+kid.AD[1]) >= 20 && (kid.AD[0]+kid.AD[1]) < 250 && \
              (mom.GQ >= 20 && mom.hom_ref) && \
              (dad.GQ >= 20 && dad.hom_ref) && \
              mom.AD[1]/(mom.AD[0]+mom.AD[1]) < 0.02 && \
              dad.AD[1]/(dad.AD[0]+dad.AD[1]) < 0.02 && \
              (mom.AD[0]+mom.AD[1]) >= 20 && (mom.AD[0]+mom.AD[1]) < 250 && \
              (dad.AD[0]+dad.AD[1]) >= 20 && (dad.AD[0]+dad.AD[1]) < 250 \
            )"
# sample  denovo
# HG002   1653
```

### Benchmark using hap.py
#### Prepare the truth file
+ ref: [DenovoCNN_Supplementary.pdf](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/50/17/10.1093_nar_gkac511/1/gkac511_supplemental_files.zip?Expires=1750873210&Signature=Ez5ByyNzjgOyIsZnPihvvRXaJbdhKykZFGkcxFqZjlmfGGVUcO3vhIrfMcCb1cFbQ2pkS6mE7diZadVH~Hri-ltzIXtbbO-bl~cusfDeR5bI51yPB1qw4t4g8qfBBBTnnEF9hB39JGJEjrMoOi1MOBFTwORMiCV7S~dbG0eHQFea0UXr46vgOIaoE7V3wMwZPQ6nPkUg1RvlqnVcejoy652vh437FvASI9qhbLsk7-VX99hF3S1ICEYPPpCuFGXrkmlETWCpHzAqLsk4fZkqk0GHzgCfrpenNizwCLxMVeIo0RrN29eE3ssidBoO9Zv218sgSaelzm-mHuhuOSUtzA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

```bash
cd /vf/users/DCEG_Chernobyl/NP0436-HE5/GIAB_DATA
mkdir truth_files

# get highconf regions of each subject
wget --directory-prefix=truth_files/ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/supplementaryFiles/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_trioinconsistent_nodenovo_slop50.bed

wget --directory-prefix=truth_files/ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh38/supplementaryFiles/HG003_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.bed

wget --directory-prefix=truth_files/ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh38/supplementaryFiles/HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.bed

# get the intersection of highconf regions of all trio members
# => AJtrio_3.3.2_hc.bed 
module load bedtools
# [+] Loading bedtools  2.31.1

cd /data/DCEG_Chernobyl/NP0436-HE5/GIAB_DATA/truth_files

multiIntersectBed -i HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf.bed HG003_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.bed HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.bed -names A B C | grep "A,B,C" > AJtrio_3.3.2_hc.bed 


### Get DNM candidates and nodenovo_slop50
wget --directory-prefix=truth_files/ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/supplementaryFiles/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_trioinconsistent_nodenovo_slop50.bed
wget --directory-prefix=truth_files/ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/supplementaryFiles/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_trioinconsistent.vcf.gz.tbi
wget --directory-prefix=truth_files/ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/supplementaryFiles/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_trioinconsistent.vcf.gz

# The purpose of nondenovo_slop50 is to exclude those 1071 cases of trioinconsistent other than "0/0+0/0->0/1"
wc -l ../../truth_files/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_trioinconsistent_nodenovo_slop50.bed
1071 ../../truth_files/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_trioinconsistent_nodenovo_slop50.bed

 zgrep -v -e "0/0+0/0->0/1" -e "^#"  ../../truth_files/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_trioinconsistent.vcf.gz |wc -l
1071

# substract nodenovo_slop50 from high conf regions 
# => AJtrio_3.3.2.truth.bed
bedtools subtract -a ../../truth_files/AJtrio_3.3.2_hc.bed -b ../../truth_files/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_trioinconsistent_nodenovo_slop50.bed > AJtrio_3.3.2.truth.bed

### AJtrio_truth.vcf.gz
13M.call_dnm.md:bcftools view -s child  -R ../../truth_files/AJtrio_3.3.2_hc.bed  ../../truth_files/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_trioinconsistent.vcf.gz -o AJtrio_truth.vcf.gz -O z

zgrep -c -v "^#" AJtrio_truth.vcf.gz
2494

### so the total number of the real DNMs is: 2494-1071=1423 
```


#### Test the DNM calling Approach 1 

```bash
module load bcftools rtg-tools hap.py
# [+] Loading bedtools  2.31.1 
# [+] Loading rtg-tools  3.12.1 
# [+] Loading hap.py  0.3.14

hap.py $GIAB/giab_output/AJtrio_truth.vcf.gz \
   HG002.dnm.vcf  \
    -f  $GIAB/giab_output/AJtrio_3.3.2.truth.bed \
    -r $GIAB/giab_output/ref/Homo_sapiens_assembly38.fasta \
    -o happy_out/test1 \
    --engine=vcfeval --threads 4  --fixchr --bcftools-norm
```

#### Test the DNM calling Approach 
```bash
hap.py $GIAB/giab_output/AJtrio_truth.vcf.gz \
   HG002_silvar_final.dnm.vcf.gz \
    -f  $GIAB/giab_output/AJtrio_3.3.2.truth.bed \
    -r $GIAB/giab_output/ref/Homo_sapiens_assembly38.fasta \
    -o happy_out/HG002_silvar_final \
    --engine=vcfeval --threads 4  --fixchr --bcftools-norm
```

#### Benchmark the results from TriosCompass
```bash
hap.py $GIAB/giab_output/AJtrio_truth.vcf.gz \
   ../GIAB/output_40X/dnm_vcf/GIAB.dnm.vcf.gz  \
    -f  $GIAB/giab_output/AJtrio_3.3.2.truth.bed \
    -r $GIAB/giab_output/ref/Homo_sapiens_assembly38.fasta \
    -o happy_out/trioscompass \
    --engine=vcfeval --threads 4  --fixchr --bcftools-norm
```
