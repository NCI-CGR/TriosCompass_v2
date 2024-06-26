<!-- omit in toc -->
# TriosCompass_v2
This is a trios analysis workflow written in Snakemake.

---
- [Introduction](#introduction)
- [Installation](#installation)
  - [Locally installed bioinformatics tools](#locally-installed-bioinformatics-tools)
    - [GangSTR](#gangstr)
    - [WhatsHap](#whatshap)
    - [HipSTR and revision (under TriosCompass\_v2 clone folder)](#hipstr-and-revision-under-trioscompass_v2-clone-folder)
  - [Preparation of the calling regions of DNMs](#preparation-of-the-calling-regions-of-dnms)
  - [The reference genome hg38](#the-reference-genome-hg38)
  - [The reference panels for GangSTR and HipSTR](#the-reference-panels-for-gangstr-and-hipstr)
- [Some details about the DNM calling](#some-details-about-the-dnm-calling)
  - [*Slivar* expression to select DNMs](#slivar-expression-to-select-dnms)
- [Development and benchmark of TriosCompass](#development-and-benchmark-of-trioscompass)
  - [Benchmark on the GIAB trios](#benchmark-on-the-giab-trios)
  - [Development and benchmark on the 8 trios from the Chernobyl data set](#development-and-benchmark-on-the-8-trios-from-the-chernobyl-data-set)
- [Applications of TriosCompass on the real data](#applications-of-trioscompass-on-the-real-data)


---

## Introduction

This Snakemake workflow takes fastq/bam files and pedigree files as input, and generate lists of de novo mutations (DNMs) for each trios in the end. We had employed NVIDIA Parabricks toolbox in the pipeline, to accelerate NGS analysis with GPUs.

Briefly, there are 3 major components in the workflow:
+ Call DNMs
  1. Variant calling by GATK HaplotypeCaller, DeepVariant (and Strelka).
  2. Identify DNMs by [Slivar](https://github.com/brentp/slivar).
  3. Ensemble call. 
  4. DNM candidates visualized by JIGV 
+ Extract parental origin of DNMs
+ Call mDNMs (microsatellites DNMs; or dnSTR for de novo short tandem repeat).


---

## Installation

The workflow has been tested under biowulf with snakemake version 7.3.7.  Snakemake is installed under conda, using a command like blow:
```bash
conda  create -c conda-forge -c bioconda -n snakemake snakemake=7.3.7
```

There are other dependencies required by the workflow, most of which are available as biowulf modules:
+ Biowulf modules:
  + bcftools/1.19
  + bedtools/2.31.1
  + fastp/0.23.2
  + GATK/4.3.0.0
  + glnexus/1.4.1
  + hipstr/0.7
  + igv/2.12.3 (optional)
  + parabricks/4.0.0
  + perl/5.38
  + picard/2.27.3
  + samtools/1.19
  + singularity/4.0.1
+ Conda packages (locally installed)
  + csvtk (V0.28.0)
  + dumpSTR (V5.0.2)
  + [slivar](https://github.com/brentp/slivar) (V0.2.8 23df117c3809a2bf57eb0dd1fffdca0df06252a3)
+ Locally installed
  + GangSTR (V2.5.0)
  + WhatsHap (V2.1)
  + HipSTR (with revision)
+ Containers
  + docker://gymreklab/monstr (sha256:f5e90770d3a9dd1220764986c11f0cd8d345ed03955e36374bcbcdabfe8cb71d)

Besides, some resources files are also needed:
+ The reference genome
  + Homo_sapiens_assembly38.fasta
+ Calling regions of DNMs
+ The reference panels for GangSTR and HipSTR
+ GRCh38GenomicSuperDup.bed.gz
  + Obtained from the UCSC Table Browser (hg38.genomicSuperDups table)

### Locally installed bioinformatics tools
#### GangSTR
```bash
wget https://github.com/gymreklab/GangSTR/releases/download/v2.5/GangSTR-2.5.0-Source.tar.gz

tar xvfz GangSTR-2.5.0-Source.tar.gz
cd GangSTR-2.5.0-Source

mkdir build
cd build
cmake ..
make 
cmake --install . --prefix $HOME
-- Install configuration: ""
-- Installing: /home/zhuw10/bin/GangSTR
-- Set runtime path of "/home/zhuw10/bin/GangSTR" to ""

which GangSTR
~/bin/GangSTR
```

#### WhatsHap
```bash
pip3 install --user whatshap
which whatshap
~/.local/bin/whatshap
```

#### HipSTR and revision (under TriosCompass_v2 clone folder)
By default, VizAln generates output html page as temporary file with random name.  We made revision so that the HTML output is written to the specific file using VizAln_rev. 
```bash
cd TrosCompass_v2
git clone https://github.com/tfwillems/HipSTR.git

cp scripts/VizAln_rev HipSTR/


```

+ diff HipSTR/VizAln_rev HipSTR/VizAln
```diff
3c3
< # Usage: VizAln alns.html.gz output_html chrom start [sample]
---
> # Usage: VizAln alns.html.gz chrom start [sample]
8c8
< if [ "$#" -eq 4 ] || [ "$#" -eq 5 ]; then
---
> if [ "$#" -eq 3 ] || [ "$#" -eq 4 ]; then
56,58c56,57
<     f=$2
<     chrom=$3
<     start=$4
---
>     chrom=$2
>     start=$3
60,61c59,60
<     # f="/tmp/vizaln.$RANDOM.html"
<     if [ "$#" -eq 4 ]
---
>     f="/tmp/vizaln.$RANDOM.html"
>     if [ "$#" -eq 3 ]
65c64
<       tabix $1 $chrom:$start-$stop -h | awk -v START=$start '$2 == START || $2 == "#"' | awk -v SAMPLE=$5 '$4 == SAMPLE || $4 == "ALL"' | cut -f 1-4 --complement | python ${parent_path}/scripts/generate_aln_html.py > $f
---
>       tabix $1 $chrom:$start-$stop -h | awk -v START=$start '$2 == START || $2 == "#"' | awk -v SAMPLE=$4 '$4 == SAMPLE || $4 == "ALL"' | cut -f 1-4 --complement | python ${parent_path}/scripts/generate_aln_html.py > $f
68c67
<     # open $f 2>/dev/null || { echo >&2 "VizAln successfully generated the alignment visualization file, but failed to open it. Please open $f using a web browser";}
---
>     open $f 2>/dev/null || { echo >&2 "VizAln successfully generated the alignment visualization file, but failed to open it. Please open $f using a web browser";}
71c70
< Usage: VizAln alns.html.gz output_html chrom start [sample]
---
> Usage: VizAln alns.html.gz chrom start [sample]
```

---

### Preparation of the calling regions of DNMs
The interval file is to include the selected regions only to call DNMs (via bcftools).  We used the file hg38.wgs_interval.bed in this study, which was converted from *resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list*. The latter is part of [resource bundle hosted by the Broad Institute](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle). There are severals way to retrieve the interval file, for example, ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/wgs_calling_regions.hg38.interval_list.

```bash
### Download calling region of hg38 from the Broad Institute 
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list -O ref/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list

### Then, convert the interval list to bed file using Picard
module load picard 

java -jar $PICARDJAR IntervalListToBed -I ref/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list -O ref/hg38.wgs_interval.bed

bgzip -c ref/hg38.wgs_interval.bed > ref/hg38.wgs_interval.bed.gz
tabix ref/hg38.wgs_interval.bed.gz


```

### The reference genome hg38
Hg38 reference genome is used here, and it was indexed by bwa, samtools (.fai) and Picard CreateSequenceDictionary (.dict).
```bash
-rw-r----- 1 zhuw10 DCEG_Trios     608628 Jul  3 08:39 Homo_sapiens_assembly38.dict
-rw-r----- 1 zhuw10 DCEG_Trios 3249912778 Jul  3 08:39 Homo_sapiens_assembly38.fasta
-rw-r----- 1 zhuw10 DCEG_Trios     487553 Jul  3 08:39 Homo_sapiens_assembly38.fasta.64.alt
-rw-r----- 1 zhuw10 DCEG_Trios      20199 Jul  3 08:39 Homo_sapiens_assembly38.fasta.64.amb
-rw-r----- 1 zhuw10 DCEG_Trios     455474 Jul  3 08:39 Homo_sapiens_assembly38.fasta.64.ann
-rw-r----- 1 zhuw10 DCEG_Trios 3217347004 Jul  3 08:39 Homo_sapiens_assembly38.fasta.64.bwt
-rw-r----- 1 zhuw10 DCEG_Trios  804336731 Jul  3 08:39 Homo_sapiens_assembly38.fasta.64.pac
-rw-r----- 1 zhuw10 DCEG_Trios 1608673512 Jul  3 08:39 Homo_sapiens_assembly38.fasta.64.sa
-rw-r----- 1 zhuw10 DCEG_Trios     160928 Jul  3 08:39 Homo_sapiens_assembly38.fasta.fai
```

### The reference panels for GangSTR and HipSTR
Both GangSTR and HipSTR use STRs specified in one reference panel, but with the different format.  One paper suggested that there is good consistence between HpSTR and GangSTR.
>Oketch, J. W., Wain, L. V, Hollox, E. J., & Hollox, E. (2022). A comparison of software for analysis of rare and common short tandem repeat (STR) variation using human genome sequences from clinical and population-based samples. 1–22. https://doi.org/10.1101/2022.05.25.493473


We would like to have the STR genotypes predicted jointly by both of the two program, therefore we need create one common reference panel file for both of GangSTR and HipSTR. 

```bash
### Get the GangSTR reference panel
wget https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver13.bed.gz -O STR/hg38_ver13.bed.gz

gzip -d STR/hg38_ver13.bed.gz

### reformat for HipSTR
awk -v OFS='\t' '{print $1,$2,$3,$4,($3-$2+1)/$4,"GangSTR_STR_"NR,$5}' STR/hg38_ver13.bed  >  STR/hg38_ver13.hipstr.bed

### HipSTR cannot take STR with unit length > 9 bp
awk -v OFS='\t' '{if($4<=9) print $0}' STR/hg38_ver13.hipstr.bed > STR/hg38_ver13.hipstr_9.bed

### We put the same length restrict to hg38_ver13.bed.gz, so that both of the two reference panels are consistent.
awk -v OFS='\t' '{if($4<=9) print $0}' STR/hg38_ver13.bed > STR/hg38_ver13.le9.bed
```

---

## Some details about the DNM calling
### [*Slivar*](https://github.com/brentp/slivar) expression to select DNMs
```bash
 "denovo:( \
      ( \
          (variant.CHROM == 'chrX' && kid.sex=='male') && \
          kid.hom_alt && kid.AB > 0.98  \
      ) || \
      ( \
          (!(variant.CHROM == 'chrX' && kid.sex=='male')) && \
          kid.het && kid.AB > 0.25 && kid.AB < 0.75 \
      ) \
      ) &&  (kid.AD[0]+kid.AD[1]) >= {params.min_dp}/(1+(variant.CHROM == 'chrX' && kid.sex == 'male' ? 1 : 0)) && \
      mom.hom_ref && dad.hom_ref \
          && (mom.AD[1] + dad.AD[1]) <= 5 \
          && kid.GQ >= {params.min_gq} && mom.GQ >= {params.min_gq} && dad.GQ >= {params.min_gq} \
          && (mom.AD[0]+mom.AD[1]) >= {params.min_dp} && (dad.AD[0]+dad.AD[1]) >= {params.min_dp}/(1+(variant.CHROM == 'chrX' ? 1 : 0))"
```

+ (!(variant.CHROM == 'chrX' && kid.sex=='male')) && kid.het && kid.AB > 0.25 && kid.AB < 0.75
  + In the normal cases, select variants with genotype of "0/1" and allele balance in the range of 0.25 and 0.75.
+ (variant.CHROM == 'chrX' && kid.sex=='male') && alt && kid.AB > 0.98
  + In the special case (that is, variant in chrX in the male offspring), select variants with the genotype of "1/1" and allele balance > 0.98.
+ (kid.AD[0]+kid.AD[1]) >= {params.min_dp}/(1+(variant.CHROM == 'chrX' && kid.sex == 'male' ? 1 : 0))
  + The kid's variant should have depth over the minimal read depth *{params.min_dp}*. If the variant is on chrX and the kid's gender is male, the minimum read depth requirement is reduced to its half.
+ (mom.AD[0]+mom.AD[1]) >= {params.min_dp} && (dad.AD[0]+dad.AD[1]) >= {params.min_dp}/(1+(variant.CHROM == 'chrX' ? 1 : 0))
  + Similar read depth requirement is also applied for the parents.
+ mom.hom_ref && dad.hom_ref
  + The genotypes of the parents should both be "0/0";
+ (mom.AD[1] + dad.AD[1]) <= 5
  + The total read support of the alternative allele (that is, the DNM) should be not bigger than 5 in parents.
+ kid.GQ >= {params.min_gq} && mom.GQ >= {params.min_gq} && dad.GQ >= {params.min_gq}
  + GQ score of the variants should be no less than {params.min_gq} in all the members of the trio.

:bookmark: In the workflow, different settings of {params.min_dp} and {params.min_gq} are assigned for the callers, based on our benchmarking:
+ min_gq=10, min_dp=20 for DeepVariant;
+ min_gq=20, min_dp=30 for GATK and Strelka.


---

## Development and benchmark of TriosCompass
### [Benchmark on the GIAB trios](https://github.com/NCI-CGR/TriosCompass_v2/tree/GIAB_Trios)

### [Development and benchmark on the 8 trios from the Chernobyl data set](https://github.com/NCI-CGR/TriosCompass_v2/tree/8trios)

We selected 8 trios from the 340 Chernobyl samples, which had been published in Yeager et al; Science 2021. This small data set has been used to develop new features and benchmarked with the original DNM predicts which has been manually curated.

---

## Applications of TriosCompass on the real data
+ [Processing 107 new CGR samples](https://github.com/NCI-CGR/TriosCompass_v2/blob/main/Process_107_new_Chernobyl_data.md)
+ [Re-processing 340 old Chernobyl samples](https://github.com/NCI-CGR/TriosCompass_v2/tree/old_chernobyl?tab=readme-ov-file)


