<!-- omit in toc -->
# TriosCompass_v2
This is a trios analysis workflow written in Snakemake.

---
- [Introduction](#introduction)
- [Installation](#installation)
  - [Preparation of the calling regions of DNMs](#preparation-of-the-calling-regions-of-dnms)
  - [The reference genome hg38](#the-reference-genome-hg38)
- [Some details about the DNM calling](#some-details-about-the-dnm-calling)
  - [Slivar expression to select DNMs](#slivar-expression-to-select-dnms)
- [New improvement for Strelka DNM calls](#new-improvement-for-strelka-dnm-calls)
- [Development and benchmark of TriosCompass](#development-and-benchmark-of-trioscompass)
  - [Benchmark on the GIAB trios](#benchmark-on-the-giab-trios)
  - [Development and benchmark on the 8 trios from the Chernobyl data set](#development-and-benchmark-on-the-8-trios-from-the-chernobyl-data-set)
- [Applications of TriosCompass on the real data](#applications-of-trioscompass-on-the-real-data)


---

## Introduction

This Snakemake workflow takes fastq/bam files and pedigree files as inupt, and generate lists of de novo mutations (DNMs) for each tios in the end. We had employed NVIDIA Parabricks toolbox in the pipeline, so as to accelerate NGS analysis with GPUs.

Briefly, there are 3 major components in the workflow:
+ Call DNMs
  1. Variant calling by GATK HaplotypeCaller, DeepVariant (and Strelka).
  2. Identify DNMs by Slivar.
  3. Ensemble call. 
  4. DNM candidates visulized by JIGV 
+ Extract parental origin of DNMs
+ Call mDNMs (microsatillites DNMs; or dnSTR for de novo short tandem repeat).


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
  + slivar (V0.2.8 23df117c3809a2bf57eb0dd1fffdca0df06252a3)
+ Locally installed
  + GangSTR (V2.5.0)
  + WhatsHap (V2.1)
+ Containers
  + docker://gymreklab/monstr (sha256:f5e90770d3a9dd1220764986c11f0cd8d345ed03955e36374bcbcdabfe8cb71d)

Besides, some resources files are also needed:
+ The reference genome
  + Homo_sapiens_assembly38.fasta
  + Specify calling regions of DNMs
    

### Preparation of the calling regions of DNMs
The interval file is to include the selected regions only to call DNMs (via bcftools).  We used the file hg38.wgs_interval.bed in this study, which was converted from *resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list*. The latter is part of [resource bundle hosted by the Broad Institute](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle). There are severals way to retrieve the interval file, for example, ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/wgs_calling_regions.hg38.interval_list.

```bash
### Download calling region of hg38 from the Broad Institute 
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list -O ref/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list

### Then, convert the interval list to bed file uisng Picard
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

---

## Some details about the DNM calling
### Slivar expression to select DNMs
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
  + The genotypeos of the parents should both be "0/0";
+ (mom.AD[1] + dad.AD[1]) <= 5
  + The total read support of the alternative allele (that is, the DNM) should be not bigger than 5 in parents.
+ kid.GQ >= {params.min_gq} && mom.GQ >= {params.min_gq} && dad.GQ >= {params.min_gq}
  + GQ score of the variants should be no less than {params.min_gq} in all the members of the trio.

:bookmark: In the workflow, different settings of {params.min_dp} and {params.min_gq} are assigned for the callers, based on our benchmarking:
+ min_gq=10, min_dp=20 for DeepVariant;
+ min_gq=20, min_dp=30 for GATK and Strelka.

---

## New improvement for Strelka DNM calls
We had explored to use Strelka to call DNMs, together with DeepVariant and HaplotypeCaller, as demonstrated in the file [Snakefile](./Snakefile).
+ [New changes in the Snakefile](https://github.com/NCI-CGR/TriosCompass_v2/commit/36a6726943af9f0473edfc3150e1121bc2c994ee#diff-47959dfd378b3cd1d39b5515418ee8e4444ab7a6036d5197a5bea82814f928a3):
  + Restrict to “pass” only.
  + Use percentage in filters to address those variants with lower depth.
  + Add new IGV snapshots.
    + New perl script [prepare_yml_from_vcf.pl](./scripts/prepare_yml_from_vcf.pl) to prepare YAML file required by igv_snapshot_maker.
    + Install igv-snapshot-maker to run the updated workflow.
      + pip install -i https://test.pypi.org/simple/ igv-snapshot-maker==1.0.0
+ Updated [summary Excel file](https://github.com/NCI-CGR/TriosCompass_v2/blob/main/cgr_summary.xlsx).
+ New output files under output/call_ism/strelka
  + Users may copy the files (output/call_ism/strelka) to local and start with the Excel files, for example, like output/call_ism/strelka/t0007c1.xlsx
  + The updated DNM calls from Strelka and JIGV should remain in the same location as before under /data/DCEG_Trios/new_cgr_data/TriosCompass_v2.
```bash
### The overall fold struture
tree -L 4  output/call_ism/strelka/t0007c1.xlsx
output/call_ism/
└── strelka
    ├── IGV_Snapshots
    │   ├── t0007c1
    │   │   ├── t0007c1.bat
    │   │   ├── t0007c1_MIE_00001.bat
    │   │   ├── t0007c1_MIE_00001.png
    │   │   ├── t0007c1_MIE_00002.bat
    │   │   ├── t0007c1_MIE_00002.png
    │   │   ├── t0007c1_MIE_00003.bat
    │   │   ├── t0007c1_MIE_00003.png
    │   │   ├── t0007c1_MIE_00004.bat
    │   │   ├── t0007c1_MIE_00004.png
    │   │   ├── t0007c1_MIE_00005.bat
    │   │   ├── t0007c1_MIE_00005.png
    │   │   ├── t0007c1_MIE_00006.bat
    │   │   ├── t0007c1_MIE_00006.png
    │   │   └── t0007c1_ROIs.bat
    │   ├── t0022c1
    │   │   ├── t0022c1.bat
    │   │   ├── t0022c1_MIE_00001.bat
    │   │   ├── t0022c1_MIE_00001.png
    │   │   ├── t0022c1_MIE_00002.bat
    │   │   ├── t0022c1_MIE_00002.png
    │   │   ├── t0022c1_MIE_00003.bat
    │   │   ├── t0022c1_MIE_00003.png
    │   │   ├── t0022c1_MIE_00004.bat
    │   │   ├── t0022c1_MIE_00004.png
    │   │   ├── t0022c1_MIE_00005.bat
    │   │   ├── t0022c1_MIE_00005.png
    │   │   └── t0022c1_ROIs.bat
    │   ├── t0042c1
    │   │   ├── t0042c1.bat
    │   │   ├── t0042c1_MIE_00001.bat
    │   │   ├── t0042c1_MIE_00001.png
    │   │   ├── t0042c1_MIE_00002.bat
    │   │   ├── t0042c1_MIE_00002.png
    │   │   ├── t0042c1_MIE_00003.bat
    │   │   ├── t0042c1_MIE_00003.png
    │   │   ├── t0042c1_MIE_00004.bat
    │   │   ├── t0042c1_MIE_00004.png
    │   │   ├── t0042c1_MIE_00005.bat
    │   │   ├── t0042c1_MIE_00005.png
    │   │   ├── t0042c1_MIE_00006.bat
    │   │   ├── t0042c1_MIE_00006.png
    │   │   ├── t0042c1_MIE_00007.bat
    │   │   ├── t0042c1_MIE_00007.png
    │   │   ├── t0042c1_MIE_00008.bat
    │   │   ├── t0042c1_MIE_00008.png
    │   │   ├── t0042c1_MIE_00009.bat
    │   │   ├── t0042c1_MIE_00009.png
    │   │   ├── t0042c1_MIE_00010.bat
    │   │   ├── t0042c1_MIE_00010.png
    │   │   ├── t0042c1_MIE_00011.bat
    │   │   ├── t0042c1_MIE_00011.png
    │   │   ├── t0042c1_MIE_00012.bat
    │   │   ├── t0042c1_MIE_00012.png
    │   │   └── t0042c1_ROIs.bat
    │   ├── t0058c1
    │   │   ├── t0058c1.bat
...
    ├── t0007c1.xlsx
    ├── t0007c1.yaml
    ├── t0022c1.xlsx
...
    └── t0766c2.yaml

41 directories, 938 files


### Example of the command to transfer data to local (run it at your laptop)
scp -r helix:/data/DCEG_Trios/new_cgr_data/TriosCompass_v2/output/call_ism/strelka ~/Downloads/strelka

### Then open ~/Downloads/strelka/t0007c1.xlsx to review the snapshots
```
---

## Development and benchmark of TriosCompass
### [Benchmark on the GIAB trios](https://github.com/NCI-CGR/TriosCompass_v2/tree/GIAB_Trios)

### [Development and benchmark on the 8 trios from the Chernobyl data set](https://github.com/NCI-CGR/TriosCompass_v2/tree/8trios)

We selected 8 trios from the 340 Chernobyl samples, which had been published in Yeager et al; Science 2021. This small data set has been used to develop new features and benchmarked with the original DNM predictes which has been manually curated.

---

## Applications of TriosCompass on the real data
+ Processing 107 new CGR samples
+ Re-processing 340 old Chernobyl samples
