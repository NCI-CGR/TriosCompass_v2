<!-- omit in TOC -->
# Snakemake workflow: TriosCompass

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.3.7-brightgreen.svg)](https://snakemake.github.io)


A Snakemake workflow for DNM (de novo mutation) calling.

---
- [Snakemake workflow: TriosCompass](#snakemake-workflow-trioscompass)
  - [Overview](#overview)
    - [I. Introduction](#i-introduction)
    - [II. Dependencies](#ii-dependencies)
    - [III. Methods](#iii-methods)
      - [A. Call DNMs](#a-call-dnms)
      - [B. Phase DNMs](#b-phase-dnms)
      - [C. Call dnSTRs](#c-call-dnstrs)
      - [D. Call dnSVs](#d-call-dnsvs)
  - [User's guides](#users-guides)
    - [I. Installation](#i-installation)
    - [II. Inputs](#ii-inputs)
      - [Reference genome](#reference-genome)
      - [The STR reference panel to call dnSTRs](#the-str-reference-panel-to-call-dnstrs)
    - [III. Outputs](#iii-outputs)
    - [IV. Run TrisCompass](#iv-run-triscompass)


---

## Overview

### I. Introduction
TriosCompass consists of 4 functional components:
+ Call DNMs (de novo mutations) using *DeepVariant* and *GATK HaplotypeCaller*
+ Phase DNMs using *whatshap*
+ Call dnSTR (de novo simple tandem repeats) using HipSTR and MonSTR.
+ Call dnSV (de novo structural variants) using *Manta*, *GraphType2* and *smoove*

---

### II. Dependencies

All required bioinformatics tools are wrapped as conda, container and etc, so there is no need to for users to install any of them.

Nevertheless, there are still some dependencies required to start Snakemake workflow, which have been specified in [environment.yaml](environment.yaml).  Users can create a new *conda* env for TriosCompassV2
```bash

mamba env create -f environment.yaml

conda activate TriosCompassV2
```

Besides, singularity needs to be installed globally. 

### III. Methods
#### A. Call DNMs

The DNM candidates are jointly called by both *DeepVariant (DV)* and *GATK HaplotypeCaller (HC)* in the specified callable regions, then filtered using [*sliver*](https://github.com/brentp/slivar):

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

:bookmark: In the workflow, the parameters {params.min_dp} and {params.min_gq} can be configured via config/config.yaml, so are the callable regions:

+ config/config.yaml
```yml
call_dnm:
  interval: "ref/hg38.wgs_interval.bed" # DNM callable regions
  dv:
    min_gq: 3
    min_dp: 20
  hc:
    min_gq: 20
    min_dp: 30
```

---
#### B. Phase DNMs
We identified parental origin of DNMs using [*WhatsHap*](https://github.com/whatshap/whatshap). 

There are two running modes in WhatsHap: individual and pedigree.  Pedigree mode is ideal to identify parental origin of variants in the child.  However, WhatsHap does not phase DNMs in the pedigree mode, as DNMs do not follow mendelian inheritance. Phasing DNMs is not supported by WhatsHap at present (see this [issue](https://github.com/whatshap/whatshap/issues/82) for the details). However, DNMs are phased in the individual mode.  Therefore, we managed to run WhatsHap in both individual mode and pedigree mode and extract parental original from the two outputs.

We first explored to run WhatsHap on the whole genome using the workflow [Snakefile_whatshap_WG](./Snakefile_whatshap_WG).  It turns out to be very slow: it took up to 7 days to process one trio family, and 2-3 days for each trio in average.  Finally,  we chose to phase a 10 Kb window around each DNM. 

![](img/Trios_workflow_rulegraph_whatshap_100K.png)

We developed [a Perl script](./scripts/extract_parental_origin.pl) to identify the parental origin, and the aggregated results are output to the file: output/phase_DNMs/{Trios_ID}.parental_origin.tab.  The output file is a headless tab delimited text file, with 9 columns including variant id of DNM, parental origin prediction.  

The perl script extracts the haplotype block containing the DNMs from the output of WhatsHap phasing child, and collects the phased variants of child in the block.  It also repeats the processing on the output of WhatsHap phasing trios.  The phase prediction from the latter is given as *F|M*, i.e, the first allele is the one inherited from the father (F) and the other in inherited from the mother (M). However, the phase prediction from the former is not certain, that is, it could be F|M or M|F. Ideally, the phases of the variants in each haplotype block from WhatsHap phasing child only should be either all identical (FM count) or all opposite (MF count) with those from WhatsHap phasing trios.  We simply assign parental origin as Not Determined (ND) if there is inconsistency (i.e., both FM count > 0 and MF count > 0). 

+ Example of the output of extract_parental_origin.pl
  
| VariantID                        | Parental Origin | Phase | Haplotype Block Size | # Informative Sites | FM count | MF count | FM status | Parent origin change |
|----------------------------------|-----------------|-------|----------------------|---------------------|----------|----------|-----------|----------------------|
| chr1:208903875:G:A               | ND              | 1\|0  | 36                   | 0                   | 25       | 11       | ND        | ND=>ND               |
| chr2:12194398:A:C                | paternal        | 1\|0  | 2                    | 2                   | 2        | 0        | FM        | paternal=>paternal   |
| chr2:14986429:A:T                | paternal        | 1\|0  | 5                    | 5                   | 5        | 0        | FM        | paternal=>paternal   |
| chr2:50677044:C:T                | ND              | 0/1   | 0                    | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr2:139840973:T:A               | paternal        | 0\|1  | 14                   | 14                  | 0        | 14       | MF        | paternal=>paternal   |
| chr2:154745781:T:C               | ND              | 0/1   | 0                    | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr2:203705880:T:TTCTTTC         | ND              | 0/1   | 0                    | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr2:220449146:A:G               | ND              | 0/1   | 0                    | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr3:14911256:G:C                | ND              | 0/1   | 0                    | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr3:76410708:A:C                | ND              | 0/1   | 0                    | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr3:185020246:C:T               | ND              | 0/1   | 0                    | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr3:194728646:G:GTGTGTGTGTGTGTC | maternal        | 0\|1  | 6                    | 6                   | 6        | 0        | FM        | maternal=>maternal   |
| chr3:197879976:C:T               | paternal        | 0\|1  | 3                    | 2                   | 0        | 1        | MF        | paternal=>paternal   |
| chr3:197879978:C:A               | paternal        | 0\|1  | 3                    | 2                   | 0        | 1        | MF        | paternal=>paternal   |
| chr4:29657055:GC:G               | ND              | 0/1   | 0                    | 0                   | 0        | 0        | ND        | ND=>ND               |

:bookmark: For most users, the first two columns provide essential information about parental origins of the predicted DNMs.

---

+ Descriptions for each output column.

| Column position | Column Header          | Description                                                   | Notes                                                                                                                |
|-----------------|------------------------|---------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------|
| 1               | VariantID              | Variant ID of DNM                                             |                                                                                                                      |
| 2               | Parental Origin        | Predicted parental origin                                     | ND (not determined)                                                                                                  |
| 3               | Phase                  | Phased genotype of DNMs (from phasing child)                  | Parental origin cannot be determined from the output from phasing child only.                                        |
| 4               | Haplotype Block Size   | Size of the haplotype block (from phasing child)              | Prediction from a large block is more reliable.                                                                      |
| 5               | # Informative Sites    | # informative sites in the haplotype block                    | Prediction is not reliable if # informative sites is 0.                                                              |
| 6               | FM count               | # phased variants supporting F\|M genotype                    |                                                                                                                      |
| 7               | MF count               | # phased variants supporting M\|F genotype                    |                                                                                                                      |
| 8               | FM status              | FM (for F\|M), MF (for M\|F), ND (not determined)             | FM status is ND if both FM and MF counts are bigger than 0.                                                          |
| 9               | Parental origin change | Prediction (wo phasing trios) => Prediction(wi phasing trios) | Generally, the prediction with phasing child only is consistent with that with both phasing child and phasing trios. |



+ config.yaml setting for phasing
```yml
phasing:
  window_size: 10000
  # the perl script path relative to the working directory
  perl_cmd: "perl TriosCompass_v2/workflow/scripts/extract_parental_origin.pl "

```

---

#### C. Call dnSTRs

---

#### D. Call dnSVs

---

## User's guides

### I. Installation 
+ Install singularity globally.
+ Install conda and mamba.
+ Create the workspace directory \$WORKSPACE and change directory to \$WORKSPACE.
+ Git clone the TriosCompass_v2 repository.
  ```bash
  git clone https://github.com/NCI-CGR/TriosCompass_v2.git
  ```
+ Create new conda environment to run Snakemake workflow
  ```bash
  mamba env create -f TriosCompass_v2/environment.yaml
  ``` 
  
Folder structure of your workspace is as follows:
```bash
$WORKSPACE
└── TriosCompass_v2
```
 
### II. Inputs
#### Reference genome

User can put the reference genome any location under the folder $WORKSPACE and specify its relative location at config/config.yaml.  

For instance, we may put the hg38 human genome *Homo_sapiens_assembly38.fasta* under the folder $WORKSPACE/ref/. 

+ config/config.yaml
```yml
ref:
  sequence: "ref/Homo_sapiens_assembly38.fasta"
  build: "hg38"
```

---



---

#### The STR reference panel to call dnSTRs

We had developed joint call of dnSTR using both GangSTR and HipSTR, so we had prepared the same STR reference panels for both of the two programs.

Both GangSTR and HipSTR use STRs specified in the STR reference panel, but with the different formats.  One paper suggested that there is good consistence between HipSTR and GangSTR.
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



:bookmark: 
+ After manual curation, we chose to use *HipSTR* only for the dnSTR prediction.
+ The STR reference panel is specified in config/config.yaml
  + ref_panel: "ref/STR/hg38_ver13.hipstr_9.bed"

---
### III. Outputs

### IV. Run TrisCompass