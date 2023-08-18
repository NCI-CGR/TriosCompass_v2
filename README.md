# TriosCompass_v2
This is a trios analysis workflow written in Snakemake.

## Introduction

This workflow takes fastq files and pedigree files as inupt and generate lists of de novo mutations (DNMs) for each tios in the end. Three variant callers are used here: DeepVariant, GATK HaplotypeCaller and Strelka.  The workflow is as described in [this diagram](./Trios_workflow_dag_CGR.pdf), and the details is availale in [this Snakefile](./Snakefile).

Briefly, there are 4 major steps in the workflow:
+ flastp: trim FASTQ reads;
+ fq2bam: align one or multipe pair FASTQ files (with trimed reads) to the reference genome;
+ gatk_markdup: mark duplicate reads;
+ Call DNMs
  1.  Variants called by deepvariant and merged by glnexus_dv by trios and then DNMs are called by slivar;
  2. DNMs called by the GATK pipeline: haplotypecaller, CombineGVCFs, genotypegvcf, CalculateGenotypePosteriors and finally followed by slivar filtering.
  3. Classify DNMs from Deepvariant and GATK in two groups: ones called by bother callers, and the others (i.e., called by only one of the two callers).
  4. Use those DNMs called by only one caller as the candidates to run Strelka and slivar.  
+ In the end, we have two sets of DNM candidates: 1) DNMs called by both DeepVariant and GATK; 2) DNMs called by (DeepVariant or GATK) and Strelka.  JIGV is utilized to generate IGV snapshots for bother of the two candidate sets.


---

## Installation

The workflow has been tested under biowulf with snakemake version 7.3.7.  Snakemake is installed under conda.

---
## Prepare input files
###  Manifest file and its schema specification

The manifest file is a csv file with typical NCI-CGR sample sheet format, with additional two columns (for the absolute paths to the paired Fastq files): R1, R2 (see https://github.com/NCI-CGR/TriosCompass_v2/blob/main/pep/manifest_fastq.csv).  The columns "CGR_ID", "INDEX", "FLOWCELL", "R1/2" are required and all the others are optional.  Of them, "CGR_ID", "INDEX", and "FLOWCELL" are used to define the read groups in the bam files, in the format:  
```
@RG\\tPL:ILLUMINA\\tID:{FLOWCELL}_{LANE}\\tSM:{CGF_ID}\\tPU:{CGF_ID}_{FLOWCELL}\\tLB:{CGF_ID}_{INDEX}
```

The format of the manifest file is specified by schemas/cgr_manifest_schema.yaml. The manifest file will be automatically validated at the beginning of the Snakemake workflow via the new Snakemake feature of PEP (protable encapsulated project).  Users may vist [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#configuring-scientific-experiments-via-peps) to learn more about this PEP feature. 

:bookmark: Note that one sample is allowed to have multiple fastq files from different flow cells.  The workflow will combine those fastq files and generate single bam file for each sample.

### Pedigree files
In additioanl to the fastq input files, another important input files are pedigress files, to specify relationship among subjects/samples.  It always consists of one child and its parents. For example, for a family with two kids, we need two peidgree files:

```bash
cat  new_cgr_pedfiles/t0679c1.ped 
t0679c1 SC742188        0       0       1       1
t0679c1 SC742277        0       0       2       1
t0679c1 SC742275        SC742188        SC742277        1       1

cat  new_cgr_pedfiles/t0679c2.ped 
t0679c2 SC742188        0       0       1       1
t0679c2 SC742277        0       0       2       1
t0679c2 SC742276        SC742188        SC742277        2       1

```  

:bookmark: Note that the samples are ordered in the pedigree file in this way: father, mother and kid.  It may affect the order of the sample tracks in the output html page of JIGV.

---

## Resource files
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

### Interval file
The interval file is to include the selected regions only to call DNMs (via bcftools).  We used the file hg38.wgs_interval.bed in this study, which was converted from *resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list*. The latter is part of [resource bundle hosted by the Broad Institute](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle). There are severals way to retrieve the interval file, for example, ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/wgs_calling_regions.hg38.interval_list.


+ Below is the bash command to convert the interval list to the bed format using *picard*.
```bash
picard IntervalListToBed -I ref/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list -O ref/hg38.wgs_interval.bed
```

### Configure file for fastq_screen
[Fastq_screen](https://stevenwingett.github.io/FastQ-Screen/) is used to identify likely sample contamination in the NGS data in this workflow. A configure file is required at ref/fastq_screen.abs.conf to have it work.

```conf
#DATABASE       Human   fasta_human/Homo_sapiens_assembly38_masked_GRC_exclusions.fasta
DATABASE        Yeast   /data/DCEG_Trios/new_cgr_data/TriosCompass_v2/ref/fasta_nonhuman/Yeast/Saccharomyces_cerevisiae.R64-1-1.fa
DATABASE        Ecoli   /data/DCEG_Trios/new_cgr_data/TriosCompass_v2/ref/fasta_nonhuman/E_coli/Ecoli.fa
DATABASE        PhiX    /data/DCEG_Trios/new_cgr_data/TriosCompass_v2/ref/fasta_nonhuman/PhiX/phi_plus_SNPs.fa
DATABASE        Lambda  /data/DCEG_Trios/new_cgr_data/TriosCompass_v2/ref/fasta_nonhuman/Lambda/Lambda.fa
DATABASE        Vectors /data/DCEG_Trios/new_cgr_data/TriosCompass_v2/ref/fasta_nonhuman/Vectors/Vectors.fa
DATABASE        Adapters        /data/DCEG_Trios/new_cgr_data/TriosCompass_v2/ref/fasta_nonhuman/Adapters/Contaminants.fa
DATABASE        Cow     /data/DCEG_Trios/new_cgr_data/TriosCompass_v2/ref/fasta_nonhuman/Cow/GCF_002263795.1_ARS-UCD1.2_genomic.fna
DATABASE        Pig     /data/DCEG_Trios/new_cgr_data/TriosCompass_v2/ref/fasta_nonhuman/Pig/GCF_000003025.6_Sscrofa11.1_genomic.fna
BWA     /usr/local/apps/bwa/0.7.17/bwa
```     
Accordingly, we need to provide non_human fasta sequence files a specified, so is the bwa binary file.

### Configure file for Snakemake workflow
Lastly but also importantly, a configure file for the Snakemake workflow is required.  The yaml file listed below serves as a both example and template.

```yaml
name: CGR_Trios_Data1

pep_version: 2.0.0
sample_table: manifest_fastq.csv
output_dir: "output"
input_bam_dir: "bam"

hg38_ref: "ref/Homo_sapiens_assembly38.fasta"
wgs_interval: "ref/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list"
ped_dir: "new_cgr_pedfiles"

# In CGR manifest file, CGF_ID + Flowcell should be unique
sample_modifiers:
  append:
    sample_name: "sn"
  derive:
    attributes: [sample_name]
    sources:
      sn: "{CGF_ID}_{FLOWCELL}"

```

---

## Get started 
Once the input files, resource data and the configure are ready, users may run the workflow using the wrapper shell script run_it_cgr.sh:

```bash
sbatch -J cgr_trios -t 200:00:00 --export=ALL --mem=12g -p norm  --wrap='./run_it_cgr.sh '
```

+ run_it_cgr.sh
  + :bookmark: Users may need change TMPDIR setting accordingly.
  + The profile workflow/profiles/biowulf is used in this example, which is tailored for the HPC system *biowulf* at NIH. 
```bash
#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH -o ${PWD}/snakemake.%j.out
#SBATCH -e ${PWD}/snakemake.%j.err


# conda activate snakemake
mkdir -p /data/DCEG_Trios/new_cgr_data/TriosCompass_v2/TMP
export TMPDIR=/data/DCEG_Trios/new_cgr_data/TriosCompass_v2/TMP

snakemake --skip-script-cleanup -k  --keep-incomplete --rerun-incomplete --profile workflow/profiles/biowulf --verbose -p --use-conda --jobs 400 --use-envmodules --latency-wait 600 -T 0 -s Snakefile
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
## Output
### Location of output files
The workflow processed the NGS data of Chernobyl Trios-Additional Families (see the details about the study in the [fogbugz: 31945: SR0436-012 Chernobyl Trios-Additional Families for 80x Germline WGS-ANALYSIS](https://cgr-bugz.nci.nih.gov/login?dest=%2ff%2fcases%2f31945)).  Both the Snakemake workflow and the output directory is under ***/data/DCEG_Trios/new_cgr_data/TriosCompass_v2*** at biowulf.

### DNM candidates
As mentioined in the introduction, there are two sets of DNM candidates: 
+ Called by DeepVariant and GATK.
  + output/GATK_DV/D_and_G.{TRIO_ID}.dnm.vcf.gz
+ Called by DeepVariant/GATK and Strelka.
  + output/slivar/strelka_{TRIO_ID}.dnm.vcf.gz

Here, {TRIO_ID} is the placeholder of the trio identifier, as specified in the pedigree file.  For instance, t0599c1, t0315c2 and etc.

### JIGV html pages
Accordinlgy, JIGV snapshots of both sets of DNMs can be found under output/call_JIGV/
+ D_and_G_{TRIO_ID}.JIGV.html
+ strelka_{TRIO_ID}.JIGV.html

### QC metrics
There are many metrics generated by NGS tools employed in this workflow:
+ fastp
+ fastq_screen
+ samtools flagstats
+ collectwgsmetrics (Picard)
+ collectmultiplemetrics (Parabricks)
+ GATK MarkDuplicate

Those output files can all be processed and summerized by *MultiQC*.  Below is a multiqc command for your example:  

```bash
multiqc --title QC --filename multiqc_report.html --outdir output_multiqc output/{collectmultiplemetrics,collectwgsmetrics,fastp,fastqc,fastq_screen,flagstats,gatk_markdup}  --interactive
|         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 4092/4092 
|            picard | Found 107 AlignmentSummaryMetrics reports
|            picard | Found 107 GcBiasMetrics reports
|            picard | Found 107 InsertSizeMetrics reports
|            picard | Found 107 MarkDuplicates reports
|            picard | Found 107 QualityByCycleMetrics reports
|            picard | Found 107 QualityScoreDistributionMetrics reports
|            picard | Found 107 QualityYieldMetrics reports
|            picard | Found 107 WgsMetrics reports
|          samtools | Found 107 flagstat reports
|          samtools | Found 107 idxstats reports
|      fastq_screen | Found 124 reports
|             fastp | Found 248 reports
|            fastqc | Found 248 reports
|           multiqc | Compressing plot data
|           multiqc | Previous MultiQC output found! Adjusting filenames..
|           multiqc | Use -f or --force to overwrite existing reports instead
|           multiqc | Report      : output_multiqc/multiqc_report_1.html
|           multiqc | Data        : output_multiqc/multiqc_report_data_1
|           multiqc | MultiQC complete


```
:bookmark: The json output from fastp needs to be named as "*fastp.json" to be recoganized by MultiQC.  We have modified the workflow accordingly.

### Summary report of new CGR run
In our latest trios analysis, 107 CGR samples from 40 trios were processed. We developed a Perl script [generate_summary_report.pl](./scripts/generate_summary_report.pl) to generate summary table in Xlxs format.

```bash
### remove the problematic trio t0588c1
ls new_cgr_pedfiles/*.ped  | grep -v t0588c1 > ped.lst

wc -l ped.lst
39 ped.lst

scripts/generate_summary_report.pl cgr_summary.xlsx ped.lst
```

+ [Summary table](./cgr_summary.xlsx)


| FamilyID | ID      | SampleID | FatherSampleID | MotherSampleID | Gender | DNM_Count_DG | DNM_Count_Strelka | JIGV_DG   | JIGV_Strelka |
| -------- | ------- | -------- | -------------- | -------------- | ------ | ------------ | ----------------- | --------- | ------------ |
| t0666    | t0666c1 | SC502245 | SC502256       | SC502234       | F      | 54           | 20                | JIGV HTML | JIGV HTML    |
| t0315    | t0315c1 | SC260714 | SC260727       | SC260729       | F      | 74           | 18                | JIGV HTML | JIGV HTML    |
| t0140    | t0140c1 | SC742196 | SC742286       | SC742197       | F      | 88           | 26                | JIGV HTML | JIGV HTML    |
| t0565    | t0565c2 | SC109437 | SC109373       | SC109368       | F      | 79           | 24                | JIGV HTML | JIGV HTML    |
| t0575    | t0575c1 | SC109390 | SC109395       | SC109405       | M      | 67           | 19                | JIGV HTML | JIGV HTML    |
| t0712    | t0712c1 | SC742313 | SC742314       | SC742315       | F      | 114          | 18                | JIGV HTML | JIGV HTML    |
| t0600    | t0600c2 | SC109495 | SC109514       | SC109500       | F      | 93           | 14                | JIGV HTML | JIGV HTML    |
| t0565    | t0565c1 | SC109438 | SC109373       | SC109368       | F      | 68           | 14                | JIGV HTML | JIGV HTML    |
| t0705    | t0705c1 | SC742298 | SC742301       | SC742299       | F      | 101          | 15                | JIGV HTML | JIGV HTML    |
| t0750    | t0750c1 | SC736755 | SC736760       | SC736736       | F      | 93           | 27                | JIGV HTML | JIGV HTML    |
| t0707    | t0707c1 | SC742305 | SC742306       | SC742307       | F      | 61           | 11                | JIGV HTML | JIGV HTML    |
| t0592    | t0592c1 | SC109501 | SC109499       | SC109498       | F      | 62           | 31                | JIGV HTML | JIGV HTML    |
| t0679    | t0679c1 | SC742275 | SC742188       | SC742277       | M      | 76           | 20                | JIGV HTML | JIGV HTML    |
| t0007    | t0007c1 | SC499427 | SC499423       | SC499428       | M      | 67           | 14                | JIGV HTML | JIGV HTML    |
| t0739    | t0739c1 | SC736720 | SC736795       | SC736726       | M      | 64           | 16                | JIGV HTML | JIGV HTML    |
| t0575    | t0575c2 | SC109406 | SC109395       | SC109405       | M      | 60           | 12                | JIGV HTML | JIGV HTML    |
| t0243    | t0243c1 | SC736787 | SC736702       | SC736772       | F      | 99           | 28                | JIGV HTML | JIGV HTML    |
| t0058    | t0058c2 | SC108472 | SC109353       | SC109336       | M      | 86           | 58                | JIGV HTML | JIGV HTML    |
| t0042    | t0042c1 | SC253873 | SC253877       | SC253872       | M      | 97           | 27                | JIGV HTML | JIGV HTML    |
| t0766    | t0766c2 | SC736756 | SC736703       | SC730933       | F      | 113          | 31                | JIGV HTML | JIGV HTML    |
| t0765    | t0765c1 | SC736700 | SC736738       | SC736739       | M      | 115          | 19                | JIGV HTML | JIGV HTML    |
| t0766    | t0766c1 | SC742179 | SC736703       | SC730933       | M      | 82           | 12                | JIGV HTML | JIGV HTML    |
| t0315    | t0315c2 | SC260715 | SC260727       | SC260729       | M      | 80           | 18                | JIGV HTML | JIGV HTML    |
| t0483    | t0483c2 | SC742320 | SC742219       | SC742321       | F      | 125          | 19                | JIGV HTML | JIGV HTML    |
| t0623    | t0623c1 | SC253881 | SC253893       | SC253894       | F      | 52           | 20                | JIGV HTML | JIGV HTML    |
| t0760    | t0760c1 | SC736729 | SC736824       | SC736803       | F      | 106          | 26                | JIGV HTML | JIGV HTML    |
| t0703    | t0703c1 | SC742295 | SC742296       | SC742297       | F      | 77           | 18                | JIGV HTML | JIGV HTML    |
| t0271    | t0271c1 | SC109417 | SC109426       | SC109450       | M      | 144          | 25                | JIGV HTML | JIGV HTML    |
| t0679    | t0679c2 | SC742276 | SC742188       | SC742277       | F      | 105          | 19                | JIGV HTML | JIGV HTML    |
| t0693    | t0693c2 | SC742220 | SC742221       | SC742222       | F      | 90           | 25                | JIGV HTML | JIGV HTML    |
| t0599    | t0599c1 | SC736788 | SC736742       | SC736743       | M      | 106          | 20                | JIGV HTML | JIGV HTML    |
| t0617    | t0617c1 | SC260701 | SC260705       | SC260697       | M      | 131          | 33                | JIGV HTML | JIGV HTML    |
| t0600    | t0600c1 | SC109494 | SC109514       | SC109500       | M      | 58           | 9                 | JIGV HTML | JIGV HTML    |
| t0749    | t0749c1 | SC736820 | SC736759       | SC736701       | F      | 90           | 21                | JIGV HTML | JIGV HTML    |
| t0022    | t0022c1 | SC742290 | SC742217       | SC742218       | M      | 64           | 12                | JIGV HTML | JIGV HTML    |
| t0280    | t0280c1 | SC109512 | SC109516       | SC109507       | M      | 89           | 12                | JIGV HTML | JIGV HTML    |
| t0594    | t0594c1 | SC736817 | SC736753       | SC736754       | F      | 110          | 17                | JIGV HTML | JIGV HTML    |
| t0058    | t0058c1 | SC109341 | SC109353       | SC109336       | F      | 73           | 42                | JIGV HTML | JIGV HTML    |
| t0209    | t0209c1 | SC742216 | SC742271       | SC742272       | F      | 87           | 21                | JIGV HTML | JIGV HTML    |


:bookmark: In the Excel file, the last two columns of "JIGV HTML" are built in with HTTP links to the local HTML pages under output/call_JIGV.  Therefore, the HTTP links will not work properly unless the report file is located in the right position. 
```bash
cgr_summary.xlsx
output/call_JIGV/
├── D_and_G_t0007c1.JIGV.html
├── D_and_G_t0022c1.JIGV.html
...
├── D_and_G_t0766c1.JIGV.html
├── D_and_G_t0766c2.JIGV.html
├── strelka_t0007c1.JIGV.html
├── strelka_t0022c1.JIGV.html
...
├── strelka_t0766c1.JIGV.html
└── strelka_t0766c2.JIGV.html
```

Alternatively, users may choose to modify scripts/generate_summary_report.pl to generate HTTP links in different ways. 