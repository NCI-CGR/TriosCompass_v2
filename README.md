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

## Output 
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