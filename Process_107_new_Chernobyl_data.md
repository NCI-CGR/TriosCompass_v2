<!-- omit in toc -->
# Processing 107 new CGR data using TriosCompass

---
- [Introduction](#introduction)
- [The pipeline](#the-pipeline)
  - [Configuration and input files](#configuration-and-input-files)
    - [Fastq input files](#fastq-input-files)
    - [Manifest file and its schema specification](#manifest-file-and-its-schema-specification)
    - [Pedigree files](#pedigree-files)
  - [Configure file for Snakemake workflow](#configure-file-for-snakemake-workflow)
- [Get started](#get-started)
- [Workflow components for Fastq input files](#workflow-components-for-fastq-input-files)
- [Output](#output)
  - [Location of output files](#location-of-output-files)
  - [DNM candidates](#dnm-candidates)
  - [JIGV html pages](#jigv-html-pages)
  - [QC metrics](#qc-metrics)
    - [Output of MultiQC](#output-of-multiqc)
    - [Problematic samples](#problematic-samples)
  - [Summary report of new CGR run](#summary-report-of-new-cgr-run)


---

## Introduction

The 107 WGS data are from Chernobyl Trios-Additional Families (see the details about the study in the [fogbugz: 31945: SR0436-012 Chernobyl Trios-Additional Families for 80x Germline WGS-ANALYSIS](https://cgr-bugz.nci.nih.gov/login?dest=%2ff%2fcases%2f31945)). 

This workflow *TriosCompass* has been generally introduced [here](https://github.com/NCI-CGR/TriosCompass_v2/blob/main/README.md).  In this study, we had call DNMs using two approaches:
1. The candidates called by at least two out of the 3 callers: DeepVariant, GATK HaplotypeCaller and Strelka (see [Snakefile](./Snakefile));
![](./img/TrioCompass_Strelka_March2024_dag.png)
2. The candidates called by both DeepVariant and GATK HaplotypeCaller (see [Snakefile_CGRv2](./Snakefile_CGRv2)).  
![](./img/TrioCompass_CGRv2_March2024_dag.png)

We took the second one as final approach to have higer precision at the trade-off of the recall.  

---

## The pipeline
### Configuration and input files
#### Fastq input files
Fastq files were transferred from /mnt/nfs/gigantor/ifs/DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_
Analysis/Data/ (under the CGR T-drive) to /data/DCEG_Trios/new_cgr_data/fastq at biowulf via globus. 


```bash
tree /data/DCEG_Trios/new_cgr_data/fastq | head -n 40
/data/DCEG_Trios/new_cgr_data/fastq
├── 220817_A00423_0187_AH2JYFDSX5
│   └── CASAVA
│       └── L1
│           └── Project_NP0436-HE5
│               ├── Sample_SC074219-CTGATCGT-GCGCATAT
│               │   ├── SC074219_CTGATCGT-GCGCATAT_L001_R1_001.fastq.gz
│               │   └── SC074219_CTGATCGT-GCGCATAT_L001_R2_001.fastq.gz
│               ├── Sample_SC108472-TACGCTAC-CGTGTGAT
│               │   ├── SC108472_TACGCTAC-CGTGTGAT_L001_R1_001.fastq.gz
│               │   └── SC108472_TACGCTAC-CGTGTGAT_L001_R2_001.fastq.gz
│               ├── Sample_SC109336-CCTTGATC-GATGGAGT
│               │   ├── SC109336_CCTTGATC-GATGGAGT_L001_R1_001.fastq.gz
│               │   └── SC109336_CCTTGATC-GATGGAGT_L001_R2_001.fastq.gz
│               ├── Sample_SC109341-TGTGACTG-AGCCTATC
│               │   ├── SC109341_TGTGACTG-AGCCTATC_L001_R1_001.fastq.gz
│               │   └── SC109341_TGTGACTG-AGCCTATC_L001_R2_001.fastq.gz
│               ├── Sample_SC109353-GAGCAGTA-TGAGCTGT
│               │   ├── SC109353_GAGCAGTA-TGAGCTGT_L001_R1_001.fastq.gz
│               │   └── SC109353_GAGCAGTA-TGAGCTGT_L001_R2_001.fastq.gz
│               ├── Sample_SC109368-TGGTAGCT-AGGTGTTG
│               │   ├── SC109368_TGGTAGCT-AGGTGTTG_L001_R1_001.fastq.gz
│               │   └── SC109368_TGGTAGCT-AGGTGTTG_L001_R2_001.fastq.gz
│               ├── Sample_SC109373-GTTGACCT-CAGTGCTT
│               │   ├── SC109373_GTTGACCT-CAGTGCTT_L001_R1_001.fastq.gz
│               │   └── SC109373_GTTGACCT-CAGTGCTT_L001_R2_001.fastq.gz
│               └── Sample_SC109390-CACCTGTT-AGGATAGC
│                   ├── SC109390_CACCTGTT-AGGATAGC_L001_R1_001.fastq.gz
│                   └── SC109390_CACCTGTT-AGGATAGC_L001_R2_001.fastq.gz
├── 220817_A00423_0188_BH2JVVDSX5
│   └── CASAVA
│       └── L1
│           └── Project_NP0436-HE5
│               ├── Sample_SC742298-CTAGGCAT-ACAGAGGT
│               │   ├── SC742298_CTAGGCAT-ACAGAGGT_L001_R1_001.fastq.gz
│               │   └── SC742298_CTAGGCAT-ACAGAGGT_L001_R2_001.fastq.gz
│               ├── Sample_SC742299-TACATCGG-TCCACGTT
│               │   ├── SC742299_TACATCGG-TCCACGTT_L001_R1_001.fastq.gz
│               │   └── SC742299_TACATCGG-TCCACGTT_L001_R2_001.fastq.gz
│               ├── Sample_SC742301-TGTGCGTT-TGAGACGA
...
```

####  Manifest file and its schema specification
The locations of each pair of fastq are specified in a manifest file.  The manifest file is a csv file with typical NCI-CGR sample sheet format, with additional two columns (for the absolute paths to the paired Fastq files): R1, R2 (see https://github.com/NCI-CGR/TriosCompass_v2/blob/main/pep/manifest_fastq.csv).  The columns "CGR_ID", "INDEX", "FLOWCELL", "R1/2" are required and all the others are optional.  Of them, "CGR_ID", "INDEX", and "FLOWCELL" are used to define the read groups in the bam files, in the format:  
```
@RG\\tPL:ILLUMINA\\tID:{FLOWCELL}_{LANE}\\tSM:{CGF_ID}\\tPU:{CGF_ID}_{FLOWCELL}\\tLB:{CGF_ID}_{INDEX}
```

The format of the manifest file is specified by schemas/cgr_manifest_schema.yaml. The manifest file will be automatically validated at the beginning of the Snakemake workflow via the new Snakemake feature of PEP (protable encapsulated project).  Users may vist [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#configuring-scientific-experiments-via-peps) to learn more about this PEP feature. 

:bookmark: Note that one sample is allowed to have multiple fastq files from different flow cells.  The workflow will combine those fastq files and generate single bam file for each sample.

#### Pedigree files
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
---
### Configure file for Snakemake workflow
Lastly but also importantly, a configure file for the Snakemake workflow is required:
```yaml
name: CGR_Trios_Data1

pep_version: 2.0.0
sample_table: manifest_fastq.csv
output_dir: "output"
input_bam_dir: "bam"

hg38_ref: "ref/Homo_sapiens_assembly38.fasta"
wgs_interval: "ref/hg38.wgs_interval.bed"
ped_dir: "new_cgr_pedfiles"

### Settings for dnSTR calling
split_total: 400
gangstr_panel: "STR/hg38_ver13.le9.bed"
hipstr_panel: "STR/hg38_ver13.hipstr_9.bed"
hipstr_filters: " --min-span-coverage 3 --min-supp-reads 3 "
gangstr_filters: " --max-perc-encl-parent 0.05 --min-encl-match 0.9 --min-total-encl 10 --gangstr "
dup_reg: "STR/GRCh38GenomicSuperDup.bed.gz" # come with GRCh38GenomicSuperDup.bed.gz.tbi 

### WhatsHap
phase_window: 10000

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
cd  /data/DCEG_Trios/new_cgr_data/TriosCompass_v2
conda activate snakemake

sbatch -J cgrv2 -t 200:00:00 --export=ALL --mem=12g -p norm  --wrap='./run_it_cgrv2.sh '
```

+ run_it_cgrv2.sh
  + :bookmark: Users may need change TMPDIR setting accordingly.
  + The profile workflow/profiles/biowulf is used in this example, which is tailored for the HPC system *biowulf* at NIH. 
```bash
#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH -o ${PWD}/snakemake.%j.out
#SBATCH -e ${PWD}/snakemake.%j.err


# conda activate snakemake
mkdir -p TMP
export TMPDIR=TMP

module load singularity 

snakemake --skip-script-cleanup -k  --keep-incomplete --rerun-incomplete --profile workflow/profiles/biowulf --verbose -p --use-conda --jobs 400 --use-singularity --use-envmodules --latency-wait 600 -T 0 -s Snakefile_CGRv2
```

---
## Workflow components for Fastq input files

This workflow starts the process from the fastq files, which is different from the one for the old Chernobyl samples.  The latter takes the exsiting bam files as the starting point.  Therefore, there are several unique components in this workflow:
+ fastp
+ fastqc
+ fastq_screen
+ fq2bam
+ gatk_markdup

Accordingly, we need provide one additional configure file for [Fastq_screen](https://stevenwingett.github.io/FastQ-Screen/).

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


#### Output of MultiQC
We have generated two sets of MultiQC report files under the folder output_multiqc/: 
+ multiqc_report.html (static mode)
+ multiqc_report_1.html (interactive mode)

#### Problematic samples
+ SC074219 (t0450c2) is a singleton in this analysis and the data of its parents are processed in the previous study.
+ ***SC736795*** has coverage about 28X, which is much lower than the others.
+	SC108472 has many unmapped reads, origin from bateria and viruses commonly found in the upper respiratory tractor.  It is confimred that the sample was procured from saliva.
+ In the trio t0588c1 (kid: SC109409; father: SC109418; mother: SC109419), the sample ***SC109418*** lacks of the expected relatedness with SC109409.

---
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

---


