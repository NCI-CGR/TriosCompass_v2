<!-- omit in toc -->
# Call DNMs in STRs (Simple Tandem Repeats)


--- 

- [Introduction](#introduction)
- [Launch the Snakemake workflow](#launch-the-snakemake-workflow)
- [Results](#results)
  - [Working space at Biowulf](#working-space-at-biowulf)
  - [VizAln output](#vizaln-output)
    - [Example of VisAln realignment of the DNM (started at chr6:38571975) in the family *t0612*](#example-of-visaln-realignment-of-the-dnm-started-at-chr638571975-in-the-family-t0612)


--- 

## Introduction

We explored to call de novo mutations in STRs using HipSTR, GangSTR and MonSTR. The Snakemake workflow is as illustrated as the diagram below: 

![](./img/Trios_workflow_dag_STR3.png)

In particular, there are several features in this workflow:
+ We prepared identical STR reference panels for both HipSTR and GangSTR. 
  + The reference panel is derived from the one from GangSTR: hg38_ver13.bed. 
  + As HipSTR cannot take repeat unit size above 9, we restricted repeat unit size $\le$ 9 and generated the file: hg38_ver13.le9.bed.
  + And we formatted hg38_ver13.le9.bed to hg38_ver13.hipstr_9.bed, for the use with HipSTR. 
+ Both of the two STR reference panels were split into 400 chunks, so as to speed up the downstream analyses.
+ DNMs in STRs were called in two ways by each chunk: 
  + HipSTR + MonSTR
  + GangSTR + MonSTR
+ The chunks of both two DNM callings were merged in *merge_monstr*, respectively.
+ We selected DNMs as called by both of the two callings in *joint_STR*.
+ Finally, we explored two ways of visualization to facilitate manual curation.
  + By *JIGV* (in call_JIGV)
  + By [*VisAln* realignment](https://github.com/tfwillems/HipSTR/tree/master#alignment-visualization)

--- 

## Launch the Snakemake workflow

Users may launch workflow with the wrapper bash script [run_str3.sh](./data/run_str3.sh) at biowulf:

```bash
### Assume you have already installed Snakemake at conda env snakemake 
conda activate snakemake

snakemake --version
7.3.7

### Singularity is required to call MonSTR container
module load singularity

cd /data/DCEG_Chernobyl/NP0436-HE5/Chernobyl_data/STR_DNM2
sbatch -J STR3 -t 200:00:00 --export=ALL --mem=12g -p norm  --wrap='./run_str3.sh '
```
---

## Results
### Working space at Biowulf
+ /data/DCEG_Chernobyl/NP0436-HE5/Chernobyl_data/STR_DNM2
  + [Snakefile_STR3](./data/Snakefile_STR3) (Snakemake workflow)
  + [run_str3.sh](./data/run_str3.sh) (the wrapper bash script to launch Snakefile_STR3 at biowulf)
  + STR/
    + hg38_ver13.le9.bed
    + hg38_ver13.hipstr_9.bed
  + ped_files/
    + {family_id}.ped
  + output/ (output of the workflow)
    + merge_monstr/
      + gangstr.all_mutations.tab ([merged calls of GangSTR+MonSTR](https://github.com/gymreklab/STRDenovoTools#outprefixall_mutationstab))
      + gangstr.locus_summary.tab ([merged calls of GangSTR+MonSTR](https://github.com/gymreklab/STRDenovoTools#outprefixlocus_summarytab))
      + hipstr.all_mutations.tab ([merged calls of HipSTR+MonSTR](https://github.com/gymreklab/STRDenovoTools#outprefixall_mutationstab))
      + hipstr.locus_summary.tab ([merged calls of HipSTR+MonSTR](https://github.com/gymreklab/STRDenovoTools#outprefixlocus_summarytab))
    + joint_STR/
      + [gangstr_hipstsr.final.tab](./data/gangstr_hipstsr.final.tab) (DNMs called by both HipSTR & GangSTR + MonSTR).
    + call_JIGV
      + both_{family_id}.JIGV.html
    + vizaln/
      + {family_id}
        + variants
          + chr10_46899557.html
          + chr14_98831450.html

### VizAln output
VizAln is a tool from the HipSTR package to visualize the MGS reads in html page. It was designed to generate HTML page with random file name.  We had to make some revisions to apply it in our workflow. Different from JIGV, each HTML output from VizAln is for one DNM variant.  

Below is the outline of the VizAln output: 
```bash
tree output/vizaln/
output/vizaln/
├── t0311
│   ├── DONE
│   ├── t0311.aln.viz.gz
│   ├── t0311.aln.viz.gz.tbi
│   ├── t0311.hipstr.log
│   ├── t0311_recall.bed
│   ├── t0311.vcf.gz
│   ├── t0311.vcf.gz.tbi
│   └── variants
│       ├── chr10_46899557.dnm
│       ├── chr10_46899557.html
│       ├── chr14_98831450.dnm
│       ├── chr14_98831450.html
│       ├── chr15_61218675.dnm
│       ├── chr15_61218675.html
│       ├── chr5_71033413.dnm
│       ├── chr5_71033413.html
│       ├── chr6_84500344.dnm
│       └── chr6_84500344.html
...
└── t0612
    ├── DONE
    ├── t0612.aln.viz.gz
    ├── t0612.aln.viz.gz.tbi
    ├── t0612.hipstr.log
    ├── t0612_recall.bed
    ├── t0612.vcf.gz
    ├── t0612.vcf.gz.tbi
    └── variants
        ├── chr6_38571975.dnm
        └── chr6_38571975.html

```

#### Example of VisAln realignment of the DNM (started at chr6:38571975) in the family *t0612*

An example of VisAln realignment is available [here](https://htmlpreview.github.io/?https://github.com/NCI-CGR/TriosCompass_v2/blob/main/data/chr6_38571975.html).

In the family *t0612*, the identifiers of child, father and mother are SC260721, SC260720 and SC260670, respectively.
```bash
cat ped_files/t0612.ped 
t0612   SC260720        0       0       1       1
t0612   SC260670        0       0       2       1
t0612   SC260721        SC260720        SC260670        2       1
```

From the file [*gangstr_hipstsr.final.tab*](./data/gangstr_hipstsr.final.tab), we have: 
| chrom | pos      | period | child    | newallele | mutsize | child_gt | mat_gt | pat_gt |
| ----- | -------- | ------ | -------- | --------- | ------- | -------- | ------ | ------ |
| 6     | 38571975 | 2      | SC260721 | 13        | 2       | 13,13    | 11,13  | 11,11  |


It suggests there is an de novo tandem repeat mutation started at the position chr6:38571975 in the subject SC260721.  The unit length of repeat (period) is 2, and has 13 repeat units (newallele) in the child *SC260721*, which is 2 more than the reference genome (mutsize; i.e., the wild type genotype is [11,11]).  The MonSTR prediction indicated that the genotypes of child, father and mother are [13,13], [11,11] and [11,13], respectively. 

In the VizAln html page, we have genotypes of the family members marked in different ways.  For example, "SC260670: 0|4" means one allele is wild type and another is 4 bp insertion. Negative values stands for deletion (read [this](https://github.com/tfwillems/HipSTR/tree/master#alignment-visualization) for some details).  Therefore, "SC260721: 4|4", "SC260720: 0|0" and  "SC260670: 0|4" well match with the MonSTR predictions: [13,13], [11,11] and [11,13] for child, father and mother respectively.


:bookmark: Please note that the alignments from VizAln were arranged by the alphabet order of the sample identifiers.  