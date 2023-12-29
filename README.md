# Snakemake workflow to process 8 trios
This is a revision of the workflow to process 8 Chernobyl trios, with the following new features:
+ Call *strelka* to refine DNM calls made by DV and GATK.
+ Use *JIGV* to generate IGV snapshots of the DNM variants.

There is no separate configure file for this workflow and all the settings have been coded directly within the file *Snakefile*.  The workflow has been run under the folder ***/data/DCEG_Chernobyl/NP0436-HE5/Chernobyl_data*** at Biowulf, where you may find all the details.

The diagram of the snakemake workflow is available [here](https://github.com/NCI-CGR/TriosCompass_v2/blob/8trios/Chernobyl_strelka_workflow_dag.pdf).

## Installation

The workflow has been tested under biowulf with snakemake version 7.3.7.  Snakemake can be installed under conda.

---

## Get started
### Set up the workspace
```bash
### Assume you have conda env dedicated to snakemake V7.3.7
conda activate snakemake

### Create folders used by the workflow (see /data/DCEG_Chernobyl/NP0436-HE5/Chernobyl_data at biowulf for the actual examples)
mkdir -p TMP logs bam ped_files ref

### Transfer bam files to the folder bam (not shown)
### Transfer ped files to the folder ped_files (not shown)
### Transfer the reference genome hg38 and WGS interval to the folder ref (not shown)
# ref/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list
# ref/Homo_sapiens_assembly38.fasta
```

### Run the wrapper script run_it_8trios.sh to launch the workflow at the slurm cluster at Biowulf.
```bash
sbatch -J 8trios -t 200:00:00 --export=ALL --mem=12g -p norm  --wrap='./run_it_8trios.sh '
```

--- 
## Recent changes

### Snakefile_samtools

A revised workflow [*Snakefile_samtools*](./Snakefile_samtools) is just added, together with the corresponding rulegraph of the workflow: [*Chernobyl_samtools_workflow_dag.pdf*](./Chernobyl_samtools_workflow_dag.pdf).  The major changes are to be highlighted, compared to [the original one](./Snakefile):
+ *samtools* is used to fix read group in the bam files and also save the output files in the cram format.
+ The downstream rules were also modifiled to use cram as inputs. 
+ The steps of *fix_dv_vcf* and *fix_gatk_vcf* in the original workflow is not needed any more.
+ The DNMs called by *DeepVariant* & *GATK HaplotypeCaller* are further phased by *unfazed*.
+ Besides, we had explored differet *slivar* filters in rule *call_dnm_strelka2*:
  + Add new PASS filter for *strelka* specifically
    + variant.FILTER == 'PASS'
  + Replace *(mom.AD[1] + dad.AD[1]) <= 5* with  <mark><i>(mom.AD[1]/(mom.AD[0]+mom.AD[1])) < 0.05 &&  (dad.AD[1]/(dad.AD[0]+dad.AD[1])) < 0.05</i></mark> for general improvement.
    + New rule is better to exclude low quality DNM candidates, for example, reads supports for the alternate allele present in parents with low read depth. 


### Phase DNMs using [WhatsHap](https://github.com/whatshap/whatshap/)
#### Install the latest version of WhatsHap
Currently, the version of WhatsHap available at Biowulf is V1.1, and the latest release is V2.1.  We installed the latest release using pip.

```bash
pip3 install --user whatshap
which whatshap
~/.local/bin/whatshap

whatshap --version
2.1
```

### Identify parental origin of DNMs using WhatsHap
There are two running modes in WhatsHap: individual and pedigree.  Pedigree mode is ideal to identify parental origin of variants in the child.  However, WhatsHap does not phase DNMs in the pedigree mode, as DNMs do not follow mendelian inheritance. Phasing DNMs is not supported by WhatsHap at present (see this [issue](https://github.com/whatshap/whatshap/issues/82) for the details). However, DNMs are phased in the individual mode.  Therefore, we managed to run WhatsHap in both individual mode and pedigree mode and extract parental original from the two outputs.

We first explored to run WhatsHap on the whole genome using the workflow [Snakefile_whatshap_WG](./Snakefile_whatshap_WG).  It turns out to be very slow: it took up to 7 days to process one trio family, and 2-3 days for each trio in average.  Then, we developed the workflow [Snakefile_whatshap_100K](./Snakefile_whatshap_100K) to phase a 200K window around each DNM:

![](img/Trios_workflow_rulegraph_whatshap_100K.png)

We developed [a Perl script](./scripts/extract_parental_origin.pl) to identify the parental origin, and the aggregated results are output to the file: output/phase_DNMs/{Trios_ID}.parental_origin.tab.  The output file is a headless tab delimited text file, with 9 columns including variant id of DNM, parental origin prediction.  

The perl script extracts the haplotype block containing the DNMs from the output of WhatsHap phasing child, and collects the phased variants of child in the block.  It also repeats the processing on the output of WhatsHap phasing trios.  The phase prediction from the latter is given as *F|M*, i.e, the first allele is the one inherited from the father (F) and the other in inherited from the mother (M). However, the phase prediction from the former is not certain, that is, it could be F|M or M|F. Ideally, the phases of the variants in each haplotype block from WhatsHap phasing child only should be either all identical (FM count) or all opposite (MF count) with those from WhatsHap phasing trios.  We simply assign parental origin as Not Determined (ND) if there is inconsistency (i.e., both FM count > 0 and MF count > 0). 

+ Example of the output of extract_parental_origin.pl
  
| VariantID                        | Parental Origin | Phase   | Haplotype Block Size | # Informative Sites | FM count | MF count | FM status | Parent origin change |
| -------------------------------- | --------------- | ---- | ---------------- | ------------------- | -------- | -------- | --------- | -------------------- |
| chr1:208903875:G:A               | ND              | 1\|0 | 36               | 0                   | 25       | 11       | ND        | ND=>ND               |
| chr2:12194398:A:C                | paternal        | 1\|0 | 2                | 2                   | 2        | 0        | FM        | paternal=>paternal   |
| chr2:14986429:A:T                | paternal        | 1\|0 | 5                | 5                   | 5        | 0        | FM        | paternal=>paternal   |
| chr2:50677044:C:T                | ND              | 0/1  | 0                | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr2:139840973:T:A               | paternal        | 0\|1 | 14               | 14                  | 0        | 14       | MF        | paternal=>paternal   |
| chr2:154745781:T:C               | ND              | 0/1  | 0                | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr2:203705880:T:TTCTTTC         | ND              | 0/1  | 0                | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr2:220449146:A:G               | ND              | 0/1  | 0                | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr3:14911256:G:C                | ND              | 0/1  | 0                | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr3:76410708:A:C                | ND              | 0/1  | 0                | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr3:185020246:C:T               | ND              | 0/1  | 0                | 0                   | 0        | 0        | ND        | ND=>ND               |
| chr3:194728646:G:GTGTGTGTGTGTGTC | maternal        | 0\|1 | 6                | 6                   | 6        | 0        | FM        | maternal=>maternal   |
| chr3:197879976:C:T               | paternal        | 0\|1 | 3                | 2                   | 0        | 1        | MF        | paternal=>paternal   |
| chr3:197879978:C:A               | paternal        | 0\|1 | 3                | 2                   | 0        | 1        | MF        | paternal=>paternal   |
| chr4:29657055:GC:G               | ND              | 0/1  | 0                | 0                   | 0        | 0        | ND        | ND=>ND               |


---

+ Descriptions for each output column.

| Column position | Column Header          | Description                                                   | Notes                                                                                                                |
| --------------- | ---------------------- | ------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| 1               | VariantID              | Variant ID of DNM                                             |                                                                                                                      |
| 2               | Parental Origin        | Predicted parental origin                                     | ND (not determined)                                                                                                  |
| 3               | Phase                     | Phased genotype of DNMs (from phasing child)                 | Parental origin cannot be determined from the output from phasing child only.                                        |
| 4               | Haplotype Block Size       | Size of the haplotype block (from phasing child)                  | Prediction from a large block is more reliable.                                                                      |
| 5               | # Informative Sites    | # informative sites in the haplotype block                        | Prediction is not reliable if # informative sites is 0.                                                             |
| 6               | FM count               | # phased variants supporting F\|M genotype                    |                                                                                                                      |
| 7               | MF count               | # phased variants supporting M\|F genotype                    |                                                                                                                      |
| 8               | FM status              | FM (for F\|M), MF (for M\|F), ND (not determined)             | FM status is ND if both FM and MF counts are bigger than 0.                                                          |
| 9               | Parental origin change | Prediction (wo phasing trios) => Prediction(wi phasing trios) | Generally, the prediction with phasing child only is consistent with that with both phasing child and phasing trios. |

