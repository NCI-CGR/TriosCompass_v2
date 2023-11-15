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