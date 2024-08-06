# Run TriosCompass (public version) for the 8-trios data set

## Create the working space
```bash
mkdir -p /data/DCEG_Trios/ChernobylTrios/test_8trios

cd /data/DCEG_Trios/ChernobylTrios/test_8trios
```

## Prepare input files
+ ref/
  + STR/ (for reference panels to call STR)
+ ped/
+ TriosCompass_v2/
+ a wrapper script
  + run_fullbam.sh

```bash
cd /data/DCEG_Trios/ChernobylTrios/test_8trios

mkdir ref
cp ../dev/ref/{hg38.wgs_interval.bed,Homo_sapiens_assembly38.fasta} ref/
cp -r ../dev/ref/STR ref/STR

tree ref
# ref
# ├── hg38.wgs_interval.bed
# ├── Homo_sapiens_assembly38.fasta
# └── STR
#     ├── GRCh38GenomicSuperDup.bed.gz
#     ├── GRCh38GenomicSuperDup.bed.gz.tbi
#     ├── hg38_ver13.bed
#     ├── hg38_ver13.hipstr_9.bed
#     └── hg38_ver13.le9.bed

mkdir fixed_rg

mv ../8trios/output/fixed-rg/* fixed_rg/

cat  ../8trios/8trios.lst | parallel 'cp ../8trios/ped_files/{}.ped ped/'

tree ped
# ped
# ├── t0311.ped
# ├── t0321.ped
# ├── t0334.ped
# ├── t0341.ped
# ├── t0402.ped
# ├── t0474.ped
# ├── t0605.ped
# └── t0612.ped

git clone git@github.com:NCI-CGR/TriosCompass_v2.git
git checkout general_use

cp ../dev/run_fullbam.sh  . 

cat run_fullbam.sh
# #!/bin/bash
# #SBATCH --time=200:00:00
# #SBATCH -o ${PWD}/snakemake.%j.out
# #SBATCH -e ${PWD}/snakemake.%j.err

# ### or make sure singularity is available in alternative way
# module load singularity

# mkdir -p TMP
# export TMPDIR=TMP

# snakemake  --profile TriosCompass_v2/workflow/profiles/fullbam(/data/zhuw10/conda_envs/TriosCompass
```

---

## Prepare configure files
+ 8trios_bam.csv
+ TriosCompass_v2/workflow/profiles/fullbam/config.yaml
  + TriosCompass_v2/config/fullbam_config.yaml
+ TriosCompass_v2/config/fullbam_config.yaml
  + pepfile: "config/fullbam_pep.yaml"
  + output_dir: "output" :+1:
+ config/fullbam_pep.yaml :+1:
  + sample_table: 8trios_bam.csv

```bash
echo "SAMPLE_ID,BAM" > TriosCompass_v2/config/8trios_bam.csv
ls  fixed_rg/*.bam | parallel --dry-run "{/.},{}" >> TriosCompass_v2/config/8trios_bam.csv

cat  TriosCompass_v2/config/8trios_bam.csv
# SAMPLE_ID,BAM
# SC056068,fixed_rg/SC056068.bam
# SC056073,fixed_rg/SC056073.bam
# SC056082,fixed_rg/SC056082.bam
# SC074198,fixed_rg/SC074198.bam
# SC074199,fixed_rg/SC074199.bam
# SC074201,fixed_rg/SC074201.bam
# SC109440,fixed_rg/SC109440.bam
# SC109446,fixed_rg/SC109446.bam
# SC109448,fixed_rg/SC109448.bam
# SC109457,fixed_rg/SC109457.bam
# SC109459,fixed_rg/SC109459.bam
# SC109515,fixed_rg/SC109515.bam
# SC260670,fixed_rg/SC260670.bam
# SC260700,fixed_rg/SC260700.bam
# SC260706,fixed_rg/SC260706.bam
# SC260709,fixed_rg/SC260709.bam
# SC260720,fixed_rg/SC260720.bam
# SC260721,fixed_rg/SC260721.bam
# SC501095,fixed_rg/SC501095.bam
# SC501096,fixed_rg/SC501096.bam
# SC501105,fixed_rg/SC501105.bam
# SC501108,fixed_rg/SC501108.bam
# SC501110,fixed_rg/SC501110.bam
# SC501111,fixed_rg/SC501111.bam

```

--- 

## Dry run 
+ :bookmark: The conda env *TriosCompassV2* was created using the file TriosCompass_v2/environment.yaml

```bash
conda activate TriosCompassV2

module load singularity

snakemake  -np --profile TriosCompass_v2/workflow/profiles/fullbam

Job stats:
job                        count    min threads    max threads
-----------------------  -------  -------------  -------------
aggregate_phase                8              1              1
aggregate_visaln               8              1              1
all                            1              1              1
bedfile_for_hipstr_call        8              1              1
call_JIGV                      8              1              1
call_dnm_dv                    8              1              1
call_dnm_gatk                  8              1              1
collectmultiplemetrics        24              1              1
collectwgsmetrics             24              1              1
deepvariant_pb                24             48             48
dnm_vcf                        8              1              1
dumpstr_call                 400              1              1
dumpstr_locus                400              1              1
flagstat                      24             16             16
gatk_cgp                       8             16             16
gatk_combine_gvcf              8             16             16
gatk_genotype_gvcf_pb          8             16             16
gatkhc_pb                     24             24             24
genome_dict                    1              1              1
genome_faidx                   1              1              1
glnexus_dv                     8              8              8
hipstr                       400              1              1
hipstr_recall                  8              1              1
merge_DV_GATK                  8              1              1
merge_monstr                   1              1              1
merge_ped                      1              1              1
monstr                       400              1              1
monstr_filter                  1              1              1
scatter_chr_pos                8              1              1
scatter_dnms                   8              1              1
split_bed_hipstr               1              1              1
vcf_index                    400              1              1
total                       2247              1             48

snakemake  -np --rulegraph --profile TriosCompass_v2/workflow/profiles/fullbam | grep -v "^Detect" | dot -Tpng > TrioCompassV2_dag_8trios.png


```


--- 

### Launch TriosCompass workflow using slurm
```bash
sbatch -J 8trios -t 200:00:00 --export=ALL --mem=12g -p norm  --wrap='./run_fullbam.sh '
```
