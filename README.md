<!-- omit in toc -->
# Process 340 Old Chernobyl sample using TriosCompass

---
- [Introduction](#introduction)
- [Get started](#get-started)
  - [Set up the workspace](#set-up-the-workspace)
  - [Run the wrapper script launch the workflow at the slurm cluster at Biowulf.](#run-the-wrapper-script-launch-the-workflow-at-the-slurm-cluster-at-biowulf)
- [Prepare 341 bam files](#prepare-341-bam-files)
  - [Transfer existing bam files of Set 1-3 from Object Store at Biowulf](#transfer-existing-bam-files-of-set-1-3-from-object-store-at-biowulf)
  - [Transfer Set4 bam files from T drive at nci-cgr to biowulf](#transfer-set4-bam-files-from-t-drive-at-nci-cgr-to-biowulf)
  - [Check integrity of the transferred bam files](#check-integrity-of-the-transferred-bam-files)
  - [Remove redundant copies of the bam files](#remove-redundant-copies-of-the-bam-files)
  - [Identify 15 bam files missed (340-325=15)](#identify-15-bam-files-missed-340-32515)
  - [Transfer 12 bam files (aligned to hg19) from T-drive at nci-cgr to biowufl](#transfer-12-bam-files-aligned-to-hg19-from-t-drive-at-nci-cgr-to-biowufl)
  - [Realign 12 bam files to hg38](#realign-12-bam-files-to-hg38)
  - [Finally, compile a list of 341 bam files](#finally-compile-a-list-of-341-bam-files)
- [Configuration of the workflow](#configuration-of-the-workflow)
  - [Notes for dnSTR calling](#notes-for-dnstr-calling)
- [MultiQC](#multiqc)
- [Post-analysis](#post-analysis)

---
## Introduction

We applied TriosCompass on the 340 old Chernobyl data, together with one new sample SC074219 (t0450c2) from recently sequenced at the CGR lab.

The pipeline is almost idential to the one to process the 107 new CGR sample, and we highlighted major differences in the following sections. 

---

## Get started
### Set up the workspace
```bash
cd /data/DCEG_Trios/ChernobylTrios
git clone git@github.com:NCI-CGR/TriosCompass_v2.git

cd /data/DCEG_Trios/ChernobylTrios/TriosCompass_v2

### Assume you have conda env dedicated to snakemake V7.3.7
conda activate snakemake

### Create folders used by the workflow (see /data/DCEG_Chernobyl/NP0436-HE5/Chernobyl_data at biowulf for the actual examples)
mkdir -p logs bchernobyl_ped2/ ref/ STR/


### Transfer bam files to the folder bam (not shown)
### Transfer ped files to the folder ped_files (not shown)
### Transfer the reference genome hg38 and WGS interval to the folder ref (not shown)
# ref/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list
# ref/Homo_sapiens_assembly38.fasta
```

### Run the wrapper script launch the workflow at the slurm cluster at Biowulf.
```bash
sbatch -J chernobyl -t 200:00:00 --export=ALL --mem=12g -p norm  --wrap='./run_it_chernobyl.sh '
```

--- 
## Prepare 341 bam files
### Transfer existing bam files of Set 1-3 from Object Store at Biowulf
```bash
cd /data/DCEG_Trios/ChernobylTrios

### Get a list of bam files in the vaults DCEG_Trios and DCEG_Chernobyl
obj_ls -v DCEG_Trios > DCEG_Trios_data.lst
obj_ls -v DCEG_Chernobyl > DCEG_Chernobyl_data.lst

grep -c  "bam$" *.lst
DCEG_Chernobyl_data.lst:123
DCEG_Trios_data.lst:107

# obj_get -v DCEG_Chernobyl -D /vf/users/DCEG_Chernobyl/NP0436-HE5/Chernobyl_data/bam 0006/SC074198.pb.realn.bam 
pcregrep -o1 ' ([\S]+.bam$)'  *.lst | sort -u | wc -l
# 230

### At least 6 bam files are duplicated
pcregrep -o1 '/(SC[0-9]{6}).+.bam$'  *.lst | sort -u | wc -l
# 224

### Prepare swarm commands to download all and figure out the duplicated later
mkdir -p bams/{DCEG_Chernobyl,DCEG_Trios} 

### Use --progressbar --checksum --verbose
# cat logs/test_f1_6951527_0.e
# No checksum found for 0002/SC056044.pb.realn.bam - was it uploaded with the -c argument to obj_put?
# Checksum comparison will be skipped.
pcregrep -o1 ' ([\S]+.ba[mi]?$)'  DCEG_Chernobyl_data.lst | sort -u | parallel --dry-run 'obj_get --progressbar --checksum --verbos
e -v DCEG_Chernobyl -D bams/DCEG_Chernobyl/ {}' > DCEG_Chernobyl_data.swarm

pcregrep -o1 ' ([\S]+.ba[mi]?$)'  DCEG_Trios_data.lst | sort -u | parallel --dry-run 'obj_get --progressbar --checksum --verbose -v
 DCEG_Trios -D bams/DCEG_Trios/ {}' > DCEG_Trios_data.swarm

wc -l DCEG*.swarm
  249 DCEG_Chernobyl_data.swarm
  214 DCEG_Trios_data.swarm
  463 total

### There is some problem in pulling data from object store, have to split the processing into several smaller batches
split --numeric-suffixes=0 -n l/3  --suffix-length=2  --additional-suffix=".swarm" DCEG_Chernobyl_data.swarm "splitted_DCEG_Chernob
yl_data."

split --numeric-suffixes=0 -n l/2  --suffix-length=2  --additional-suffix=".swarm" DCEG_Trios_data.swarm "splitted_DCEG_Trios_data.
"

grep -c ".bam$" *.swarm
DCEG_Chernobyl_data.swarm:123
DCEG_Trios_data.swarm:107
splitted_DCEG_Chernobyl_data.00.swarm:42
splitted_DCEG_Chernobyl_data.01.swarm:40
splitted_DCEG_Chernobyl_data.02.swarm:41
splitted_DCEG_Trios_data.00.swarm:54
splitted_DCEG_Trios_data.01.swarm:53

### Download
swarm -J chernobyl_p0 --logdir logs --time=100:00:00 splitted_DCEG_Chernobyl_data.00.swarm
swarm -J chernobyl_p1 --logdir logs --time=100:00:00 splitted_DCEG_Chernobyl_data.01.swarm
swarm -J chernobyl_p2 --logdir logs --time=100:00:00 splitted_DCEG_Chernobyl_data.02.swarm


swarm -J trios_p0 --logdir logs --time=100:00:00 splitted_DCEG_Trios_data.00.swarm
swarm -J trios_p1 --logdir logs --time=100:00:00 splitted_DCEG_Trios_data.01.swarm

```

###  Transfer Set4 bam files from T drive at nci-cgr to biowulf
```bash
### @ Biowulf
mkdir -p /data/DCEG_Trios/ChernobylTrios/bams/set4

### @CGR
cd /DCEG/Scimentis/DNM/data/BATCH2_b38

ls -al *.bam |wc -l
126 

### transfer via rsync
screen
rsync -av --progress *.bai *.bam helix.nih.gov:/data/DCEG_Trios/ChernobylTrios/bams/set4
```

### Check integrity of the transferred bam files
```bash
cd /data/DCEG_Trios/ChernobylTrios

du -hs bams/*
25T     bams/DCEG_Chernobyl
21T     bams/DCEG_Trios
12T     bams/set4

module load samtools


mkdir quickcheck
find bams -name "*.bam" | parallel -j 2  "samtools quickcheck -v {} > quickcheck/{/.}.bad_bams.fofn   && echo 'all ok' || echo 'som
e files failed check, see quickcheck/{/.}.bad_bams.fofn' "
```

### Remove redundant copies of the bam files
```bash
cd /data/DCEG_Trios/ChernobylTrios/bams

find . -name "*.bam" | pcregrep -o1 '/(SC[0-9]{6}).+.bam$' | sort | uniq -D| sort -u > bam_multi.lst

### compile a list of redundant copies of the bam files to review manually
mkdir dup
cat bam_multi.lst | parallel ' grep {} bam.lst | xargs ls -al > dup/{}.lst ' 

### remove all 21 bam files under ./DCEG_Trios/0009_calling_reruns
# and 7 additional ones
rm -fr ./DCEG_Trios/0009_calling_reruns ./DCEG_Trios/reruns0001/SC056063.pb.realn.bam* ./DCEG_Chernobyl/0003/SC056107.pb.realn.bam* ./DCEG_Chernobyl/0003/SC074124.pb.realn.bam* ./DCEG_Trios/batch0005/SC074132.pb.realn.bam* ./DCEG_Chernobyl/0004/SC074135.pb.realn.bam* ./DCEG_Chernobyl/0004/SC074137.pb.realn.bam* ./DCEG_Chernobyl/0004/SC074141.pb.realn.bam*

### Review the list of bam files and we got 325 unique bam files
find . -name "*.bam" | pcregrep -o1 '/(SC[0-9]{6}).*.bam$' | sort -u |wc -l
# 325
```

### Identify 15 bam files missed (340-325=15)
+ Found 3 of them in the object store
+ Keep the remaining 12 in the file *12missing.lst*
```bash
cd /data/DCEG_Trios/ChernobylTrios
mkdir -p manifest

scp Files/Trios_Master_CCAD_Biowulf_merged.csv helix:/data/DCEG_Trios/ChernobylTrios/manifest/

### Replace ^M (ctrl-v/ctrl-M) to new line
# https://www.unix.com/unix-for-dummies-questions-and-answers/41277-how-convert-m-appearing-end-line-unix-newline.html
sed 's/^M/\
/g' manifest/Trios_Master_CCAD_Biowulf_merged.csv > Trios_Master_CCAD_Biowulf_merged.clean.csv

### There are 339 (not 340 files) in CSV file
# wc count the line characcter
wc -l Trios_Master_CCAD_Biowulf_merged.clean.csv 
340 Trios_Master_CCAD_Biowulf_merged.clean.csv

### but 340 is listed in Excel file
wc -l trios340.lst 
340 trios340.lst

### Identical! 
cut -d, -f 1 Trios_Master_CCAD_Biowulf_merged.clean.csv  | sort -u > f1
sort -u trios340.lst > f2
comm -23 f2 f1

### list the samples we got already
find bams -name "*.bam" |  pcregrep -o1 '/(SC[0-9]{6}).*.bam$' | sort -u > f3

comm -23 f1 f3 | sed -n '1!p' > 15missing.lst
wc -l 15missing.lst

cat 15missing.lst 
SC056076 x
SC074149 x
SC074157 N
SC074166 N
SC074169 N
SC074172 N
SC074182 x
SC074210 N
SC074211 N
SC074212 N
SC074214 N
SC074222 N
SC074223 N
SC074227 N
SC074237 N

### 3 out of 15 are still available in object store
grep -h -f 15missing.lst *data.swarm > last3.swarm

swarm -J last3 --logdir logs --time=100:00:00 last3.swarm

ls -al bams/DCEG_Trios/rebc_RERUNS_20APR2022_0000/SC074182.pb.realn.bam  bams/DCEG_Chernobyl/0003/SC056076.pb.realn.bam bams/DCEG_Chernobyl/0004/SC074149.pb.realn.bam

### get the list of 12 missing
find bams -name "*.bam" |  pcregrep -o1 '/(SC[0-9]{6}).*.bam$' | sort -u > f3
comm -23 f1 f3 | sed -n '1!p' > 12missing.lst

grep -f 12missing.lst Trios_Master_CCAD_Biowulf_merged.clean.csv
```

### Transfer 12 bam files (aligned to hg19) from T-drive at nci-cgr to biowufl
```bash
### @nci-cgr
cat 12missing.lst | parallel 'printf "{}.bam\n{}.bam.bai\n"' > 12bams_location.txt

rsync -av --progress  --files-from=12bams_location.txt  /DCEG/Projects/Exome/SequencingData/to_archieve/BAM_Utrios helix.nih.gov:/data/DCEG_Trios/ChernobylTrios/bams/12bams_hg19
```

### Realign 12 bam files to hg38
```bash
### @biowulf
cd /data/DCEG_Trios/ChernobylTrios/

mkdir -p ref

### Transfer the reference genome hg19 used by the original bam files
### @nci-cgr
rsync -av --progress /DCEG/CGF/Bioinformatics/Production/Bari/refGenomes/Homo_sapiens_assembly19.* helix.nih.gov:/data/DCEG_Trios/ChernobylTrios/ref/

### prepare hg38 @biowulf
cp ../new_cgr_data/TriosCompass_v2/ref/Homo_sapiens_assembly38* ref/
cp /data/DCEG_Trios/refs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz* ref/

### Launch the snakemake workflow realign12bams.snakefile using the wrapper script ./run_12bams.sh
codna activate snakemake 
sbatch -J 12bams --time=200:00:00  --export=ALL --mem=12g -p norm  --wrap='./run_12bams.sh '
```

### Finally, compile a list of 341 bam files
```bash
find /data/DCEG_Trios/ChernobylTrios/bams/{set4,12bams_hg38,DCEG_Chernobyl,DCEG_Trios} -name "*.bam" > bam.lst
find /data/DCEG_Trios/new_cgr_data/TriosCompass_v2/output/fq2bam -name "SC074219.bam" >> bam.lst
```

Accordingly, we had also compiled a list of [pedigree files]() for 131 trios from the 341 samples.

---
## Configuration of the workflow

We simply put the configuration at the beginning of [the workflow file](./Snakefile_chernobyl), with the bam file locaitons extracted from the file [bam.lst](./data/bam.lst).

### Notes for dnSTR calling
We ran into ["ERROR: Not enough reads for SC056045" issue in GangSTR calling](https://github.com/gymreklab/GangSTR/i
ssues/86), and we added extra command-line argument " --insertmean 382.6 --insertsdev 135.1 " to fix the issue.

Besides, it took very long time to process two chunks: gangstr_00292.bed and gangstr_00166.bed. To speed it up, we had further split the chunks into smaller pieces: <=2 per peince.  We run all the peices using *parallel* and then merge the GangSTR predictions back.   
```bash

wc -l output/splitted_panel/gangstr_00{166,292}.bed
  1886 output/splitted_panel/gangstr_00166.bed
  1925 output/splitted_panel/gangstr_00292.bed
  3811 total

mkdir split_beds
split -l 2 --numeric-suffixes=0 --suffix-length=5 output/splitted_panel/gangstr_00292.bed --additional-suffix=".bed" split_beds/000
292_ 

split -l 2 --numeric-suffixes=0 --suffix-length=5 output/splitted_panel/gangstr_00166.bed --additional-suffix=".bed" split_beds/001
66_


mkdir gangstr_call
ls split_beds/*.bed | parallel --dry-run 'GangSTR --bam output/cram/SC109338.cram,output/cram/SC109339.cram,output/cram/SC109349.cram,output/cram/SC109351.cram,output/cram/SC109357.cram,output/cram/SC109371.cram,output/cram/SC109382.cram,output/cram/SC109384.cram,output/cram/SC109385.cram,output/cram/SC109387.cram,output/cram/SC109389.cram,output/cram/SC109391.cram,output/cram/SC109393.cram,output/cram/SC109399.cram,output/cram/SC109404.cram,output/cram/SC109407.cram,output/cram/SC109415.cram,output/cram/SC109424.cram,output/cram/SC109432.cram,output/cram/SC109441.cram,output/cram/SC109442.cram,output/cram/SC109447.cram,output/cram/SC109449.cram,output/cram/SC109451.cram,output/cram/SC109452.cram,output/cram/SC109454.cram,output/cram/SC109455.cram,output/cram/SC109456.cram,output/cram/SC109460.cram,output/cram/SC109464.cram,output/cram/SC109467.cram,output/cram/SC109505.cram,output/cram/SC109509.cram,output/cram/SC109511.cram,output/cram/SC253859.cram,output/cram/SC253860.cram,output/cram/SC253866.cram,output/cram/SC253895.cram,output/cram/SC253899.cram,output/cram/SC253900.cram,output/cram/SC260632.cram,output/cram/SC260645.cram,output/cram/SC260649.cram,output/cram/SC260653.cram,output/cram/SC260654.cram,output/cram/SC260657.cram,output/cram/SC260669.cram,output/cram/SC260671.cram,output/cram/SC260687.cram,output/cram/SC260692.cram,output/cram/SC260696.cram,output/cram/SC260702.cram,output/cram/SC260703.cram,output/cram/SC260704.cram,output/cram/SC260713.cram,output/cram/SC260716.cram,output/cram/SC260719.cram,output/cram/SC260723.cram,output/cram/SC260728.cram,output/cram/SC260734.cram,output/cram/SC499415.cram,output/cram/SC499416.cram,output/cram/SC499417.cram,output/cram/SC499419.cram,output/cram/SC499424.cram,output/cram/SC499425.cram,output/cram/SC499429.cram,output/cram/SC499431.cram,output/cram/SC499432.cram,output/cram/SC499433.cram,output/cram/SC499434.cram,output/cram/SC499435.cram,output/cram/SC499436.cram,output/cram/SC499437.cram,output/cram/SC499438.cram,output/cram/SC501089.cram,output/cram/SC501091.cram,output/cram/SC501092.cram,output/cram/SC501093.cram,output/cram/SC501095.cram,output/cram/SC501096.cram,output/cram/SC501097.cram,output/cram/SC501098.cram,output/cram/SC501099.cram,output/cram/SC501100.cram,output/cram/SC501101.cram,output/cram/SC501102.cram,output/cram/SC501104.cram,output/cram/SC501105.cram,output/cram/SC501106.cram,output/cram/SC501107.cram,output/cram/SC501108.cram,output/cram/SC501109.cram,output/cram/SC501110.cram,output/cram/SC501111.cram,output/cram/SC501112.cram,output/cram/SC501186.cram,output/cram/SC501187.cram,output/cram/SC501195.cram,output/cram/SC501198.cram,output/cram/SC501199.cram,output/cram/SC501211.cram,output/cram/SC502199.cram,output/cram/SC502200.cram,output/cram/SC502201.cram,output/cram/SC502204.cram,output/cram/SC502205.cram,output/cram/SC502206.cram,output/cram/SC502207.cram,output/cram/SC502208.cram,output/cram/SC502211.cram,output/cram/SC502214.cram,output/cram/SC502218.cram,output/cram/SC502219.cram,output/cram/SC502235.cram,output/cram/SC502236.cram,output/cram/SC502237.cram,output/cram/SC502247.cram,output/cram/SC502249.cram,output/cram/SC502250.cram,output/cram/SC502251.cram,output/cram/SC502252.cram,output/cram/SC502253.cram,output/cram/SC502259.cram,output/cram/SC502262.cram,output/cram/SC502268.cram,output/cram/SC074157.cram,output/cram/SC074166.cram,output/cram/SC074169.cram,output/cram/SC074172.cram,output/cram/SC074210.cram,output/cram/SC074211.cram,output/cram/SC074212.cram,output/cram/SC074214.cram,output/cram/SC074222.cram,output/cram/SC074223.cram,output/cram/SC074227.cram,output/cram/SC074237.cram,output/cram/SC056044.cram,output/cram/SC056045.cram,output/cram/SC056047.cram,output/cram/SC056048.cram,output/cram/SC056049.cram,output/cram/SC056050.cram,output/cram/SC056051.cram,output/cram/SC056052.cram,output/cram/SC056053.cram,output/cram/SC056055.cram,output/cram/SC056056.cram,output/cram/SC056057.cram,output/cram/SC056058.cram,output/cram/SC056059.cram,output/cram/SC056060.cram,output/cram/SC056061.cram,output/cram/SC260700.cram,output/cram/SC260706.cram,output/cram/SC260707.cram,output/cram/SC260708.cram,output/cram/SC260709.cram,output/cram/SC260710.cram,output/cram/SC260720.cram,output/cram/SC260721.cram,output/cram/SC056062.cram,output/cram/SC056064.cram,output/cram/SC056066.cram,output/cram/SC056067.cram,output/cram/SC056069.cram,output/cram/SC056072.cram,output/cram/SC056073.cram,output/cram/SC056074.cram,output/cram/SC056076.cram,output/cram/SC056077.cram,output/cram/SC056078.cram,output/cram/SC056079.cram,output/cram/SC056080.cram,output/cram/SC056081.cram,output/cram/SC056082.cram,output/cram/SC074123.cram,output/cram/SC074125.cram,output/cram/SC074126.cram,output/cram/SC074127.cram,output/cram/SC074130.cram,output/cram/SC074131.cram,output/cram/SC074132.cram,output/cram/SC074136.cram,output/cram/SC074138.cram,output/cram/SC074143.cram,output/cram/SC074144.cram,output/cram/SC074147.cram,output/cram/SC074149.cram,output/cram/SC074150.cram,output/cram/SC074151.cram,output/cram/SC074152.cram,output/cram/SC074153.cram,output/cram/SC074154.cram,output/cram/SC074128.cram,output/cram/SC074133.cram,output/cram/SC074135.cram,output/cram/SC074137.cram,output/cram/SC074141.cram,output/cram/SC074145.cram,output/cram/SC074156.cram,output/cram/SC074158.cram,output/cram/SC074159.cram,output/cram/SC074160.cram,output/cram/SC074161.cram,output/cram/SC074162.cram,output/cram/SC074163.cram,output/cram/SC074165.cram,output/cram/SC074167.cram,output/cram/SC074168.cram,output/cram/SC074170.cram,output/cram/SC074173.cram,output/cram/SC074175.cram,output/cram/SC074176.cram,output/cram/SC074177.cram,output/cram/SC074178.cram,output/cram/SC074179.cram,output/cram/SC074180.cram,output/cram/SC074181.cram,output/cram/SC074183.cram,output/cram/SC074184.cram,output/cram/SC074185.cram,output/cram/SC074186.cram,output/cram/SC074187.cram,output/cram/SC074188.cram,output/cram/SC074189.cram,output/cram/SC074190.cram,output/cram/SC074191.cram,output/cram/SC074192.cram,output/cram/SC074193.cram,output/cram/SC074194.cram,output/cram/SC074195.cram,output/cram/SC074196.cram,output/cram/SC074198.cram,output/cram/SC074200.cram,output/cram/SC074201.cram,output/cram/SC074202.cram,output/cram/SC074203.cram,output/cram/SC074204.cram,output/cram/SC074205.cram,output/cram/SC074206.cram,output/cram/SC074208.cram,output/cram/SC074209.cram,output/cram/SC074213.cram,output/cram/SC074216.cram,output/cram/SC074217.cram,output/cram/SC074218.cram,output/cram/SC074220.cram,output/cram/SC074221.cram,output/cram/SC074224.cram,output/cram/SC074225.cram,output/cram/SC074226.cram,output/cram/SC074228.cram,output/cram/SC074230.cram,output/cram/SC074231.cram,output/cram/SC253880.cram,output/cram/SC260627.cram,output/cram/SC260629.cram,output/cram/SC260630.cram,output/cram/SC260631.cram,output/cram/SC260636.cram,output/cram/SC260637.cram,output/cram/SC260639.cram,output/cram/SC260643.cram,output/cram/SC260644.cram,output/cram/SC260670.cram,output/cram/SC260682.cram,output/cram/SC260698.cram,output/cram/SC260699.cram,output/cram/SC074232.cram,output/cram/SC074233.cram,output/cram/SC074234.cram,output/cram/SC074235.cram,output/cram/SC074236.cram,output/cram/SC074238.cram,output/cram/SC074240.cram,output/cram/SC074241.cram,output/cram/SC074242.cram,output/cram/SC074243.cram,output/cram/SC074244.cram,output/cram/SC074245.cram,output/cram/SC109343.cram,output/cram/SC109347.cram,output/cram/SC109354.cram,output/cram/SC109359.cram,output/cram/SC109363.cram,output/cram/SC109375.cram,output/cram/SC109383.cram,output/cram/SC109394.cram,output/cram/SC109397.cram,output/cram/SC109398.cram,output/cram/SC109402.cram,output/cram/SC109416.cram,output/cram/SC109434.cram,output/cram/SC109440.cram,output/cram/SC109445.cram,output/cram/SC109446.cram,output/cram/SC109448.cram,output/cram/SC109453.cram,output/cram/SC109457.cram,output/cram/SC109459.cram,output/cram/SC109506.cram,output/cram/SC109510.cram,output/cram/SC109515.cram,output/cram/SC109518.cram,output/cram/SC253857.cram,output/cram/SC253862.cram,output/cram/SC253869.cram,output/cram/SC253876.cram,output/cram/SC253878.cram,output/cram/SC253882.cram,output/cram/SC253884.cram,output/cram/SC253890.cram,output/cram/SC253892.cram,output/cram/SC253897.cram,output/cram/SC253901.cram,output/cram/SC260624.cram,output/cram/SC074164.cram,output/cram/SC074171.cram,output/cram/SC074174.cram,output/cram/SC056063.cram,output/cram/SC074140.cram,output/cram/SC074146.cram,output/cram/SC074182.cram,output/cram/SC074197.cram,output/cram/SC074199.cram,output/cram/SC074207.cram,output/cram/SC074215.cram,output/cram/SC074229.cram,output/cram/SC074239.cram,output/cram/SC056054.cram,output/cram/SC056068.cram,output/cram/SC056070.cram,output/cram/SC056071.cram,output/cram/SC056075.cram,output/cram/SC056107.cram,output/cram/SC074124.cram,output/cram/SC074129.cram,output/cram/SC074139.cram,output/cram/SC074219.cram             --ref ref/Homo_sapiens_assembly38.fasta             --regions {}             --include-ggl             --bam-samps SC109338,SC109339,SC109349,SC109351,SC109357,SC109371,SC109382,SC109384,SC109385,SC109387,SC109389,SC109391,SC109393,SC109399,SC109404,SC109407,SC109415,SC109424,SC109432,SC109441,SC109442,SC109447,SC109449,SC109451,SC109452,SC109454,SC109455,SC109456,SC109460,SC109464,SC109467,SC109505,SC109509,SC109511,SC253859,SC253860,SC253866,SC253895,SC253899,SC253900,SC260632,SC260645,SC260649,SC260653,SC260654,SC260657,SC260669,SC260671,SC260687,SC260692,SC260696,SC260702,SC260703,SC260704,SC260713,SC260716,SC260719,SC260723,SC260728,SC260734,SC499415,SC499416,SC499417,SC499419,SC499424,SC499425,SC499429,SC499431,SC499432,SC499433,SC499434,SC499435,SC499436,SC499437,SC499438,SC501089,SC501091,SC501092,SC501093,SC501095,SC501096,SC501097,SC501098,SC501099,SC501100,SC501101,SC501102,SC501104,SC501105,SC501106,SC501107,SC501108,SC501109,SC501110,SC501111,SC501112,SC501186,SC501187,SC501195,SC501198,SC501199,SC501211,SC502199,SC502200,SC502201,SC502204,SC502205,SC502206,SC502207,SC502208,SC502211,SC502214,SC502218,SC502219,SC502235,SC502236,SC502237,SC502247,SC502249,SC502250,SC502251,SC502252,SC502253,SC502259,SC502262,SC502268,SC074157,SC074166,SC074169,SC074172,SC074210,SC074211,SC074212,SC074214,SC074222,SC074223,SC074227,SC074237,SC056044,SC056045,SC056047,SC056048,SC056049,SC056050,SC056051,SC056052,SC056053,SC056055,SC056056,SC056057,SC056058,SC056059,SC056060,SC056061,SC260700,SC260706,SC260707,SC260708,SC260709,SC260710,SC260720,SC260721,SC056062,SC056064,SC056066,SC056067,SC056069,SC056072,SC056073,SC056074,SC056076,SC056077,SC056078,SC056079,SC056080,SC056081,SC056082,SC074123,SC074125,SC074126,SC074127,SC074130,SC074131,SC074132,SC074136,SC074138,SC074143,SC074144,SC074147,SC074149,SC074150,SC074151,SC074152,SC074153,SC074154,SC074128,SC074133,SC074135,SC074137,SC074141,SC074145,SC074156,SC074158,SC074159,SC074160,SC074161,SC074162,SC074163,SC074165,SC074167,SC074168,SC074170,SC074173,SC074175,SC074176,SC074177,SC074178,SC074179,SC074180,SC074181,SC074183,SC074184,SC074185,SC074186,SC074187,SC074188,SC074189,SC074190,SC074191,SC074192,SC074193,SC074194,SC074195,SC074196,SC074198,SC074200,SC074201,SC074202,SC074203,SC074204,SC074205,SC074206,SC074208,SC074209,SC074213,SC074216,SC074217,SC074218,SC074220,SC074221,SC074224,SC074225,SC074226,SC074228,SC074230,SC074231,SC253880,SC260627,SC260629,SC260630,SC260631,SC260636,SC260637,SC260639,SC260643,SC260644,SC260670,SC260682,SC260698,SC260699,SC074232,SC074233,SC074234,SC074235,SC074236,SC074238,SC074240,SC074241,SC074242,SC074243,SC074244,SC074245,SC109343,SC109347,SC109354,SC109359,SC109363,SC109375,SC109383,SC109394,SC109397,SC109398,SC109402,SC109416,SC109434,SC109440,SC109445,SC109446,SC109448,SC109453,SC109457,SC109459,SC109506,SC109510,SC109515,SC109518,SC253857,SC253862,SC253869,SC253876,SC253878,SC253882,SC253884,SC253890,SC253892,SC253897,SC253901,SC260624,SC074164,SC074171,SC074174,SC056063,SC074140,SC074146,SC074182,SC074197,SC074199,SC074207,SC074215,SC074229,SC074239,SC056054,SC056068,SC056070,SC056071,SC056075,SC056107,SC074124,SC074129,SC074139,SC074219             --samp-sex M,M,F,F,F,M,F,F,M,F,M,F,M,M,F,F,M,M,F,M,M,M,M,M,M,M,F,F,M,M,F,F,M,F,F,F,M,M,F,M,F,M,F,M,M,F,F,M,F,M,F,M,M,M,M,M,M,F,M,F,M,F,M,F,M,F,F,F,M,M,M,F,F,F,F,F,M,M,M,M,F,M,M,M,M,F,F,F,F,M,M,M,F,F,M,F,F,F,M,F,M,F,M,F,M,F,F,M,M,F,F,M,F,F,M,F,M,M,F,M,M,M,F,F,F,M,F,M,M,F,F,F,F,F,F,M,M,F,F,F,M,F,F,M,M,F,F,M,M,F,F,F,F,M,M,F,M,F,F,F,M,F,M,M,F,F,M,F,F,F,M,F,F,M,M,M,M,M,F,F,M,M,M,M,M,M,M,M,F,M,F,F,F,F,F,F,F,M,F,F,F,F,M,M,M,F,M,F,M,M,F,M,F,F,F,M,F,M,F,F,F,M,F,F,M,M,F,F,F,M,F,M,F,F,M,M,M,M,M,M,M,F,M,F,M,M,F,M,M,M,M,F,F,F,M,M,M,F,M,F,M,M,F,M,F,M,F,M,M,F,M,F,M,M,F,M,M,F,M,M,M,M,F,M,M,F,F,M,F,F,M,M,M,F,F,M,M,F,F,M,M,M,M,F,F,M,M,F,F,F,M,M,M,M,M,M,F,M,M,F,M,M,F,F,F,M,F,F,M,F,M,M,F,F,F,M,F,M,M,M,F              --insertmean 382.6 --insertsdev 135.1              --out gangstr_call/{/.}' > gangstr_00292.cmd

### Merge back so as to skip it in Snakefile_chernobyl
csvtk concat -E -t gangstr_call/00166_*.tab > output/gangstr/00166.samplestats.tab
csvtk concat -E -t gangstr_call/000292_*.tab > output/gangstr/00292.samplestats.tab

```

In the future, it would be recommended to use change *split_total* from 400 to 2000, or to call GangSTR by families.   

Lastly, we got "memory overflow" issue in *DumperSTR*.  In turns out that GangSTR predicted very complicated variants in the "alt" column in certain location, which crashed *DumperSTR* in the subsequent filtering.  We had to manually removed those variants from the GangSTR output.
```bash
###  Issue found in the chunks 00000.vcf and 00175.vcf
awk '{if( length($5)>3000) print $1,$2,length($5)}' output/gangstr/00000.vcf
chr1 7862994 24803

awk '{if( length($5)>3000) print $1,$2,length($5)}' output/gangstr/00175.vcf
chr7 23231799 18878

awk '{if( length($5)>3000) print $1,$2,length($5)}' output/gangstr/{00292,00166}.vcf

### remove the two lines 
awk '{if( length($5)<3000) print $0}' output/gangstr/00000.vcf > 00000.vcf
awk '{if( length($5)<3000) print $0}' output/gangstr/00175.vcf > 00175.vcf

### overwritten the original files after double checking
wc -l output/gangstr/00000.vcf 00000.vcf output/gangstr/00175.vcf  00175.vcf
cp 00000.vcf  00175.vcf output/gangstr/

```
---

## MultiQC
```bash
module load multiqc

mkdir output_multiqc
multiqc --title QC --filename multiqc_report.html --outdir output_multiqc output/{collectmultiplemetrics,collectwgsmetrics}  --interactive

  /// MultiQC 🔍 | v1.20

|           multiqc | MultiQC Version v1.21 now available!
|           multiqc | Report title: QC
|           multiqc | Search path : /vf/users/DCEG_Trios/ChernobylTrios/TriosCompass_v2/output/collectmultiplemetrics
|           multiqc | Search path : /vf/users/DCEG_Trios/ChernobylTrios/TriosCompass_v2/output/collectwgsmetrics
|         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 8184/8184 
|            picard | Found 341 AlignmentSummaryMetrics reports
|            picard | Found 682 BaseDistributionByCycleMetrics reports
|            picard | Found 341 GcBiasMetrics reports
|            picard | Found 341 InsertSizeMetrics reports
|            picard | Found 341 QualityByCycleMetrics reports
|            picard | Found 341 QualityScoreDistributionMetrics reports
|            picard | Found 1 QualityYieldMetrics reports
|            picard | Found 341 WgsMetrics reports
|           multiqc | Report      : output_multiqc/multiqc_report.html
|           multiqc | Data        : output_multiqc/multiqc_report_data
|           multiqc | MultiQC complete

### You may transfer the files to your local computer and view the QC report
scp -r helix:/data/DCEG_Trios/ChernobylTrios/TriosCompass_v2/output_multiqc chernobyl_multiqc

open chernobyl_multiqc/multiqc_report.html
```

---
## Post-analysis
The output of workflow is identical to [the one for 107 CGR samples](https://github.com/NCI-CGR/TriosCompass_v2/blob/main/Process_107_new_Chernobyl_data.md), so is [the post-analysis part](https://github.com/NCI-CGR/TriosCompass_v2/blob/main/Process_107_new_Chernobyl_data.md#post-analyses). 