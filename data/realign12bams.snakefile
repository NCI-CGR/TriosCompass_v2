import glob

HG19="ref/Homo_sapiens_assembly19.fasta"
HG38="ref/Homo_sapiens_assembly38.fasta"
BAM_HG19_DIR="bams/12bams_hg19/"
BAM_HG38_DIR="bams/12bams_hg38/"
OUTPUT_DIR="output/"
knownSites="ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# bams/12bams_hg19/SC074157.bam  bams/12bams_hg19/SC074210.bam  bams/12bams_hg19/SC074222.bam
# bams/12bams_hg19/SC074166.bam  bams/12bams_hg19/SC074211.bam  bams/12bams_hg19/SC074223.bam
# bams/12bams_hg19/SC074169.bam  bams/12bams_hg19/SC074212.bam  bams/12bams_hg19/SC074227.bam
# bams/12bams_hg19/SC074172.bam  bams/12bams_hg19/SC074214.bam  bams/12bams_hg19/SC074237.bam
bam_files = glob.glob(BAM_HG19_DIR + "**/*.bam", recursive=True)

### , is required and we will get Wildcards(id=['SC074157', 'SC074166', 'SC074169', 'SC074172', 'SC074210', 'SC074211', 'SC074212', 'SC074214', 'SC074222', 'SC074223', 'SC074227', 'SC074237'] (which is not a list as espected) otherwise.
ids, = glob_wildcards(BAM_HG19_DIR + "{id,SC[0-9]+}.bam", bam_files)
# print (ids)
# exit(0)

rule all:
    input: 
        # expand(OUTPUT_DIR + "bam2fq/{id}_{reads}.fastq.gz", id=ids, reads=["R1", "R2"])
        expand(BAM_HG38_DIR + "{id}.bam", id=ids)

# prefix + 1.fastq.gz
rule bam2fq:
    input: 
        bam = BAM_HG19_DIR + "{id}.bam",
        ref = HG19
    output: 
        R1= OUTPUT_DIR + "bam2fq/{id}_1.fastq.gz",
        R2= OUTPUT_DIR + "bam2fq/{id}_2.fastq.gz"
    params:
        prefix = OUTPUT_DIR + "bam2fq/{id}"
    benchmark: 
        OUTPUT_DIR  + "benchmark/bam2fq/{id}.tsv"
    threads: 48
    resources: 
        threads=48,
        mem_mb = 360000,
        disk_mb = 1000000,
        runtime= "8h",
        partition="gpu",
        slurm="gres=gpu:v100x:2",
        tmpdir=OUTPUT_DIR  + "TMP"
    envmodules: "parabricks"
    shell: """
        pbrun bam2fq \
            --in-bam {input.bam} \
            --ref {input.ref} \
            --out-prefix {params.prefix} 
    """

# "@RGtID:footLB:lib1tPL:bartSM:sampletPU:unit1"
# "'@RG\\tPL:ILLUMINA\\tID:{id}\\tSM:{id}\\tPU:{id}\\tLB:LIB1'"
# Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
rule fq2bam:
    input: 
        R1= OUTPUT_DIR + "bam2fq/{id}_1.fastq.gz",
        R2= OUTPUT_DIR + "bam2fq/{id}_2.fastq.gz",
        ref=HG38,
        knownSites=knownSites
    output: 
        bam=BAM_HG38_DIR + "{id}.bam",
        bqsr=BAM_HG38_DIR + "{id}.BQSR-REPORT.txt"
    params: 
        RG="'@RG\\tPL:ILLUMINA\\tID:{id}\\tSM:{id}\\tPU:{id}\\tLB:LIB1'"
    benchmark: 
        OUTPUT_DIR + "benchmark/fq2bam/{id}.tsv"
    threads: 48
    resources: 
        threads=48,
        mem_mb = 360000,
        disk_mb = 1000000,
        runtime= "8h",
        partition="gpu",
        slurm="gres=gpu:v100x:2",
        tmpdir=OUTPUT_DIR + "/TMP"
    envmodules: "parabricks"
    shell: """
        pbrun fq2bam \
            --ref {input.ref} \
            --in-fq {input.R1} {input.R2} {params.RG} \
            --out-bam {output.bam} \
            --knownSites {input.knownSites} \
            --out-recal-file {output.bqsr}          
    """