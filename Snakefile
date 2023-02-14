import os

pepfile: "pep/config.yaml"
pepschema: "workflow/schemas/manifest_schema.yaml"

samples = pep.sample_table

# print(samples)
# print(pep.config.output_dir)

# get unique value of a list
# print(list(set(samples["FLOWCELL"].tolist())))

flowcells = list(set(samples["FLOWCELL"].tolist()))
subjs = list(set(samples["SAMPLEID"].tolist()))

subj_flowcell_dict=samples[['SAMPLEID','FLOWCELL']].drop_duplicates().groupby('SAMPLEID')['FLOWCELL'].apply(list).to_dict()

ref_dict=os.path.splitext(pep.config.hg38_ref)[0]+'.dict'





def get_bam(wildcards):
    # ../bam/AH8VC6ADXX/HG002_NA24385_son_2A1_L1.bam
    df=samples.loc[samples['FLOWCELL'] == wildcards.flowcell]
    bam_files=[ "%s/%s/%s.bam" % (pep.config.input_bam_dir,row['FLOWCELL'], row['bam_name']) for index,row in df.iterrows()]
    return bam_files

def get_bam_by_subj(wildcards):
    flowcells=subj_flowcell_dict[wildcards.subj]
    return expand(pep.config.output_dir +"/gatk_markdup/{flowcell}.dedup.bam", flowcell=flowcells)

COVERAGES=[5,80]

trio_dict = {
    "trio":{
        "child": "HG002_NA24385_son",
        "father": "HG003_NA24149_father",
        "mother": "HG004_NA24143_mother" 
    }
}


TRIOS=list(trio_dict.keys())
TRIO_MEMBERS=["child", "father", "mother"]

def get_trio_hc_gvcf(wildcards):
    rv = {
        "child": pep.config.output_dir +"/gatk_hc2/"+trio_dict[wildcards.trio]["child"]+"_" + str(wildcards.cov) + ".gatk.g.vcf",
        "mother": pep.config.output_dir +"/gatk_hc2/"+trio_dict[wildcards.trio]["mother"]+"_" + str(wildcards.cov) + ".gatk.g.vcf",
        "father": pep.config.output_dir +"/gatk_hc2/"+trio_dict[wildcards.trio]["father"]+"_" + str(wildcards.cov) + ".gatk.g.vcf",
        "ref": pep.config.hg38_ref
    }
    return rv

def get_trio_bam(wildcards):
    # pep.config.output_dir+"downsample_bam/"+trio_dict[{trio}]["child"]+"_{cov}X.bam"
    rv = { key: pep.config.output_dir+"/downsample_bam/%s_%sX.bam" % (trio_dict[wildcards.trio][key], wildcards.cov) for key in TRIO_MEMBERS}

    return rv

# print(subj_flowcell_dict)

rule all:
    input: 
        expand(pep.config.output_dir +"/final_bam_metrics/{subj}.flagstat", subj=subjs),
        expand(pep.config.output_dir +"/final_bam_metrics/{subj}.stats", subj=subjs),
        expand(pep.config.output_dir +"/merge_by_subj/{subj}.merged.bam", subj=subjs),
        expand(pep.config.output_dir +"/bammetrics/{flowcell}.bammetrics.txt", flowcell=flowcells ),
        expand(pep.config.output_dir +"/collectmultiplemetrics/{flowcell}/sequencingArtifact.pre_adapter_summary_metrics.txt", flowcell=flowcells ),
        expand(pep.config.output_dir +"/downsample_bam/{subj}_{cov}X.fract", subj=subjs, cov=COVERAGES),
        expand(pep.config.output_dir +"/downsample_bam/{subj}_{cov}X.bam", subj=subjs, cov=COVERAGES),
        expand(pep.config.output_dir +"/gatk_hc2/{subj}_{cov}.gatk.g.vcf", subj=subjs, cov=COVERAGES),
        expand(pep.config.output_dir+"/gatk_cgp/{trio}_{cov}X.cgp.g.vcf.gz", trio=TRIOS , cov=COVERAGES),
        expand(pep.config.output_dir+"/deeptrio_gpu/{trio}_{cov}X.{member}.deeptrio.g.vcf.gz",trio=TRIOS , cov=COVERAGES, member=TRIO_MEMBERS)

# ./HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R1.fq.gz
rule fastp: 
    input: 
        R1= lambda w:  pep.config.input_fq_dir + "/%s" % samples.loc[w.flowcell+'_'+w.bam]["R1"],
        R2= lambda w:  pep.config.input_fq_dir + "/%s" % samples.loc[w.flowcell+'_'+w.bam]["R2"]
    output: 
        R1=pep.config.output_dir+"/fastp/{flowcell}/{bam}.R1.fastp.fastq.gz",
        R2=pep.config.output_dir+"/fastp/{flowcell}/{bam}.R2.fastp.fastq.gz",
        html=pep.config.output_dir+"/fastp/{flowcell}/{bam}.html",
        json=pep.config.output_dir+"/fastp/{flowcell}/{bam}.json"
    benchmark:
        pep.config.output_dir + "/benchmark/fastp/{flowcell}_{bam}.tsv"
    threads: 8
    resources:
        mem_mb=40000,
        runtime=600
    envmodules: "fastp"
    shell: """
        fastp \
                --in1 {input.R1} \
                --in2 {input.R2} \
                --out1 {output.R1} \
                --out2 {output.R2} \
                --report_title {wildcards.bam} \
                --json {output.json} \
                --html {output.html} \
                --thread {threads} \
                --qualified_quality_phred 25 --n_base_limit 10 --average_qual 25 \
                --length_required 50 --low_complexity_filter
    """

# bam/BHAC63ADXX/HG004_NA24143_mother_4F2_L1.bam
rule bwa:
    input: 
        R1=pep.config.output_dir+"/fastp/{flowcell}/{bam}.R1.fastp.fastq.gz",
        R2=pep.config.output_dir+"/fastp/{flowcell}/{bam}.R2.fastp.fastq.gz",
        ref=pep.config.hg38_ref
    output: pep.config.input_bam_dir + "/{flowcell}/{bam}.bam"
    params: RG=lambda w: "'@RG\\tPL:ILLUMINA\\tID:{FLOWCELL}.{LANE}\\tSM:{SAMPLEID}\\tPU:{FLOWCELL}.{LANE}.{INDEX}\\tLB:{SAMPLEID}_{INDEX}'" .format (**samples.loc[w.flowcell+'_'+w.bam])
    benchmark: 
        pep.config.output_dir + "/benchmark/bwa/{flowcell}_{bam}.tsv"
    threads: 8
    resources:
        mem_mb=40000,
        runtime=600
    envmodules: "samtools", "bwa"
    shell: """
        bwa mem -t {threads} {input.ref} {input.R1} {input.R2} -R {params.RG} | \
        samtools sort -O BAM -o {output} && \
        samtools index {output}
    """


rule merge_bam_by_flowcell: 
    input: get_bam
    output: 
        bam=pep.config.output_dir +"/merge_by_flowcell/{flowcell}.merged.bam",
        bai=pep.config.output_dir +"/merge_by_flowcell/{flowcell}.merged.bam.bai"
    threads: 16
    resources :
        mem_mb=40000,
        runtime=600
    envmodules: "samtools"
    shell: """
        samtools merge --threads {threads} {output.bam} {input}
        samtools index {output.bam}
    """

rule gatk_markdup: 
    input: pep.config.output_dir +"/merge_by_flowcell/{flowcell}.merged.bam"
    output: 
        bam=pep.config.output_dir +"/gatk_markdup/{flowcell}.dedup.bam",
        metrics=pep.config.output_dir +"/gatk_markdup/{flowcell}.dedup.metrics"
    threads: 1
    resources :
        mem_mb=40000,
        runtime=600
    envmodules: "GATK/4.2.1.0"
    shell: """
        gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics}
    """

rule bammetrics:
    input: pep.config.output_dir +"/gatk_markdup/{flowcell}.dedup.bam"
    output: pep.config.output_dir +"/bammetrics/{flowcell}.bammetrics.txt"
    params: 
        ref=pep.config.hg38_ref
    resources :
        threads=32,
        mem_mb= 80000,
        runtime=1200,
        partition="gpu",
        slurm="gres=gpu:v100x:2"
    envmodules: "parabricks/4.0.0"
    shell: """
        pbrun bammetrics \
        --ref {params.ref} \
        --bam {input} \
        --out-metrics-file {output} \
        --num-threads {threads} 
    """

rule collectmultiplemetrics:
    input: pep.config.output_dir +"/gatk_markdup/{flowcell}.dedup.bam"
    output: pep.config.output_dir +"/collectmultiplemetrics/{flowcell}/sequencingArtifact.pre_adapter_summary_metrics.txt"
    params: 
        ref=pep.config.hg38_ref
    resources :
        threads=32,
        mem_mb= 80000,
        runtime=1800,
        partition="gpu",
        slurm="gres=gpu:v100x:2"
    envmodules: "parabricks/4.0.0"
    shell: '''
        pbrun collectmultiplemetrics \
        --ref {params.ref} \
        --bam {input} \
        --out-qc-metrics-dir "$(dirname {output[0]})" \
        --gen-all-metrics 
    '''

rule merge_bam_by_subj: 
    input: get_bam_by_subj
    output: 
        bam=pep.config.output_dir +"/merge_by_subj/{subj}.merged.bam",
        bai=pep.config.output_dir +"/merge_by_subj/{subj}.merged.bam.bai"
    threads: 16
    resources :
        mem_mb=80000,
        runtime=2880
    envmodules: "samtools"
    shell: """
        samtools merge --threads {threads} {output.bam} {input}
        samtools index {output.bam}
    """

rule final_stats: 
    input: pep.config.output_dir +"/merge_by_subj/{subj}.merged.bam"
    output: pep.config.output_dir +"/final_bam_metrics/{subj}.stats"
    threads: 16
    resources :
        mem_mb=80000,
        runtime=600
    envmodules: "samtools"
    shell: """
        samtools stats -@ {threads} {input} > {output}
    """

rule final_flagstats: 
    input: pep.config.output_dir +"/merge_by_subj/{subj}.merged.bam"
    output: 
        flagstat=pep.config.output_dir +"/final_bam_metrics/{subj}.flagstat",
        idxstats=pep.config.output_dir +"/final_bam_metrics/{subj}.idxstats"
    threads: 16
    resources :
        mem_mb=80000,
        runtime=600
    envmodules: "samtools"
    shell: """
        samtools flagstat -@ {threads} {input} > {output.flagstat}
        samtools idxstat {input} > {output.idxstats}
    """

rule calc_frac:
    input: 
        idxstats=pep.config.output_dir +"/final_bam_metrics/{subj}.idxstats"
    output: 
        frac=pep.config.output_dir +"/downsample_bam/{subj}_{cov}X.fract"
    shell: """
        x=$(cut -f 3 {input} | paste -sd+ | bc); echo "scale=4;3.3*10^9 * {wildcards.cov}/148/$x" | bc > {output} 
    """

rule downsample_bam: 
    input: 
        fract=pep.config.output_dir +"/downsample_bam/{subj}_{cov}X.fract",
        bam=pep.config.output_dir +"/merge_by_subj/{subj}.merged.bam"
    output: 
        pep.config.output_dir +"/downsample_bam/{subj}_{cov}X.bam"
    threads: 16
    resources :
        mem_mb=80000,
        runtime=1200
    envmodules: "sambamba", "samtools"
    shell: """
        sambamba view -h -t {threads} -s $(cat {input.fract}) -f bam --subsampling-seed=123 {input.bam} -o {output}
        samtools index {output}
    """

rule create_dict:
    input: pep.config.hg38_ref
    output: ref_dict
    resources: 
        mem_mb=4096
    wrapper: 
        "v1.22.0/bio/picard/createsequencedictionary"

# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/gatk/haplotypecaller.html
rule gatk_hc2:
    input:
        bam=pep.config.output_dir +"/downsample_bam/{subj}_{cov}X.bam",
        ref=pep.config.hg38_ref,
        ref_dict=ref_dict
    output:
        gvcf=pep.config.output_dir +"/gatk_hc2/{subj}_{cov}.gatk.g.vcf"
    threads: 8
    benchmark:
        pep.config.output_dir +"/benchmark/gatk_hc2/{cov}/{subj}.tsv"
    resources: 
        mem_mb = 40000,
        runtime= "10d"
    wrapper:
        "v1.22.0/bio/gatk/haplotypecaller"

rule gatk_combine_gvcf:
    input: unpack(get_trio_hc_gvcf)
    output: pep.config.output_dir+"/gatk_combine_gvcf/{trio}_{cov}X.combine_gvcf.g.vcf.gz"
    threads: 16
    envmodules: "GATK/4.3.0.0"
    benchmark:
        pep.config.output_dir + "/benchmark/gatk_combine_gvcf/{trio}_{cov}X.tsv"
    resources: 
        mem_mb = 40000,
        runtime= "10d"
    shell: """
        gatk CombineGVCFs \
            -V {input.child} \
            -V {input.mother} \
            -V {input.father} \
            -R {input.ref} \
            -O {output} 
    """

rule gatk_genotype_gvcf:
    input: 
        gvcf=pep.config.output_dir+"/gatk_combine_gvcf/{trio}_{cov}X.combine_gvcf.g.vcf.gz",
        ref=pep.config.hg38_ref
    output:
        pep.config.output_dir+"/gatk_genotype_gvcf/{trio}_{cov}X.genotype_gvcf.g.vcf.gz"
    threads: 16
    envmodules: "GATK/4.3.0.0"
    benchmark:
        pep.config.output_dir + "/benchmark/gatk_genotype_gvcf/{trio}_{cov}X.tsv"
    resources: 
        mem_mb = 40000,
        runtime= "10d"
    shell: """
        gatk GenotypeGVCFs \
            -V {input.gvcf} \
            -R {input.ref} \
            -O {output} 
    """

# CalculateGenotypePosteriors
rule gatk_cgp:
    input: 
        gvcf=pep.config.output_dir+"/gatk_genotype_gvcf/{trio}_{cov}X.genotype_gvcf.g.vcf.gz",
        ped="pep/{trio}.ped"
    output: pep.config.output_dir+"/gatk_cgp/{trio}_{cov}X.cgp.g.vcf.gz"
    threads: 16
    envmodules: "GATK/4.3.0.0"
    benchmark:
        pep.config.output_dir + "/benchmark/gatk_cgp_gvcf/{trio}_{cov}X.tsv"
    resources: 
        mem_mb = 40000,
        runtime= "10d"
    shell: """
        gatk CalculateGenotypePosteriors \
            -V {input.gvcf} \
            -ped {input.ped} --skip-population-priors \
            -O {output} 
    """

rule deeptrio_gpu:
    input: 
        unpack(get_trio_bam),
        ref=pep.config.hg38_ref
    output:
        gvcf=expand(pep.config.output_dir+"/deeptrio_gpu/{{trio}}_{{cov}}X.{member}.deeptrio.g.vcf.gz", member=TRIO_MEMBERS),
        vcf=expand(pep.config.output_dir+"/deeptrio_gpu/{{trio}}_{{cov}}X.{member}.deeptrio.vcf.gz", member=TRIO_MEMBERS),
        tmpdir=directory(pep.config.output_dir+"/intermediate_results_dir_{trio}_{cov}X") 
    threads: 60
    envmodules: "singularity"
    params:
        name=lambda wildcards: [trio_dict[wildcards.trio][x] for x in TRIO_MEMBERS]
    benchmark:
        pep.config.output_dir + "/benchmark/gatk_cgp_gvcf/{trio}_{cov}X.tsv"
    resources: 
        mem_mb = 120000,
        runtime= "10d",
        partition="gpu",
        slurm="gres=gpu:v100x:4",
        tmpdir=pep.config.output_dir + "/TMP"
    shell: """
        singularity  run --pwd /mnt --nv -B /usr/lib/locale/:/usr/lib/locale/ -B "$PWD":/mnt  \
            docker://google/deepvariant:deeptrio-1.1.0-rc20201125-gpu \
            /opt/deepvariant/bin/deeptrio/run_deeptrio \
            --model_type WGS \
            -v 1 \
            --ref {input.ref} \
            --reads_child {input.child} \
            --reads_parent1 {input.father} \
            --reads_parent2 {input.mother} \
            --output_vcf_child {output.vcf[0]} \
            --output_vcf_parent1 {output.vcf[1]} \
            --output_vcf_parent2 {output.vcf[2]} \
            --sample_name_child {params.name[0]} \
            --sample_name_parent1 {params.name[1]} \
            --sample_name_parent2 {params.name[2]} \
            --num_shards {threads}  \
            --intermediate_results_dir {output.tmpdir} \
            --output_gvcf_child {output.gvcf[0]} \
            --output_gvcf_parent1 {output.gvcf[1]} \
            --output_gvcf_parent2 {output.gvcf[2]} 
    """

