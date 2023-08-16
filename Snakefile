import os
import uuid
import peds
import glob
from itertools import compress

pepfile: "pep/cgr_config.yaml"
pepschema: "schemas/cgr_manifest_schema.yaml"

samples = pep.sample_table
# samples = samples.loc[ samples['CGF_ID'].isin( ['SC108472','SC730933'])]

# print(samples)
# print(pep.config.output_dir)


# get unique value of a list
# print(list(set(samples["FLOWCELL"].tolist())))

flowcells = list(set(samples["FLOWCELL"].tolist()))
subjs = list(set(samples["CGF_ID"].tolist()))
ids = samples["sample_name"].tolist()
# print(len(ids))
# print(subjs)
# print(len(list(set(ids))))
# exit()

# print (samples)


# 'SC742296': ['BH2HYTDSX5', 'BH2JN2DSX5']
subj_flowcell_dict=samples[['CGF_ID','FLOWCELL']].drop_duplicates().groupby('CGF_ID')['FLOWCELL'].apply(list).to_dict()

ref_dict=os.path.splitext(pep.config.hg38_ref)[0]+'.dict'

# print (subj_flowcell_dict)


### Define trios
ped_files = glob.glob(pep.config.ped_dir + "/*.ped")

### make more compatible
output_dir = pep.config.output_dir

### It seems order of the family member is as the order of entries in the pedigree file
# see https://github.com/jeremymcrae/peds/blob/master/peds/ped.py
# Therefore, it is likely to be the order: father, mother and kid
# t0334   SC074198        0       0       1       1
# t0334   SC074199        0       0       2       1
# t0334   SC074201        SC074198        SC074199        1       1
families = {}
for fn in ped_files:
    f=peds.open_ped(fn)[0]
    families[f.id]=f

fam_ids = list(families.keys()) 

# print(fam_ids)


# print(subj_flowcell_dict)

rule all:
    input: 
       expand(pep.config.output_dir +"/collectmultiplemetrics/{subj}/sequencingArtifact.pre_adapter_summary_metrics.txt", subj=subjs),
       expand(pep.config.output_dir+"/fastqc/{id}.{reads}.html", id=ids, reads=['R1','R2']),
       expand(pep.config.output_dir+"/fastq_screen/{id}.fastq_screen.png", id= ids),
       expand(pep.config.output_dir +"/flagstats/{subj}.flagstats",subj=subjs),
       # expand(pep.config.output_dir +"/bammetrics/{subj}.bammetrics.txt", subj=subjs),
       expand(pep.config.output_dir + "/collectwgsmetrics/{subj}.collect_wgs_metrics.txt", subj=subjs),
       expand(pep.config.output_dir +"/call_JIGV/{caller}_{fam}.JIGV.html", fam=fam_ids, caller=["strelka", "D_and_G"]),
       expand(pep.config.output_dir +"/deepvariant_pb/SC074219.deepvariant.g.vcf.gz.tbi"),
       expand(output_dir +"/gatkhc_pb/SC074219.gatk.g.vcf.gz.tbi")
       



# ./HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R1.fq.gz
rule fastp: 
    input: 
        R1= lambda w:  samples.loc[w.sample_name]["R1"],
        R2= lambda w:  samples.loc[w.sample_name]["R2"]
    output: 
        R1=pep.config.output_dir+"/fastp/{sample_name}.R1.fastp.fastq.gz",
        R2=pep.config.output_dir+"/fastp/{sample_name}.R2.fastp.fastq.gz",
        html=pep.config.output_dir+"/fastp/{sample_name}.fastp.html",
        json=pep.config.output_dir+"/fastp/{sample_name}.fastp.json"
    benchmark:
        pep.config.output_dir + "/benchmark/fastp/{sample_name}.tsv"
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
                --report_title {wildcards.sample_name} \
                --json {output.json} \
                --html {output.html} \
                --thread {threads} \
                --qualified_quality_phred 25 --n_base_limit 10 --average_qual 25 \
                --length_required 50 --low_complexity_filter
    """

### fastqc on fastp output
# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastqc.html
# pep.config.output_dir+"/fastp/{sample_name}.R1.fastp.fastq.gz
rule fastqc:
    input: pep.config.output_dir+"/fastp/{id}.fastp.fastq.gz"
    output:
        html=pep.config.output_dir+"/fastqc/{id}.html",
        zip=pep.config.output_dir+"/fastqc/{id}_fastqc.zip" 
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{id}.log"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=600
    wrapper:
        "v2.1.1/bio/fastqc"

### fastqscreen on fastp R1 output
# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastq_screen.html
rule fastq_screen:
    input:
        pep.config.output_dir+"/fastp/{id}.R1.fastp.fastq.gz"
    output:
        txt=pep.config.output_dir+"/fastq_screen/{id}.fastq_screen.txt",
        png=pep.config.output_dir+"/fastq_screen/{id}.fastq_screen.png"
    params:
        fastq_screen_config="ref/fastq_screen.abs.conf",
        subset=100000,
        aligner='bwa'
    threads: 8
    resources:
        mem_mb=50000,
        runtime=600
    wrapper:
        "v2.1.1/bio/fastq_screen"

# https://docs.nvidia.com/clara/parabricks/4.0.1/documentation/tooldocs/man_fq2bam.html#man-fq2bam
# For the same sample, Read Groups should have the same sample name (SM) and a different ID and PU. (default: None)
rule fq2bam:
    input: 
        R1s=lambda w: expand(pep.config.output_dir+"/fastp/{id}.R1.fastp.fastq.gz", id=[ w.CGF_ID + '_'+ flowcell for flowcell in subj_flowcell_dict[w.CGF_ID] ]),
        R2s=lambda w: expand(pep.config.output_dir+"/fastp/{id}.R2.fastp.fastq.gz", id=[ w.CGF_ID + '_'+ flowcell for flowcell in subj_flowcell_dict[w.CGF_ID] ]),
        ref=pep.config.hg38_ref
    output: pep.config.output_dir + "/fq2bam/{CGF_ID}.bam"
    params: 
       fqs=lambda w, input: " --in-fq ".join([r1+' '+r2+' ' + rg for r1,r2, rg in zip(input.R1s, input.R2s, ["'@RG\\tPL:ILLUMINA\\tID:{FLOWCELL}_{LANE}\\tSM:{CGF_ID}\\tPU:{CGF_ID}_{FLOWCELL}\\tLB:{CGF_ID}_{INDEX}'" .format (**samples.loc[id]) for id in list(compress(ids, [ id == w.CGF_ID for id in samples['CGF_ID']]))  ])])
    benchmark: 
        pep.config.output_dir + "/benchmark/fq2bam/{CGF_ID}.tsv"
    threads: 48
    resources: 
        threads=48,
        mem_mb = 360000,
        disk_mb = 1000000,
        runtime= "8h",
        partition="gpu",
        slurm="gres=gpu:v100x:2",
        tmpdir=pep.config.output_dir + "/TMP"
    envmodules: "parabricks"
    shell: """
        pbrun fq2bam \
            --ref {input.ref} \
            --in-fq {params.fqs} \
            --out-bam {output} 
    """

### subj === CGF_ID
rule gatk_markdup: 
    input: pep.config.output_dir + "/fq2bam/{subj}.bam"
    output: 
        bam=pep.config.output_dir +"/gatk_markdup/{subj}.dedup.bam",
        metrics=pep.config.output_dir +"/gatk_markdup/{subj}.dedup.metrics"
    threads: 4
    resources :
        mem_mb=80000,
        runtime="24h"
    envmodules: "GATK/4.2.1.0", "samtools"
    shell: """
        gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics}
        samtools index -@ {threads} {output.bam}
    """



rule flagstats: 
    input: pep.config.output_dir +"/gatk_markdup/{subj}.dedup.bam"
    output: 
        flagstat=pep.config.output_dir +"/flagstats/{subj}.flagstats",
        idxstats=pep.config.output_dir +"/flagstats/{subj}.idxstats"
    threads: 16
    resources :
        mem_mb=80000,
        runtime=600
    envmodules: "samtools"
    shell: """
        samtools flagstat -@ {threads} {input} > {output.flagstat}
        samtools idxstat {input} > {output.idxstats}
    """



### Collect metrics is slow and repalced with collectwgsmetrics
# rule bammetrics:
#     input: pep.config.output_dir +"/fq2bam/{subj}.bam"
#     output: pep.config.output_dir +"/bammetrics/{subj}.bammetrics.txt"
#     params: 
#         ref=pep.config.hg38_ref
#     resources :
#         threads=32,
#         mem_mb= 120000,
#         disk_mb = 1000000,
#         runtime="3d",
#         partition="gpu",
#         slurm="gres=gpu:v100x:2"
#     benchmark:
#         pep.config.output_dir +"/benchmark/bammetrics/{subj}.tsv"
#     envmodules: "parabricks"
#     shell: """
#         pbrun bammetrics \
#         --ref {params.ref} \
#         --bam {input} \
#         --out-metrics-file {output} \
#         --num-threads {threads} 
#     """

rule collectwgsmetrics:
    input: 
        bam= pep.config.output_dir +"/gatk_markdup/{subj}.dedup.bam",
        ref=pep.config.hg38_ref
    output: pep.config.output_dir + "/collectwgsmetrics/{subj}.collect_wgs_metrics.txt"
    resources :
        mem_mb= 80000,
        runtime="1d"
    benchmark:
        pep.config.output_dir +"/benchmark/collectwgsmetrics/{subj}.tsv"
    envmodules: "picard/2.27.3"
    shell: """
        java -jar $PICARDJAR CollectWgsMetrics \
            I={input.bam} \
            O={output} \
            R={input.ref} 
    """

# qualityscore.png and qualityscore.pdf are redundant

rule collectmultiplemetrics:
    input: pep.config.output_dir +"/gatk_markdup/{subj}.dedup.bam"
    output: pep.config.output_dir +"/collectmultiplemetrics/{subj}/sequencingArtifact.pre_adapter_summary_metrics.txt"
    benchmark:
        pep.config.output_dir +"/benchmark/collectmultiplemetrics/{subj}.tsv"
    params: 
        ref=pep.config.hg38_ref,
        prefix= pep.config.output_dir +"/collectmultiplemetrics/{subj}"
    resources :
        threads=32,
        mem_mb= 80000,
        runtime="6h",
        partition="gpu",
        slurm="gres=gpu:v100x:2"
    envmodules: "parabricks"
    shell: '''
        pbrun collectmultiplemetrics \
        --ref {params.ref} \
        --bam {input} \
        --out-qc-metrics-dir {params.prefix} \
        --gen-all-metrics 
    '''

rule gatkhc_pb:
    input:
        bam=pep.config.output_dir +"/gatk_markdup/{subj}.dedup.bam",
        ref=pep.config.hg38_ref,
        ref_dict=ref_dict
    output:
        gvcf=output_dir +"/gatkhc_pb/{subj}.gatk.g.vcf",
        gz=output_dir +"/gatkhc_pb/{subj}.gatk.g.vcf.gz",
        tbi=output_dir +"/gatkhc_pb/{subj}.gatk.g.vcf.gz.tbi"
    threads: 24
    benchmark:
        pep.config.output_dir +"/benchmark/gatkhc_pb/{subj}.tsv"
    resources: 
        mem_mb = 240000,
        runtime= "1d",
        partition="gpu",
        slurm="gres=gpu:v100x:2",
        tmpdir=pep.config.output_dir + "/TMP"
    envmodules: "parabricks", "bcftools"
    shell: """
        pbrun haplotypecaller \
            --gvcf \
            --ref {input.ref} \
            --in-bam {input.bam} \
            --out-variants {output.gvcf}   
        bgzip -c {output.gvcf} > {output.gz}
        tabix {output.gz}     
    """
 
### Call DV
rule deepvariant_pb:
    input:
        bam=pep.config.output_dir +"/gatk_markdup/{subj}.dedup.bam",
        ref=pep.config.hg38_ref
    output:
        gvcf=pep.config.output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf",
        gz=pep.config.output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf.gz",
        tbi=pep.config.output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf.gz.tbi"
    threads: 48
    benchmark:
        pep.config.output_dir +"/benchmark/deepvariant_pb/{subj}.tsv"
    envmodules: "parabricks", "bcftools"
    resources: 
        mem_mb = 180000,
        runtime= "1d",
        partition="gpu",
        slurm="gres=gpu:v100x:2",
        tmpdir=output_dir + "/TMP"
    shell: """
        pbrun deepvariant \
        --gvcf  \
        --ref {input.ref} \
        --in-bam {input.bam} \
        --out-variants {output.gvcf} 
        bgzip -c {output.gvcf} > {output.gz}
        tabix {output.gz}
    """

### as RG tag is not properly set in bam
# rule fix_dv_vcf:
#     input: output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf"
#     output: 
#         vcf=output_dir +"/deepvariant_pb/{subj}.deepvariant.fixed.vcf.gz",
#         tbi=output_dir +"/deepvariant_pb/{subj}.deepvariant.fixed.vcf.gz.tbi"
#     # params: id=lambda w: w.id
#     envmodules: "bcftools"
#     threads: 8
#     resources: 
#         mem_mb = 180000,
#         runtime= "1d"
#     shell: """
#         bcftools reheader -s  <(echo -e "{wildcards.subj}")  {input} -T {threads} | bgzip > {output.vcf}
#         tabix {output.vcf}
#     """


# use rule fix_dv_vcf as fix_gatk_vcf with:
#     input: 
#         output_dir +"/gatkhc_pb/{id}.gatk.g.vcf"
#     output: 
#         vcf=output_dir +"/gatkhc_pb/{id}.gatk.fixed.vcf.gz",
#         tbi=output_dir +"/gatkhc_pb/{id}.gatk.fixed.vcf.gz.tbi"

### vcf index is required here
rule gatk_combine_gvcf:
    input: 
        vcfs=lambda w: expand(output_dir + "/gatkhc_pb/{id}.gatk.g.vcf.gz", id= [person.id for person in families[w.fam]]),
        tbis=lambda w: expand(output_dir + "/gatkhc_pb/{id}.gatk.g.vcf.gz.tbi", id= [person.id for person in families[w.fam]]),
        ref=pep.config.hg38_ref
    output: output_dir+"/gatk_combine_gvcf/{fam}.combine_gvcf.g.vcf.gz"
    threads: 16
    envmodules: "GATK/4.3.0.0"
    benchmark:
        output_dir + "/benchmark/gatk_combine_gvcf/{fam}.tsv"
    resources: 
        mem_mb = 40000,
        runtime= "10d"
    params: 
        v=lambda w, input: " -V ".join(input.vcfs)
    shell: """
        gatk CombineGVCFs \
            -V {params.v} \
            -R {input.ref} \
            -O {output} 
    """

### 
rule gatk_genotype_gvcf_pb:
    input: 
        gvcf=output_dir+"/gatk_combine_gvcf/{fam}.combine_gvcf.g.vcf.gz",
        ref=pep.config.hg38_ref
    output:
        vcf=output_dir+"/gatk_genotype_gvcf_pb/{fam}.genotype_gvcf.g.vcf",
        gz=output_dir+"/gatk_genotype_gvcf_pb/{fam}.genotype_gvcf.g.vcf.gz",
    threads: 16
    envmodules: "parabricks", "bcftools"
    benchmark:
        output_dir + "/benchmark/gatk_genotype_gvcf_pb/{fam}.tsv"
    resources: 
        mem_mb = 80000,
        runtime= "1d",
        partition="gpu",
        slurm="gres=gpu:v100x:2",
        tmpdir=output_dir + "/TMP"
    shell: """
        pbrun genotypegvcf \
            --in-gvcf {input.gvcf} \
            --ref {input.ref} \
            --out-vcf {output.vcf} 
        bgzip -@ {threads} -c {output.vcf}> {output.gz} && \
        tabix {output.gz}
    """

rule gatk_cgp:
    input: 
        gvcf=output_dir+"/gatk_genotype_gvcf_pb/{fam}.genotype_gvcf.g.vcf.gz",
        ped=pep.config.ped_dir+"/{fam}.ped",
        ref=pep.config.hg38_ref
    output: 
        gvcf=output_dir+"/gatk_cgp/{fam}.cgp.g.vcf.gz",
        norm=output_dir+"/gatk_cgp/{fam}.cgp_norm.vcf.gz",
    threads: 16
    envmodules: "GATK/4.3.0.0", "bcftools"
    benchmark:
        output_dir + "/benchmark/gatk_cgp/{fam}.tsv"
    resources: 
        mem_mb = 40000,
        runtime= "10d"
    shell: """
        gatk CalculateGenotypePosteriors \
            -V {input.gvcf} \
            -ped {input.ped} --skip-population-priors \
            -O {output.gvcf} 
        bcftools norm -f {input.ref} -m -  {output.gvcf} | bcftools view -i'ALT!="*"' -Oz -o {output.norm}
    """

rule glnexus_dv: 
    input:
        gvcf=lambda w: expand(output_dir +"/deepvariant_pb/{id}.deepvariant.g.vcf.gz", id=[person.id for person in families[w.fam]]),
        tbi=lambda w: expand(output_dir +"/deepvariant_pb/{id}.deepvariant.g.vcf.gz.tbi", id=[person.id for person in families[w.fam]]),
        ref=pep.config.hg38_ref
    output:
        vcf = output_dir +"/glnexus/{fam}.dv_combined.vcf.gz"
    threads: 8
    benchmark:
        output_dir +"/benchmark/glnexus_dv/dv_{fam}.tsv"
    resources: 
        mem_mb = 100*1000,
        runtime= "1d"
    envmodules: "glnexus", "bcftools"
    params: 
        tempdir=temp(directory(output_dir +"/glnexus/GLnexus.DB_{}")).format(uuid.uuid4()) 
    shell: """
        glnexus_cli -t {threads} --config DeepVariantWGS --dir {params.tempdir} {input.gvcf} |  bcftools norm -f {input.ref} -m - -O z -o {output} 
        rm -fr {params.tempdir}
    """

### fix DP (for strelka) and AD
rule call_dnm_dv: 
    input: 
        vcf = output_dir +"/glnexus/{fam}.dv_combined.vcf.gz",
        ped=pep.config.ped_dir+"/{fam}.ped",
        interval="ref/hg38.wgs_interval.bed"
    output:
        vcf= output_dir +"/slivar/DV_{fam}.dnm.vcf",
        norm_vcf= temp(output_dir +"/slivar/DV_{fam}.tmp0.vcf.gz"),
        tmp = temp(output_dir +"/slivar/DV_{fam}.tmp.vcf.gz"),
        gz = output_dir +"/slivar/DV_{fam}.dnm.vcf.gz"
    benchmark:
        output_dir +"/benchmark/slivar/DV_{fam}.tsv"
    resources: 
        mem_mb = 20*1000,
        runtime= "1d"
    params: min_gq=10, min_dp=20
    envmodules: "bcftools"
    conda:
        "workflow/envs/slivar.yaml"
    shell: """
        zcat {input.vcf} | sed -e 's/ID=AD,Number=\./ID=AD,Number=R/' |bcftools norm -f ref/Homo_sapiens_assembly38.fasta -m - -O z -o {output.norm_vcf} 
        tabix {output.norm_vcf}
        slivar expr  \
            --vcf {output.norm_vcf} \
            --ped  {input.ped} \
            --pass-only \
            --out-vcf {output.vcf} \
            --trio "denovo:( \
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
        bgzip -c {output.vcf} > {output.tmp}
        tabix {output.tmp}
        bcftools view -R {input.interval} {output.tmp} -O z -o {output.gz}
        tabix {output.gz}
    """


# use rule call_dnm_dv as call_dnm_dv2 with: 
#     output:
#         vcf= output_dir +"/slivar/DV2_{fam}.dnm.vcf",
#         norm_vcf= temp(output_dir +"/slivar/DV2_{fam}.tmp0.vcf.gz"),
#         tmp = temp(output_dir +"/slivar/DV2_{fam}.tmp.vcf.gz"),
#         gz = output_dir +"/slivar/DV2_{fam}.dnm.vcf.gz"
#     benchmark:
#         output_dir +"/benchmark/slivar/DV2_{fam}.tsv"
#     params: min_gq=10, min_dp=20


### Call dnm from gatk
use rule  call_dnm_dv as call_dnm_gatk with:
    input: 
        vcf=output_dir+"/gatk_cgp/{fam}.cgp_norm.vcf.gz",
        ped=pep.config.ped_dir + "/{fam}.ped",
        interval="ref/hg38.wgs_interval.bed"
    output: 
       vcf= output_dir +"/slivar/GATK_{fam}.dnm.vcf",
       norm_vcf= temp(output_dir +"/slivar/GATK_{fam}.tmp0.vcf.gz"),
       tmp = temp(output_dir +"/slivar/GATK_{fam}.tmp.vcf.gz"),
       gz = output_dir +"/slivar/GATK_{fam}.dnm.vcf.gz"
    params: 
       min_gq=20, 
       min_dp=30
    benchmark:
        output_dir +"/benchmark/slivar/GATK_{fam}.tsv"

rule merge_DV_GATK:
    input: 
        DV=output_dir +"/slivar/DV_{fam}.dnm.vcf.gz",
        GATK=output_dir +"/slivar/GATK_{fam}.dnm.vcf.gz"
    output:
        gz=output_dir +"/GATK_DV/{fam}.merge.dnm.vcf.gz",
        tbi=output_dir +"/GATK_DV/{fam}.merge.dnm.vcf.gz.tbi",
        both=output_dir +"/GATK_DV/D_and_G.{fam}.dnm.vcf.gz",
        one=output_dir +"/GATK_DV/D_or_G.{fam}.dnm.vcf.gz"
    envmodules: "bcftools"
    shell: """
        bcftools merge --force-samples --threads 2 -m none {input.DV} {input.GATK} -Oz -o {output.gz}
        tabix {output.gz}
        bcftools filter -i "N_MISSING == 0" -Oz -o {output.both} {output.gz}
        bcftools filter -i "N_MISSING > 0" -Oz -o {output.one} {output.gz}
        tabix {output.both}
        tabix {output.one}
    """



### prepare the bed regions cover dnm calls from both GATK and DV for each trio
rule get_dnm_regions:
    input: 
        output_dir +"/GATK_DV/D_or_G.{fam}.dnm.vcf.gz",
    output:
        bed=output_dir +"/strelka/call_regioins.{fam}.bed.gz"
    envmodules: "bedtools", "bcftools"
    shell: """
        mergeBed -i {input} | bgzip -c - > {output.bed}
        tabix {output.bed}
    """


### Use call_regioins2 and DV2_{fam}.dnm.vcf.gz instead
rule config_strelka:
    input: 
        bams= lambda w: expand(pep.config.output_dir +"/gatk_markdup/{subj}.dedup.bam", subj= [person.id for person in families[w.fam]]),
        bed= output_dir +"/strelka/call_regioins.{fam}.bed.gz",
        ref=pep.config.hg38_ref,
        vcf=output_dir +"/GATK_DV/D_or_G.{fam}.dnm.vcf.gz"
    output:
        conf=output_dir + "/strelka/{fam}/runWorkflow.py"
    envmodules: "strelka"
    params: 
        bams=lambda w, input: " --bam ".join(input.bams),
        runDir=output_dir + "/strelka/{fam}"
    resources: mem_mb=10000 
    shell: """
        configureStrelkaGermlineWorkflow.py \
        --forcedGT {input.vcf} \
        --bam {params.bams} \
        --referenceFasta {input.ref} \
        --runDir {params.runDir} \
        --callRegions {input.bed}
    """

rule run_strelka:
    input: cmd=output_dir + "/strelka/{fam}/runWorkflow.py"
    output:
        vcf=output_dir + "/strelka/{fam}/results/variants/variants.vcf.gz"
    envmodules: "strelka"
    threads: 16
    benchmark:
        output_dir +"/benchmark/run_strelka/strelka_{fam}.tsv"
    resources: 
        mem_mb = 100*1000,
        runtime= "1d"
    shell: """
        {input.cmd} -m local -j {threads} 
    """

use rule  call_dnm_dv as call_dnm_strelka with:
    input: 
        vcf=output_dir+"/strelka/{fam}/results/variants/variants.vcf.gz",
        ped=pep.config.ped_dir + "/{fam}.ped",
        interval="ref/hg38.wgs_interval.bed"
    output: 
       vcf= output_dir +"/slivar/strelka_{fam}.dnm.vcf",
       norm_vcf= temp(output_dir +"/slivar/strelka_{fam}.tmp0.vcf.gz"),
       tmp = temp(output_dir +"/slivar/strelka_{fam}.tmp.vcf.gz"),
       gz = output_dir +"/slivar/strelka_{fam}.dnm.vcf.gz"
    params: 
       min_gq=20, 
       min_dp=30
    benchmark:
        output_dir +"/benchmark/slivar/strelka_{fam}.tsv"

### Call JIGV
rule call_JIGV:
    input:
        bam=lambda w: expand(pep.config.output_dir +"/gatk_markdup/{subj}.dedup.bam", subj=[person.id for person in families[w.fam]]),
        ped=pep.config.ped_dir + "/{fam}.ped",
        ref=pep.config.hg38_ref,
        sites= output_dir +"/slivar/{caller}_{fam}.dnm.vcf.gz"
    output:
        html= output_dir +"/call_JIGV/{caller}_{fam}.JIGV.html"
    benchmark:
        output_dir +"/benchmark/call_JIGV/{caller}_{fam}.tsv"
    resources: 
        mem_mb = 60*1000,
        runtime= "3d"
    params: 
        proband=lambda w: [person.id for person in families[w.fam] if families[w.fam].get_father(person) ][0]
    shell: """
        jigv  \
            --fasta {input.ref} \
            --sample {params.proband} \
            --ped {input.ped} \
            --sites {input.sites} \
            {input.bam} > {output.html}
    """

### Call JIGV for DV&GATK
use rule call_JIGV as call_JIGV_for_DG with:
    input:
        bam=lambda w: expand(pep.config.output_dir +"/gatk_markdup/{subj}.dedup.bam", subj=[person.id for person in families[w.fam]]),
        ped=pep.config.ped_dir + "/{fam}.ped",
        ref=pep.config.hg38_ref,
        sites= output_dir +"/GATK_DV/{caller}.{fam}.dnm.vcf.gz"
  
# # Escape certain characters, such as \t by \\t, $ by \$, and { by {{.
# rule true_sites_bed:
#     input: "ChernobylTriosTruth/all_clean/{fam}c1.csv"
#     output:
#         output_dir + "true_sites_bed/truth_{fam}c1.bed"
#     shell: """
#          awk -v FS=',' -v OFS='\t' "{{ if(NR>1) print \$1,\$2-1,\$2}}" {input} | sort -u > {output}
#     """


### Call JIGV for true sites
# use rule call_JIGV as call_JIGV_for_truth with:
#     input:
#         bam=lambda w: expand(output_dir +"/fixed-rg/{id}.bam", id=[person.id for person in families[w.fam]]),
#         ped="ped_files/{fam}.ped",
#         ref=hg38_ref,
#         sites= output_dir + "true_sites_bed/{caller}_{fam}c1.bed"
  
### Call JIGV for curation
# true_sites_not_predicted
# use rule call_JIGV as call_JIGV_for_true_sites_not_predicted with:
#     input:
#         bam=lambda w: expand(output_dir +"/fixed-rg/{id}.bam", id=[person.id for person in families[w.fam]]),
#         ped="ped_files/{fam}.ped",
#         ref=hg38_ref,
#         sites= "curation/{caller}.{fam}.bed"
#     output:
#         html= "curation/JIGV/{caller}.{fam}.JIGV.html"

# use rule call_JIGV as call_JIGV_for_new_candidates with:
#     input:
#         bam=lambda w: expand(output_dir +"/fixed-rg/{id}.bam", id=[person.id for person in families[w.fam]]),
#         ped="ped_files/{fam}.ped",
#         ref=hg38_ref,
#         sites= "curation/{caller}.{fam}.vcf.gz"
#     output:
#         html= "curation/JIGV/{caller}.{fam}.JIGV.html"



