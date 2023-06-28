
import glob
import peds 
import uuid

# bam/SC501095.recal.bam
bam_files = glob.glob("bam/**/*.bam", recursive=True)
_,ids,_= glob_wildcards("{path}/{id,SC[0-9]+}.{any,.*}bam", bam_files)


BAM_DICT=dict(zip(ids,bam_files))
# print(bam_files)
# print(BAM_DICT)
# exit(0)

### Define trios
ped_files = glob.glob("ped_files/*.ped")

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

# fam_ids = [f.id for f in famimiles]
# 8trios.lst
fam_ids = [line.strip() for line in open("8trios.lst", "r")]

callers=['DV', 'GATK', 'strelka']

# print(fam_ids)
# print(BAM_DICT['SC501095'])
###

### configure
output_dir="output"
wgs_interval= "ref/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list"
hg38_ref="ref/Homo_sapiens_assembly38.fasta"
ref_dict=os.path.splitext(hg38_ref)[0]+'.dict'


### 

rule all: 
    input: 
        expand(output_dir +"/collectmultiplemetrics/{id}/sequencingArtifact.pre_adapter_summary_metrics.txt", id=ids),
        expand(output_dir + "/collectwgsmetrics/{id}.collect_wgs_metrics.txt", id=ids),
        expand(output_dir +"/gatkhc_pb/{id}.gatk.g.vcf", id=ids),
        expand(output_dir +"/deepvariant_pb/{id}.deepvariant.g.vcf", id=ids),
        expand(output_dir +"/deepvariant_pb/{id}.deepvariant.fixed.vcf.gz", id=ids),
        expand(output_dir +"/gatkhc_pb/{id}.gatk.fixed.vcf.gz", id=ids),
        expand(output_dir+"/gatk_cgp/{fam}.cgp.g.vcf.gz", fam=fam_ids),
        expand(output_dir +"/glnexus/{fam}.dv_combined.vcf.gz", fam=fam_ids),
        expand(output_dir +"/slivar/{caller}_{fam}.dnm.vcf.gz", fam=fam_ids, caller=callers),
        expand(output_dir + "/fixed-rg/{id}.bam", id=ids),
        expand(output_dir +"/call_JIGV/{caller}_{fam}.JIGV.html", fam=fam_ids, caller=["strelka", "D_and_G","truth"])

rule replace_rg:
    input: lambda w: BAM_DICT[w.id]
    output:
        output_dir + "/fixed-rg/{id}.bam"
    params:
        extra="--RGLB lib1 --RGPL illumina --RGPU {id} --RGSM {id} --CREATE_INDEX true "
    benchmark:
        output_dir +"/benchmark/replace_rg/{id}.tsv"
    resources:
        mem_mb=60000,
        runtime='5d',
        threads=8,
        tmpdir=output_dir+"/TMP"
    wrapper:
        "v1.25.0/bio/picard/addorreplacereadgroups"

rule collectmultiplemetrics:
    input: lambda w: BAM_DICT[w.id]
    output: output_dir +"/collectmultiplemetrics/{id}/sequencingArtifact.pre_adapter_summary_metrics.txt"
    params: 
        ref=hg38_ref
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

# picard/gatk CollectWGSMetrics is equivalent to bammetrics
# https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard-
# https://hpc.nih.gov/apps/picard.html
rule collectwgsmetrics:
    input: 
        bam=lambda w: BAM_DICT[w.id],
        ref=hg38_ref
    output: output_dir + "/collectwgsmetrics/{id}.collect_wgs_metrics.txt"
    resources :
        mem_mb= 80000,
        runtime="3d"
    benchmark:
        output_dir +"/benchmark/collectwgsmetrics/{id}.tsv"
    envmodules: "picard/2.27.3"
    shell: """
        java -jar $PICARDJAR CollectWgsMetrics \
            I={input.bam} \
            O={output} \
            R={input.ref} 
    """

rule gatkhc_pb:
    input:
        bam=lambda w: BAM_DICT[w.id],
        ref=hg38_ref,
        ref_dict=ref_dict
    output:
        gvcf=output_dir +"/gatkhc_pb/{id}.gatk.g.vcf"
    threads: 24
    benchmark:
        output_dir +"/benchmark/gatkhc_pb/{id}.tsv"
    resources: 
        mem_mb = 180000,
        runtime= "1d",
        partition="gpu",
        slurm="gres=gpu:v100x:4",
        tmpdir=output_dir + "/TMP"
    envmodules: "parabricks/4.0.0"
    shell: """
        pbrun haplotypecaller \
            --gvcf \
            --ref {input.ref} \
            --in-bam {input.bam} \
            --out-variants {output.gvcf}        
    """
 
### Call DV
rule deepvariant_pb:
    input:
        bam=lambda w: BAM_DICT[w.id],
        ref=hg38_ref
    output:
        gvcf=output_dir +"/deepvariant_pb/{id}.deepvariant.g.vcf"
    threads: 48
    benchmark:
        output_dir +"/benchmark/deepvariant_pb/{id}.tsv"
    envmodules: "parabricks/4.0.0", "bcftools"
    resources: 
        mem_mb = 180000,
        runtime= "1d",
        partition="gpu",
        slurm="gres=gpu:v100x:4",
        tmpdir=output_dir + "/TMP"
    shell: """
        pbrun deepvariant \
        --gvcf  \
        --ref {input.ref} \
        --in-bam {input.bam} \
        --out-variants {output.gvcf} 
    """

### as RG tag is not properly set in bam
rule fix_dv_vcf:
    input: output_dir +"/deepvariant_pb/{id}.deepvariant.g.vcf"
    output: 
        vcf=output_dir +"/deepvariant_pb/{id}.deepvariant.fixed.vcf.gz",
        tbi=output_dir +"/deepvariant_pb/{id}.deepvariant.fixed.vcf.gz.tbi"
    # params: id=lambda w: w.id
    envmodules: "bcftools"
    threads: 8
    resources: 
        mem_mb = 180000,
        runtime= "1d"
    shell: """
        bcftools reheader -s  <(echo -e "{wildcards.id}")  {input} -T {threads} | bgzip > {output.vcf}
        tabix {output.vcf}
    """


use rule fix_dv_vcf as fix_gatk_vcf with:
    input: 
        output_dir +"/gatkhc_pb/{id}.gatk.g.vcf"
    output: 
        vcf=output_dir +"/gatkhc_pb/{id}.gatk.fixed.vcf.gz",
        tbi=output_dir +"/gatkhc_pb/{id}.gatk.fixed.vcf.gz.tbi"

### vcf index is required here
rule gatk_combine_gvcf:
    input: 
        vcfs=lambda w: expand(output_dir + "/gatkhc_pb/{id}.gatk.fixed.vcf.gz", id= [person.id for person in families[w.fam]]),
        tbis=lambda w: expand(output_dir + "/gatkhc_pb/{id}.gatk.fixed.vcf.gz.tbi", id= [person.id for person in families[w.fam]]),
        ref=hg38_ref
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
        ref=hg38_ref
    output:
        vcf=output_dir+"/gatk_genotype_gvcf_pb/{fam}.genotype_gvcf.g.vcf",
        gz=output_dir+"/gatk_genotype_gvcf_pb/{fam}.genotype_gvcf.g.vcf.gz",
    threads: 16
    envmodules: "parabricks/4.0.0", "bcftools"
    benchmark:
        output_dir + "/benchmark/gatk_genotype_gvcf_pb/{fam}.tsv"
    resources: 
        mem_mb = 80000,
        runtime= "1d",
        partition="gpu",
        slurm="gres=gpu:v100x:4",
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
        ped="ped_files/{fam}.ped",
        ref=hg38_ref
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
        gvcf=lambda w: expand(output_dir +"/deepvariant_pb/{id}.deepvariant.fixed.vcf.gz", id=[person.id for person in families[w.fam]]),
        tbi=lambda w: expand(output_dir +"/deepvariant_pb/{id}.deepvariant.fixed.vcf.gz.tbi", id=[person.id for person in families[w.fam]]),
        ref=hg38_ref
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
        ped="ped_files/{fam}.ped",
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
        ped="ped_files/{fam}.ped",
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
        bams= lambda w: expand(output_dir + "/fixed-rg/{id}.bam", id= [person.id for person in families[w.fam]]),
        bed= output_dir +"/strelka/call_regioins.{fam}.bed.gz",
        ref=hg38_ref,
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
        ped="ped_files/{fam}.ped",
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
        bam=lambda w: expand(output_dir +"/fixed-rg/{id}.bam", id=[person.id for person in families[w.fam]]),
        ped="ped_files/{fam}.ped",
        ref=hg38_ref,
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
        bam=lambda w: expand(output_dir +"/fixed-rg/{id}.bam", id=[person.id for person in families[w.fam]]),
        ped="ped_files/{fam}.ped",
        ref=hg38_ref,
        sites= output_dir +"/GATK_DV/{caller}.{fam}.dnm.vcf.gz"
  
# Escape certain characters, such as \t by \\t, $ by \$, and { by {{.
rule true_sites_bed:
    input: "ChernobylTriosTruth/all_clean/{fam}c1.csv"
    output:
        output_dir + "true_sites_bed/truth_{fam}c1.bed"
    shell: """
         awk -v FS=',' -v OFS='\t' "{{ if(NR>1) print \$1,\$2-1,\$2}}" {input} | sort -u > {output}
    """


### Call JIGV for true sites
use rule call_JIGV as call_JIGV_for_truth with:
    input:
        bam=lambda w: expand(output_dir +"/fixed-rg/{id}.bam", id=[person.id for person in families[w.fam]]),
        ped="ped_files/{fam}.ped",
        ref=hg38_ref,
        sites= output_dir + "true_sites_bed/{caller}_{fam}c1.bed"
  
