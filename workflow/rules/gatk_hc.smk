rule gatkhc_pb:
    input:
        bam=get_bam_by_subj,
        ref=genome,
        ref_dict=genome_dict
    output:
        gvcf=output_dir +"/gatkhc_pb/{subj}.gatk.g.vcf",
        gz=output_dir +"/gatkhc_pb/{subj}.gatk.g.vcf.gz",
        tbi=output_dir +"/gatkhc_pb/{subj}.gatk.g.vcf.gz.tbi"
    threads: config["threads"]["gatkhc_pb"]
    conda: "../envs/tabix.yaml"
    params: 
        singularity_cmd = config["parabricks"]["singularity_cmd"]
    benchmark:
        output_dir +"/benchmark/gatkhc_pb/{subj}.tsv"
    shell: """
        {params.singularity_cmd} \
        pbrun haplotypecaller \
            --gvcf \
            --ref {input.ref} \
            --in-bam {input.bam} \
            --out-variants {output.gvcf}   
        bgzip -c {output.gvcf} > {output.gz}
        tabix {output.gz}     
    """

rule gatk_combine_gvcf:
    input: 
        vcfs=lambda w: expand(output_dir + "/gatkhc_pb/{id}.gatk.g.vcf.gz", id= [person.id for person in families[w.fam]]),
        tbis=lambda w: expand(output_dir + "/gatkhc_pb/{id}.gatk.g.vcf.gz.tbi", id= [person.id for person in families[w.fam]]),
        ref=genome,
        fai=genome_fai,
    output: output_dir+"/gatk_combine_gvcf/{fam}.combine_gvcf.g.vcf.gz"
    threads: config["threads"]["gatk_combine_gvcf"]
    conda: "../envs/gatk4.yaml"
    benchmark:
        output_dir + "/benchmark/gatk_combine_gvcf/{fam}.tsv"
    params: 
        v=lambda w, input: " -V ".join(input.vcfs)
    shell: """
        gatk CombineGVCFs \
            -V {params.v} \
            -R {input.ref} \
            -O {output} 
    """

rule gatk_genotype_gvcf_pb:
    input: 
        gvcf=output_dir+"/gatk_combine_gvcf/{fam}.combine_gvcf.g.vcf.gz",
        ref=genome
    output:
        vcf=output_dir+"/gatk_genotype_gvcf_pb/{fam}.genotype_gvcf.g.vcf",
        gz=output_dir+"/gatk_genotype_gvcf_pb/{fam}.genotype_gvcf.g.vcf.gz",
    threads: config["threads"]["gatk_genotype_gvcf_pb"]
    conda: "../envs/tabix.yaml"
    benchmark:
        output_dir + "/benchmark/gatk_genotype_gvcf_pb/{fam}.tsv"
    params: singularity_cmd = config["parabricks"]["singularity_cmd"]
    shell: """
        {params.singularity_cmd} \
        pbrun genotypegvcf \
            --in-gvcf {input.gvcf} \
            --ref {input.ref} \
            --out-vcf {output.vcf} 
        bgzip -@ {threads} -c {output.vcf}> {output.gz} && \
        tabix {output.gz}
    """

# drop CalculateGenotypePosteriors completely
# gatk CalculateGenotypePosteriors \
        #     -V {input.gvcf} \
        #     -O {output.gvcf} 
        #     -ped {input.ped} --skip-population-priors


if config["gatk_hc"]["gatk_hard_filter"]["enable"]:
    rule gatk_cgp:
        input: 
            gvcf=output_dir+"/gatk_genotype_gvcf_pb/{fam}.genotype_gvcf.g.vcf.gz",
            # ped=config["ped_dir"]+"/{fam}.ped",
            ref=genome
        output: 
            # gvcf=output_dir+"/gatk_cgp/{fam}.cgp.g.vcf.gz",
            vcf=output_dir+"/gatk_cgp/{fam}.cgp_norm.vcf.gz",
            tbi=output_dir+"/gatk_cgp/{fam}.cgp_norm.vcf.gz.tbi",
        threads: config["threads"]["gatk_cgp"]
        conda: "../envs/bcftools.yaml"
        benchmark:
            output_dir + "/benchmark/gatk_cgp/{fam}.tsv"
        params: filter=config["gatk_hc"]["gatk_hard_filter"]["filter"]
        resources: 
        shell: """
            bcftools norm -f {input.ref} -m -  {input.gvcf} | bcftools view -i'ALT!="*"' | bcftools filter -i '{params.filter}' -Oz -o {output.vcf}
            tabix -p vcf {output.vcf}
        """

else: 
    rule gatk_cgp:
        input: 
            gvcf=output_dir+"/gatk_genotype_gvcf_pb/{fam}.genotype_gvcf.g.vcf.gz",
            ped=config["ped_dir"]+"/{fam}.ped",
            ref=genome
        output: 
            vcf=output_dir+"/gatk_cgp/{fam}.cgp_norm.vcf.gz",
            tbi=output_dir+"/gatk_cgp/{fam}.cgp_norm.vcf.gz.tbi",
        threads: config["threads"]["gatk_cgp"]
        conda: "../envs/bcftools.yaml"
        benchmark:
            output_dir + "/benchmark/gatk_cgp/{fam}.tsv"
        resources: 
        shell: """
            bcftools norm -f {input.ref} -m -  {output.gvcf} | bcftools view -i'ALT!="*"' -Oz -o {output.vcf}
            tabix -p vcf {output.vcf}
        """