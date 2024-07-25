import uuid

rule deepvariant_pb:
    input:
        bam=get_bam,
        ref=genome
    output:
        gvcf=output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf",
        gz=output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf.gz",
        tbi=output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf.gz.tbi"
    threads: config["threads"]["deepvariant_pb"]
    benchmark:
        output_dir +"/benchmark/deepvariant_pb/{subj}.tsv"
    conda: "../envs/tabix.yaml"
    params: singularity_cmd = config["parabricks"]["singularity_cmd"]
    shell: """
        {params.singularity_cmd} \
        pbrun deepvariant \
        --gvcf  \
        --ref {input.ref} \
        --in-bam {input.bam} \
        --out-variants {output.gvcf} 
        bgzip -c {output.gvcf} > {output.gz}
        tabix {output.gz}
    """

rule glnexus_dv: 
    input:
        gvcf=lambda w: expand(output_dir +"/deepvariant_pb/{id}.deepvariant.g.vcf.gz", id=[person.id for person in families[w.fam]]),
        tbi=lambda w: expand(output_dir +"/deepvariant_pb/{id}.deepvariant.g.vcf.gz.tbi", id=[person.id for person in families[w.fam]]),
        ref=genome
    output:
        vcf = output_dir +"/glnexus/{fam}.dv_combined.vcf.gz"
    threads: config["threads"]["glnexus_dv"]
    benchmark:
        output_dir +"/benchmark/glnexus_dv/dv_{fam}.tsv"
    conda: "../envs/glnexus.yaml"
    params: 
        tempdir=temp(directory(output_dir +"/glnexus/GLnexus.DB_{}")).format(uuid.uuid4()),
        config=config["deepvariant"]["glnexus_config"] 
    shell: """
        glnexus_cli -t {threads} --config {params.config} --dir {params.tempdir} {input.gvcf} |  bcftools norm -f {input.ref} -m - -O z -o {output} 
        rm -fr {params.tempdir}
    """