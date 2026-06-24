import uuid

rule deepvariant_pb:
    input:
        bam=get_bam_by_subj,
        ref=genome
    output:
        gvcf=output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf"
    threads: config["threads"]["deepvariant_pb"]
    benchmark:
        output_dir +"/benchmark/deepvariant_pb/{subj}.tsv"
    params: singularity_cmd = config["parabricks"]["singularity_cmd"]
    shell: """
        {params.singularity_cmd} \
        pbrun deepvariant \
            --gvcf \
            --ref {input.ref} \
            --in-bam {input.bam} \
            --out-variants {output.gvcf}
    """

rule deepvariant_pb_index:
    input:
        gvcf=output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf"
    output:
        gz=output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf.gz",
        tbi=output_dir +"/deepvariant_pb/{subj}.deepvariant.g.vcf.gz.tbi"
    container: CONTAINERS["tabix"]
    shell: """
        bgzip -c {input.gvcf} > {output.gz}
        tabix {output.gz}
    """

rule glnexus_dv_raw:
    input:
        gvcf=lambda w: expand(output_dir +"/deepvariant_pb/{id}.deepvariant.g.vcf.gz", id=[person.id for person in families[w.fam]]),
        tbi=lambda w: expand(output_dir +"/deepvariant_pb/{id}.deepvariant.g.vcf.gz.tbi", id=[person.id for person in families[w.fam]]),
    output:
        bcf=temp(output_dir +"/glnexus/{fam}.dv_combined.bcf")
    threads: config["threads"]["glnexus_dv"]
    benchmark:
        output_dir +"/benchmark/glnexus_dv/dv_{fam}.tsv"
    container: CONTAINERS["glnexus"]
    params:
        tempdir=temp(directory(output_dir +"/glnexus/GLnexus.DB_{}")).format(uuid.uuid4()),
        config=config["deepvariant"]["glnexus_config"]
    shell: """
        glnexus_cli -t {threads} --config {params.config} --dir {params.tempdir} {input.gvcf} > {output.bcf}
        rm -fr {params.tempdir}
    """

rule glnexus_dv:
    input:
        bcf=output_dir +"/glnexus/{fam}.dv_combined.bcf",
        ref=genome
    output:
        vcf = output_dir +"/glnexus/{fam}.dv_combined.vcf.gz",
        tbi = output_dir +"/glnexus/{fam}.dv_combined.vcf.gz.tbi"
    container: CONTAINERS["bcftools"]
    shell: """
        bcftools norm -f {input.ref} -m - {input.bcf} -O z -o {output.vcf}
        tabix -p vcf {output.vcf}
    """
