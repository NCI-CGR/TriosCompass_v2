rule config_manta_joint_call:
    input:
        bams = get_bams_by_family,
        ref = genome
    output:
        output_dir + '/manta_joint/{fam}/runWorkflow.py',
        output_dir + '/manta_joint/{fam}/runWorkflow.py.config.pickle'
    params:
        prefix = output_dir + '/manta_joint/{fam}',
        bams = lambda w, input: " --bam ".join(input.bams)
    singularity: 'docker://szarate/manta:v1.6.0'
    shell:
        'configManta.py \
            --bam {params.bams} \
            --referenceFasta {input.ref} \
            --runDir {params.prefix}'

rule run_manta_joint_call:
    input:
        cmd = output_dir + '/manta_joint/{fam}/runWorkflow.py',
        pickle = output_dir + '/manta_joint/{fam}/runWorkflow.py.config.pickle'
    output:
        multiext(output_dir + '/manta_joint/{fam}/results/variants/', 'candidateSmallIndels.vcf.gz','candidateSmallIndels.vcf.gz.tbi', 'candidateSV.vcf.gz', 'candidateSV.vcf.gz.tbi','diploidSV.vcf.gz','diploidSV.vcf.gz.tbi'),
        multiext(output_dir + '/manta_joint/{fam}/results/stats/', 
        'alignmentStatsSummary.txt', 'svCandidateGenerationStats.tsv','svCandidateGenerationStats.xml','svLocusGraphStats.tsv')
    threads: config["threads"]["manta_call"]
    benchmark:
        "benchmarks/run_manta_joint_call/{fam}.tsv"
    singularity: 'docker://szarate/manta:v1.6.0'
    shell:
        '{input.cmd} -m local -j {threads}'

rule dnSV_manta:
    input:
        vcf = output_dir + '/manta_joint/{fam}/results/variants/diploidSV.vcf.gz',
        ped = config["ped_dir"]+"/{fam}.ped"
    output:
        vcf=report(
            output_dir + '/dnSV_manta/{fam}.dnSV.manta.vcf',
            caption="../report/dnSV_manta.rst",
            category="De novo SVs/Manta",
            subcategory="Predictions",
            labels={
                "Family": "{fam}",
                "Desc": "dnSV",
                "File type": "VCF"
            }
        )
    conda: "../envs/slivar.yaml"
    shell: """
        slivar expr  \
            --vcf {input.vcf} \
            --ped {input.ped}  \
            --pass-only \
            --out-vcf {output.vcf} \
            --trio "denovo:( \
                variant.FILTER == 'PASS' && \
                ( dad.hom_ref && mom.hom_ref && (kid.het || kid.hom_alt)) && \
                ( ( ('IMPRECISE' in INFO) && dad.PR[1]<=1 && mom.PR[1]<=1 && kid.PR[1]>=4) || \
                ( !('IMPRECISE' in INFO) && dad.SR[1]<=1 && mom.SR[1]<=1 && kid.SR[1]>=4)) \
            )"
    """

rule dnSV_manta_summary:
    input: expand(output_dir + "/dnSV_manta/{fam}.dnSV.manta.vcf", fam=fam_ids)
    output:
        report(
            output_dir + "/dnSV_manta_summary/dnSV_summary.txt",
            category="De novo SVs/Manta",
            subcategory="Summary",
            labels={
                "Desc": "Count of dnSVs",    
                "File type": "txt",
            }
        )
    shell: """
        # BND is counted twice in manta otuput, so we need to divide by 2
        (echo -e "TrioID\tdnSTR_Cnt" && ls {input} | parallel 'id=$(basename {{/}} .dnSV.manta.vcf);cnt=$(echo $(grep -v -c "^#" {{}}) - $(grep -v "^#" {{}} | grep -c "BND")/2 | bc ); echo -e "$id\t$cnt" ')  > {output}
    """

optional_output.append(expand(output_dir + "/dnSV_manta/{fam}.dnSV.manta.vcf", fam=fam_ids))
optional_output.append(output_dir + "/dnSV_manta_summary/dnSV_summary.txt")
