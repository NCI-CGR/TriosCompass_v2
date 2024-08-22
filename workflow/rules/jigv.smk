### Call JIGV
rule call_JIGV:
    input:
        bam=get_bams_by_family,
        ped = config["ped_dir"]+"/{fam}.ped",
        ref=genome,
        sites= output_dir +"/GATK_DV/D_and_G.{fam}.dnm.vcf.gz"
    output:
        html= report(
            output_dir +"/call_JIGV/{fam}.JIGV.html",
            caption="../report/JIGV.rst",
            category="De novo mutations",
            subcategory="Visulization",
            labels={
                "Family": "{fam}",
                "Desc": "JIGV snapshots",
                "File type": "html",
            }
        )
    benchmark:
        output_dir +"/benchmark/call_JIGV/{fam}.tsv"
    params: 
        proband=lambda w: [person.id for person in families[w.fam] if families[w.fam].get_father(person) ][0]
    singularity: "docker://chrisamiller/jigv:0.1.10"
    shell: """
        # jigv is at / in chrisamiller/jigv:0.1.10
        /jigv  \
            --fasta {input.ref} \
            --sample {params.proband} \
            --ped {input.ped} \
            --sites {input.sites} \
            {input.bam} > {output.html}
    """