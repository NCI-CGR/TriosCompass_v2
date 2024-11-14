### fix DP (for strelka) and AD
rule call_dnm_dv: 
    input: 
        vcf = output_dir +"/glnexus/{fam}.dv_combined.vcf.gz",
        ped = config["ped_dir"]+"/{fam}.ped",
        interval= config["call_dnm"]["interval"]
    output:
        vcf= output_dir +"/slivar/DV_{fam}.dnm.vcf",
        tmp = temp(output_dir +"/slivar/DV_{fam}.tmp.vcf.gz"),
        gz = output_dir +"/slivar/DV_{fam}.dnm.vcf.gz"
    benchmark:
        output_dir +"/benchmark/slivar/DV_{fam}.tsv"
    params: 
        # params: min_gq=20, min_parent_gq=20, min_dp=20, max_dp=250
        min_00_gq=config["call_dnm"]["dv"]["min_00_gq"], 
        min_01_gq=config["call_dnm"]["dv"]["min_01_gq"], 
        min_dp=config["call_dnm"]["dv"]["min_dp"],
        max_dp=config["call_dnm"]["dv"]["max_dp"],
        max_err=config["call_dnm"]["max_err"]
    conda: "../envs/slivar.yaml"
    shell: """
        
        slivar expr  \
            --vcf {input.vcf} \
            --ped  {input.ped} \
            --pass-only \
            --out-vcf {output.vcf} \
            --trio "denovo:( \
                ( \
                    (variant.CHROM == 'chrX' && kid.sex=='male') && \
                    kid.PL[0]>20 && kid.PL[1]>20 && kid.PL[2]==0 && kid.AB > 0.98  \
                ) || \
                ( \
                    (!(variant.CHROM == 'chrX' && kid.sex=='male')) && \
                    kid.PL[0]>20 && kid.PL[1]==0 && kid.PL[2]>20 && kid.AB > 0.25 && kid.AB < 0.75 \
                ) \
                ) &&  (kid.AD[0]+kid.AD[1]) >= {params.min_dp}/(1+(variant.CHROM == 'chrX' && kid.sex == 'male' ? 1 : 0)) && \
                (kid.AD[0]+kid.AD[1]) < {params.max_dp}/(1+(variant.CHROM == 'chrX' && kid.sex == 'male' ? 1 : 0)) && \
                mom.PL[0]==0 && mom.PL[1]>20 && mom.PL[2]>20 && dad.PL[0]==0 && dad.PL[1]>20 && dad.PL[2]>20 \
                    &&  (mom.AD[1]/(mom.AD[0]+mom.AD[1])) < {params.max_err} \
                    &&  (dad.AD[1]/(dad.AD[0]+dad.AD[1])) < {params.max_err} \
                    && kid.GQ >= {params.min_01_gq} && mom.GQ >= {params.min_00_gq} && dad.GQ >= {params.min_00_gq} \
                    && (mom.AD[0]+mom.AD[1]) >= {params.min_dp} && (mom.AD[0]+mom.AD[1]) < {params.max_dp} && (dad.AD[0]+dad.AD[1]) >= {params.min_dp}/(1+(variant.CHROM == 'chrX' ? 1 : 0)) && (dad.AD[0]+dad.AD[1]) < {params.max_dp}/(1+(variant.CHROM == 'chrX' ? 1 : 0)) "
        bgzip -c {output.vcf} > {output.tmp}
        tabix {output.tmp}
        bcftools view -R {input.interval} {output.tmp} -O z -o {output.gz}
        tabix {output.gz}
    """

    #     shell: """
        
    #     slivar expr  \
    #         --vcf {input.vcf} \
    #         --ped  {input.ped} \
    #         --pass-only \
    #         --out-vcf {output.vcf} \
    #         --trio "denovo:( \
    #             ( \
    #                 (variant.CHROM == 'chrX' && kid.sex=='male') && \
    #                 kid.hom_alt && kid.AB > 0.98  \
    #             ) || \
    #             ( \
    #                 (!(variant.CHROM == 'chrX' && kid.sex=='male')) && \
    #                 kid.het && kid.AB > 0.25 && kid.AB < 0.75 \
    #             ) \
    #             ) &&  (kid.AD[0]+kid.AD[1]) >= {params.min_dp}/(1+(variant.CHROM == 'chrX' && kid.sex == 'male' ? 1 : 0)) && \
    #             (kid.AD[0]+kid.AD[1]) < {params.max_dp}/(1+(variant.CHROM == 'chrX' && kid.sex == 'male' ? 1 : 0)) && \
    #             mom.hom_ref && dad.hom_ref \
    #                 &&  (mom.AD[1]/(mom.AD[0]+mom.AD[1])) < {params.max_err} \
    #                 &&  (dad.AD[1]/(dad.AD[0]+dad.AD[1])) < {params.max_err} \
    #                 && kid.GQ >= {params.min_01_gq} && mom.GQ >= {params.min_00_gq} && dad.GQ >= {params.min_00_gq} \
    #                 && (mom.AD[0]+mom.AD[1]) >= {params.min_dp} && (mom.AD[0]+mom.AD[1]) < {params.max_dp} && (dad.AD[0]+dad.AD[1]) >= {params.min_dp}/(1+(variant.CHROM == 'chrX' ? 1 : 0)) && (dad.AD[0]+dad.AD[1]) < {params.max_dp}/(1+(variant.CHROM == 'chrX' ? 1 : 0)) "
    #     bgzip -c {output.vcf} > {output.tmp}
    #     tabix {output.tmp}
    #     bcftools view -R {input.interval} {output.tmp} -O z -o {output.gz}
    #     tabix {output.gz}
    # """

### Call dnm from gatk
use rule  call_dnm_dv as call_dnm_gatk with:
    input: 
        vcf=output_dir+"/gatk_cgp/{fam}.cgp_norm.vcf.gz",
        ped = config["ped_dir"]+"/{fam}.ped",
        interval= config["call_dnm"]["interval"]
    output: 
       vcf= output_dir +"/slivar/GATK_{fam}.dnm.vcf",
       tmp = temp(output_dir +"/slivar/GATK_{fam}.tmp.vcf.gz"),
       gz = output_dir +"/slivar/GATK_{fam}.dnm.vcf.gz"
    params: 
        min_00_gq=config["call_dnm"]["hc"]["min_00_gq"], 
        min_01_gq=config["call_dnm"]["hc"]["min_01_gq"], 
        min_dp=config["call_dnm"]["hc"]["min_dp"],
        max_dp=config["call_dnm"]["hc"]["max_dp"],
        max_err=config["call_dnm"]["max_err"]
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
    conda: "../envs/bcftools.yaml"
    shell: """
        bcftools merge --force-samples --threads 2 -m none {input.DV} {input.GATK} -Oz -o {output.gz}
        tabix {output.gz}
        bcftools filter -i "N_MISSING == 0" -Oz -o {output.both} {output.gz}
        bcftools filter -i "N_MISSING > 0" -Oz -o {output.one} {output.gz}
        tabix {output.both}
        tabix {output.one}
    """

### generate vcf file of the proband
rule dnm_vcf:
    input: output_dir +"/GATK_DV/D_and_G.{fam}.dnm.vcf.gz"
    output: 
        report(
            output_dir +"/dnm_vcf/{fam}.dnm.vcf.gz",
            caption="../report/dnm.rst",
            category="De novo mutations",
            subcategory="Prediction",
            labels={
                "Family": "{fam}",
                "File type": "VCF.gz",
                "Desc": "DNMs"
            }
        )
    conda: "../envs/bcftools.yaml"
    params: 
        proband = lambda w: CHILD_DICT[w.fam]
    shell: """
        bcftools view -s {params.proband} {input} | bcftools annotate -x ID  -I +"%CHROM:%POS:%REF:%ALT" -Oz -o {output}
        tabix -p vcf {output}
    """

rule dnm_vcf_summary:
    input:
        expand(output_dir +"/dnm_vcf/{fam}.dnm.vcf.gz", fam=fam_ids)
    output:
        report(
            output_dir + "/dnm_vcf_summary/DNM_summary.txt",
            category="De novo mutations",
            subcategory="Summary",
            labels={
                "Desc": "Count of DNMs",    
                "File type": "txt",
            }
        )
    shell: """
        (echo -e "TrioID\tDNM_Cnt" && ls {input} | parallel 'id=$(basename {{/}} .dnm.vcf.gz);cnt=$(zgrep -v "^#" {{}} |wc -l); echo -e "$id\t$cnt" ')  > {output}
    """