########################################################
# append rules for dnSTR calling
########################################################
split_total = config["dnSTR"]["split_n"]
CHUNKS =[str(x).zfill(5) for x in range(split_total)]

### split bed files first
# skip chrX for the time being
if config["dnSTR"]["hipstr"]["enable"]:
    include: "hipstr.smk"

    ### Note that gangstr output is .vcf (not vcf.gz)
    rule dumpstr_call:
        input: output_dir+"/{caller}/{chunk}.vcf.gz"
        output:
            multiext(output_dir+ "/dumpstr_call/{caller}_{chunk}", ".vcf", ".loclog.tab", ".samplog.tab")
        params:
            prefix = output_dir+"/dumpstr_call/{caller}_{chunk}",
            args =  lambda w: config["dnSTR"]["hipstr"]["dumpstr_call_args"] if (w.caller=='hipstr') else "",
        benchmark:
            output_dir + "/benchmark/dumpstr_call/{caller}_{chunk}.tsv"
        conda: "../envs/trtools.yaml"
        # container: "docker://quay.io/biocontainers/trtools:latest"
        shell: """
            mkdir -p $(dirname {params.prefix})	
            dumpSTR \
                --vcf {input} \
                --out {params.prefix} \
                {params.args}
        """

    rule merge_ped:
        input:
            expand(ped_dir + "/{fam}.ped", fam=fam_ids)
        output:
            output_dir + "/merge_ped/all.ped"
        shell: """
            cat {input} | awk '{{if($3!=0 && $4!=0) print $0}}'   | sed -e 's/^\(t....\)c./\\1/' > {output}
        """

    ### The same rule for both hipstr and gangtr
    rule dumpstr_locus:
        input: 
            vcf= output_dir + "/dumpstr_call/{caller}_{chunk}.vcf",
            dup_reg = config["dnSTR"]["dup_reg"]
        output:
            multiext(output_dir + "/dumpstr_locus/{caller}_{chunk}", ".vcf", ".loclog.tab", ".samplog.tab")
        params:
            prefix = output_dir + "/dumpstr_locus/{caller}_{chunk}"
        conda: "../envs/trtools.yaml"
        # container: "docker://quay.io/biocontainers/trtools:latest"
        shell: """
            dumpSTR --min-locus-hwep 0.00001 --min-locus-callrate 0.8 \
                --filter-regions {input.dup_reg} \
                    --filter-regions-names SEGDUP \
                    --vcf {input.vcf} \
                    --out {params.prefix}
        """

    rule vcf_index: 
        input: output_dir + "/dumpstr_locus/{caller}_{chunk}.vcf"
        output: 
            gz= output_dir + "/dumpstr_locus/{caller}_{chunk}.vcf.gz",
            tbi= output_dir + "/dumpstr_locus/{caller}_{chunk}.vcf.gz.tbi"
        conda: "../envs/bcftools.yaml"
        shell: """
            bcftools sort {input} | bgzip -c  > {output.gz}
            tabix -p vcf {output.gz}
        """

    rule monstr:
        input: 
            vcf = output_dir + "/dumpstr_locus/{caller}_{chunk}.vcf.gz",
            ped = output_dir + "/merge_ped/all.ped"
        output: 
            multiext( output_dir + "/monstr/{caller}_{chunk}", ".all_mutations.tab", ".locus_summary.tab")
        benchmark:
            output_dir + "/benchmark/monstr/{caller}_{chunk}.tsv"
        singularity: "docker://gymreklab/monstr"
        params:
            prefix=output_dir + "/monstr/{caller}_{chunk}",
            arg = lambda w: config["dnSTR"]["hipstr"]["monstr_filter"] if (w.caller=='hipstr') else " ",
        shell: """
            # https://unix.stackexchange.com/questions/330660/preventing-grep-from-causing-premature-termination-of-bash-e-script
            empty_vcf=$(zgrep -v "^#"  {input.vcf} || true; )  
            
            if [ -z "$empty_vcf" ]
            then
                mkdir -p $(dirname {params.prefix})
                touch {output}
            else  
                MonSTR  \
                    --strvcf {input.vcf} \
                    --fam {input.ped} \
                    --out {params.prefix} \
                    --min-score 0.8 \
                    --min-coverage 10 \
                    {params.arg} 
            fi
                
        """

    rule merge_monstr:
        input: 
            mutation=expand(output_dir + "/monstr/{{caller}}_{chunk}.all_mutations.tab", chunk=CHUNKS),
            summary=expand(output_dir + "/monstr/{{caller}}_{chunk}.locus_summary.tab", chunk=CHUNKS)
        output: 
            mutation=output_dir + "/merge_monstr/{caller}.all_mutations.tab",
            summary=output_dir + "/merge_monstr/{caller}.locus_summary.tab"
        conda: "../envs/csvtk.yaml"
        shell: """
            csvtk concat -E -t {input.mutation} | csvtk sort -t -k chrom:N -k pos:n -k child -E - | csvtk filter -t -f "poocase>1" > {output.mutation}
            csvtk concat -E -t {input.summary} | csvtk sort -t -k chrom:N -k pos:n  -E - | csvtk filter -t -f "total_mutations>0" > {output.summary}
        """


    rule monstr_filter:
        input: output_dir + "/merge_monstr/{caller}.all_mutations.tab"
        output: 
            filtered=report(
                output_dir + "/monstr_filter/{caller}.filtered.tab",
                caption="../report/dnSTR.rst",
                category="De novo STRs",
                subcategory="MonSTR predictions",
                labels={
                    "STR caller": "{caller}",
                    "Desc": "MonSTR predictions for all trios",
                    "File type": "tab/tsv"
                }
            ),
            log=output_dir + "/monstr_filter/{caller}.filtered.log"
        singularity: "docker://gymreklab/monstr"
        shell: """
            python3 /STRDenovoTools/scripts/qc_denovos.py --all-mutations-file {input} \
                    --filtered-mutations-file {output.filtered} \
                    --log-file {output.log} \
                    --filter-denovos-child 5 \
                    --filter-loc-denovos 5 \
                    --filter-posterior 0.8 
        """

    rule childDict2table:
        output: output_dir + "/dnSTR_summary/childDict.csv"
        params: dict = CHILD_DICT
        run: 
            import csv

            with open(str(output), 'w', newline='') as csvfile:
                fieldnames = ['child', 'trio_id']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

                writer.writeheader()
                for key in params.dict:
                    writer.writerow({'trio_id': key, 'child': params.dict[key]})

    rule dnSTR_summary:
        input: 
            tab=output_dir + "/monstr_filter/{caller}.filtered.tab",
            csv=output_dir + "/dnSTR_summary/childDict.csv"
        output: 
            report(
                output_dir + "/dnSTR_summary/{caller}.dnSTR_summary.txt",
                caption="../report/dnSTR.rst",
                category="De novo STRs",
                subcategory="Summary of MonSTR predictions",
                labels={
                    "STR caller": "{caller}",
                    "Desc": "Count of dnSTRs",
                    "File type": "txt"
                }
            )
        conda: "../envs/csvtk.yaml"
        shell: """
            csvtk join -t --na 0 -L -f child <(csvtk csv2tab {input.csv}) <(csvtk summary -t -g child -f pos:countn  {input.tab} | csvtk rename -t -f 1,2 -n child,countn) > {output}
        """

    ### Join hipstr and gangstr if both called
    # rule joint_STR: 
    #     input: expand(output_dir + "/monstr_filter/{caller}.filtered.tab", caller=['hipstr','gangstr'])
    #     output: output_dir + "/joint_STR/gangstr_hipstsr.final.tab"
    #     conda: "../envs/csvtk.yaml"
    #     shell: """
    #         csvtk join -t -f "chrom,pos,child,child_gt,mat_gt,pat_gt" -p {input} > {output}
    #     """

    optional_output.append( expand(output_dir + "/dnSTR_summary/{caller}.dnSTR_summary.txt", caller=['hipstr']) )    

