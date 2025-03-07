### add optional qc 
if config["bam_qc"]["flagstat"]["enable"]:
    rule flagstat: 
        input: get_bam_by_subj
        output: 
            flagstat=output_dir +"/flagstat/{subj}.flagstat",
            idxstat=output_dir +"/flagstat/{subj}.idxstat"
        threads: config["threads"]["flagstat"]
        # conda:
        #     "../envs/samtools.yaml"
        singularity: "docker://euformatics/samtools:1.19.2"
        shell: """
            samtools flagstat -@ {threads} {input} > {output.flagstat}
            samtools idxstat {input} > {output.idxstat}
        """
    
    optional_output.append( expand(output_dir+"/flagstat/{subj}.flagstat", subj = subjs)) 
    optional_output.append( expand(output_dir+"/flagstat/{subj}.idxstat", subj = subjs))    
    qc_output.append(expand(output_dir+"/flagstat/{subj}.idxstat", subj = subjs))

if config["bam_qc"]["collectwgsmetrics"]["enable"]:
    rule collectwgsmetrics:
        input: 
            bam= get_bam_by_subj,
            ref= genome
        output: output_dir + "/collectwgsmetrics/{subj}.collect_wgs_metrics.txt"
        threads: config["threads"]["collectwgsmetrics"]
        benchmark:
            output_dir +"/benchmark/collectwgsmetrics/{subj}.tsv"
        singularity: "docker://quay.io/biocontainers/picard:2.27.3--hdfd78af_0"
        shell: """
            export  JAVA_TOOL_OPTIONS="-XX:+UseSerialGC -Xmx48G -Xms48G -XX:ParallelGCThreads=4 -XX:MaxRAMPercentage=90.0" 
            picard CollectWgsMetrics \
                I={input.bam} \
                O={output} \
                R={input.ref} 
        """
    
    optional_output.append( expand(output_dir+"/collectwgsmetrics/{subj}.collect_wgs_metrics.txt", subj = subjs))  
    qc_output.append(expand(output_dir+"/collectwgsmetrics/{subj}.collect_wgs_metrics.txt", subj = subjs))  

if config["bam_qc"]["collectmultiplemetrics"]["enable"]:
    rule collectmultiplemetrics:
        input: 
            bam = get_bam_by_subj,
            ref = genome
        output: output_dir +"/collectmultiplemetrics/{subj}/sequencingArtifact.pre_adapter_summary_metrics.txt"
        threads: config["threads"]["collectmultiplemetrics"]
        benchmark:
            output_dir +"/benchmark/collectmultiplemetrics/{subj}.tsv"
        params: 
            prefix= output_dir +"/collectmultiplemetrics/{subj}",
            singularity_cmd = config["parabricks"]["singularity_cmd"],

        shell: '''
            {params.singularity_cmd} \
            pbrun collectmultiplemetrics \
            --ref {input.ref} \
            --bam {input.bam} \
            --out-qc-metrics-dir {params.prefix} \
            --gen-all-metrics 
        '''
    
    optional_output.append( expand(output_dir+"/collectmultiplemetrics/{subj}/sequencingArtifact.pre_adapter_summary_metrics.txt", subj = subjs)) 
    qc_output.append( expand(output_dir+"/collectmultiplemetrics/{subj}", subj = subjs))
