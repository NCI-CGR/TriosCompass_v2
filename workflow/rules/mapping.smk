### One sample may have data from multiple flowcells 
rule fq2bam:
    input: 
        #R1s=lambda w: expand(output_dir+"/fastp/{id}.R1.fastp.fastq.gz", id=[ w.SAMPLE_ID + '_'+ flowcell + '_L' + str(lane) for flowcell,lane in subj_flowcell_dict[w.SAMPLE_ID]]),
        R1s = lambda w: expand(output_dir+"/fastp/{id}.R1.fastp.fastq.gz", 
                       id=[ w.SAMPLE_ID + '_'+ flowcell + '_L' + str(lane) 
                            for flowcell, lane in subj_flowcell_dict[w.SAMPLE_ID]]),
        #R2s=lambda w: expand(output_dir+"/fastp/{id}.R2.fastp.fastq.gz", id=[ w.SAMPLE_ID + '_'+ flowcell + '_L' + str(lane) for flowcell,lane in subj_flowcell_dict[w.SAMPLE_ID]]),
        R2s = lambda w: expand(output_dir+"/fastp/{id}.R2.fastp.fastq.gz", 
                       id=[ w.SAMPLE_ID + '_'+ flowcell + '_L' + str(lane) 
                            for flowcell, lane in subj_flowcell_dict[w.SAMPLE_ID]]),
        ref=genome,
        idx=rules.bwa_index.output
    output: output_dir + "/fq2bam/{SAMPLE_ID}.bam"
    threads: config["threads"]["fq2bam"]
    params: 
       fqs=lambda w, input: " --in-fq ".join([r1+' '+r2+' ' + rg for r1,r2, rg in zip(input.R1s, input.R2s, ["'@RG\\tPL:ILLUMINA\\tID:{FLOWCELL}_{LANE}\\tSM:{SAMPLE_ID}\\tPU:{SAMPLE_ID}_{FLOWCELL}\\tLB:{SAMPLE_ID}_{INDEX}'" .format (**samples.loc[id]) for id in list(compress(ids, [ id == w.SAMPLE_ID for id in samples['SAMPLE_ID']]))  ])]),
       singularity_cmd = config["parabricks"]["singularity_cmd"]
    benchmark: 
        output_dir + "/benchmark/fq2bam/{SAMPLE_ID}.tsv"
    shell: """
        {params.singularity_cmd} \
        pbrun fq2bam \
            --ref {input.ref} \
            --in-fq {params.fqs} \
            --out-bam {output} 
    """

def get_bam(subj):
    return output_dir +"/fq2bam/{}.bam".format(subj) 

