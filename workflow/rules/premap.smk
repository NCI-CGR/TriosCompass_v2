import peppy
from itertools import compress
import os

# Process fastq sequences using fastp, fastqc, and fastq_screen

pepfile: os.path.join(workflow.basedir, "..", config["pepfile"])
# validate_project() got an unexpected keyword argument 'exclude_case'
pepschema: config["pepschema"]


samples = pep.sample_table
flowcells = list(set(samples["FLOWCELL"].tolist()))
subjs = list(set(samples["SAMPLE_ID"].tolist())) # 
ids = samples["sample_name"].tolist() # ids for each fastq file
subj_flowcell_dict=samples[['SAMPLE_ID','FLOWCELL']].drop_duplicates().groupby('SAMPLE_ID')['FLOWCELL'].apply(list).to_dict()
 

# fq2bam                  - Run bwa mem, co-ordinate sorting, marking duplicates, and Base Quality Score Recalibration
rule fastp: 
    input: 
        R1= lambda w:  samples.loc[w.sample_name]["R1"],
        R2= lambda w:  samples.loc[w.sample_name]["R2"]
    output: 
        R1=output_dir+"/fastp/{sample_name}.R1.fastp.fastq.gz",
        R2=output_dir+"/fastp/{sample_name}.R2.fastp.fastq.gz",
        html=output_dir+"/fastp/{sample_name}.html",
        json=output_dir+"/fastp/{sample_name}.json"
    benchmark:
        output_dir + "/benchmark/fastp/{sample_name}.tsv"
    threads: config["threads"]["fastp"]
    conda:
        "../envs/fastp.yaml"
    shell: """
        fastp \
                --in1 {input.R1} \
                --in2 {input.R2} \
                --out1 {output.R1} \
                --out2 {output.R2} \
                --report_title {wildcards.sample_name} \
                --json {output.json} \
                --html {output.html} \
                --thread {threads} \
                --qualified_quality_phred 25 --n_base_limit 10 --average_qual 25 \
                --length_required 50 --low_complexity_filter
    """

if config["fastq_input"]["fastqc"]["enable"]:
    rule fastqc:
        input: output_dir+"/fastp/{fid}.fastp.fastq.gz"
        output:
            html=output_dir+"/fastqc/{fid}.html",
            zip=output_dir+"/fastqc/{fid}_fastqc.zip" 
        params:
            extra = "--quiet"
        log:
            "logs/fastqc/{fid}.log"
        threads: 1
        resources:
            mem_mb=10000,
            runtime=600
        wrapper:
            "v2.1.1/bio/fastqc"     
    
    optional_output.append( expand(output_dir+"/fastqc/{id}.{reads}_fastqc.zip", id =ids, reads=['R1','R2']))


#  fastq_screen: TODO

# if config["pre_mapping"]["fastq_screen"]["enable"]:
#     rule fastq_screen:
#         input:
#             output_dir+"/fastp/{id}.R1.fastp.fastq.gz"
#         output:
#             txt=output_dir+"/fastq_screen/{id}.fastq_screen.txt",
#             png=output_dir+"/fastq_screen/{id}.fastq_screen.png"
#         params:
#             fastq_screen_config="ref/fastq_screen.abs.conf",
#             subset=100000,
#             aligner='bwa'
#         wrapper:
#             "v2.1.1/bio/fastq_screen"
    
#     optional_output.append(rule.fastq_screen.output)