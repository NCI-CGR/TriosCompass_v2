# ref: https://github.com/snakemake-workflows/dna-seq-varlociraptor/blob/master/workflow/rules/ref.smk



rule bwa_index:
    input:
        genome,
    output:
        idx=multiext(genome, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    threads:
        config["threads"]["bwa_index"]
    wrapper:
        "v2.3.2/bio/bwa/index"

rule genome_dict:
    input:
        genome,
    output:
        genome_dict,
    # conda:
    #     "../envs/samtools.yaml"
    singularity: "docker://euformatics/samtools:1.19.2"
    shell:
        "samtools dict {input} > {output} "

rule genome_faidx:
    input:
        genome,
    output:
        genome_fai,
    # conda:
    #     "../envs/samtools.yaml"
    singularity: "docker://euformatics/samtools:1.19.2"
    shell:
        "samtools faidx {input}  "