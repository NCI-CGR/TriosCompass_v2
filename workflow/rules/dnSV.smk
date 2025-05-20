# call de novo SV using manta+GraphTyper2 and smoove
EXCLUDE_BED=config["dnSV"]["GraphTyper2+smoove"]["exclude_bed"]
SV_CALLERS=["smoove", "gt2"]

rule manta_create_run_script:
    input:
        bam = get_bam_by_subj,
        ref = genome
    output:
        output_dir + '/manta/{subj}/runWorkflow.py',
        output_dir + '/manta/{subj}/runWorkflow.py.config.pickle'
    params:
        prefix = output_dir + '/manta/{subj}'
    singularity: 'docker://szarate/manta:v1.6.0'
    shell:
        'configManta.py \
            --bam {input.bam} \
            --referenceFasta {input.ref} \
            --runDir {params.prefix}'

rule manta_call:
    input:
        cmd = output_dir + '/manta/{sample}/runWorkflow.py',
        pickle = output_dir + '/manta/{sample}/runWorkflow.py.config.pickle'
    output:
        multiext(output_dir + '/manta/{sample}/results/variants/', 'candidateSmallIndels.vcf.gz','candidateSmallIndels.vcf.gz.tbi', 'candidateSV.vcf.gz', 'candidateSV.vcf.gz.tbi','diploidSV.vcf.gz','diploidSV.vcf.gz.tbi'),
        multiext(output_dir + '/manta/{sample}/results/stats/', 
        'alignmentStatsSummary.txt', 'svCandidateGenerationStats.tsv','svCandidateGenerationStats.xml','svLocusGraphStats.tsv')
    threads: config["threads"]["manta_call"]
    benchmark:
        "benchmarks/manta_call/{sample}.tsv"
    singularity: 'docker://szarate/manta:v1.6.0'
    shell:
        '{input.cmd} -m local -j {threads}'

### make uniq SV id in manta
rule uniq_svid: 
    input:
        output_dir + '/manta/{sample}/results/variants/diploidSV.vcf.gz'
    output: 
        gz = output_dir + "/manta_sv/{sample}.vcf.gz",
        tbi = output_dir + "/manta_sv/{sample}.vcf.gz.tbi"
    threads: config["threads"]["uniq_svid"]
    conda: "../envs/bcftools.yaml"
    shell: """
        bcftools annotate --threads {threads} -I '{wildcards.sample}:%ID' {input} -Oz -o {output.gz}
        tabix -p vcf {output.gz}
    """

### make manta output list in the order: father, mother and child (the same order as the samples listed in the ped file)
rule trio_manta_list:
    input: 
        vcfs=lambda w: expand(output_dir + "/manta_sv/{id}.vcf.gz", id= [person.id for person in families[w.fam]])
    output: 
        output_dir + "/trio_manta_list/{fam}.manta_vcf.lst"
    shell: """
        # ls {input.vcfs} > {output}
        echo {input.vcfs} | awk '{{OFS="\\n"; $1=$1}}1' > {output}
    """



###  merge manta SV by svimmer
rule svimmer:
    input: 
        output_dir + "/trio_manta_list/{fam}.manta_vcf.lst"
    output:
        vcf = output_dir + "/svimmer/{fam}.merged.vcf",
    threads: config["threads"]["svimmer"]
    benchmark:
        "benchmarks/svimmer/{fam}.tsv"
    # params: 
    #     cmd = workflow.basedir + "/scripts/svimmer"
    singularity: "docker://bioinformatics/svimmer:latest"
    shell: """
        svimmer {input} chr{{1..22}} --ids --max_distance 50 --max_size_difference 100  > {output.vcf}
    """

rule vcfgz:
    input: output_dir + "/svimmer/{fam}.merged.vcf"
    output: 
        vcf = output_dir + "/svimmer/{fam}.merged.vcf.gz",
        tbi = output_dir + "/svimmer/{fam}.merged.vcf.gz.tbi"
    conda: "../envs/tabix.yaml"
    shell: """
        bgzip -c {input} > {output.vcf}
        tabix -p vcf {output.vcf}
    """
    
rule graphtyper:
    input: 
        bams = get_bams_by_family,
        ref = genome,
        vcf = output_dir + "/svimmer/{fam}.merged.vcf.gz"
    output:
        bam_lst = output_dir+"/graphtyper/{fam}/{fam}.bam.lst",
        vcf_lst = output_dir+"/graphtyper/{fam}/{fam}.vcf.lst",
        vcf = output_dir+"/graphtyper/{fam}.gt2.vcf.gz",
        tbi = output_dir+"/graphtyper/{fam}.gt2.vcf.gz.tbi"
    threads: config["threads"]["graphtyper"]
    benchmark:
        "benchmarks/graphtyper/{fam}.tsv"
    conda: "../envs/graphtyper.yaml"
    params: 
        prefix = output_dir + "/graphtyper/{fam}"
    shell: """
        # ls {input.bams} > {output.bam_lst}
        # To keep the sample in the order:father, mother, child
        echo {input.bams} | tr ' ' '\\n' > {output.bam_lst}

        seq 1 22 | sed -e "s/^/chr/" | parallel  'graphtyper genotype_sv {input.ref} {input.vcf} --sams={output.bam_lst} --threads=2 --region={{}} --output={params.prefix}'

        echo chr{{1..22}} | tr ' ' '\\n' | while read chrom; do if [[ ! -d {params.prefix}/${{chrom}} ]]; then continue; fi; find {params.prefix}/${{chrom}} -name "*.vcf.gz" | sort; done > {output.vcf_lst}

        bcftools concat --naive --file-list {output.vcf_lst} -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
    """

rule graphtyper_filter: 
    input: 
        vcf = output_dir+"/graphtyper/{fam}.gt2.vcf.gz"
    output:
        vcf = output_dir+"/graphtyper_filter/{fam}.gt2_filter.vcf.gz"
    threads: 4
    resources:
        runtime="4h",
        mem_mb=20000
    benchmark:
        "benchmarks/graphtyper_filter/{fam}.tsv"
    conda: "../envs/vcflib.yaml"
    shell: """
        vcffilter -f "( SVTYPE = BND & SVMODEL = AGGREGATED & QD > 20 & ( ABHet > 0.30 | ABHet < 0 ) & ( AC / NUM_MERGED_SVS ) < 10 & PASS_AC > 0 & PASS_ratio > 0.1 ) | ( SVTYPE = DEL & SVMODEL = AGGREGATED & QD > 12 & ( ABHet > 0.30 | ABHet < 0 ) & ( AC / NUM_MERGED_SVS ) < 25 & PASS_AC > 0 & PASS_ratio > 0.1 ) | ( SVTYPE = DUP & SVMODEL = AGGREGATED & QD > 5 & PASS_AC > 0 & ( AC / NUM_MERGED_SVS ) < 25 ) | ( SVTYPE = INS & SVMODEL = AGGREGATED & PASS_AC > 0 & ( AC / NUM_MERGED_SVS ) < 25 & PASS_ratio > 0.1 & ( ABHet > 0.25 | ABHet < 0 ) & MaxAAS > 4 ) | ( SVTYPE = INV & PASS_AC > 0 & ( AC / NUM_MERGED_SVS ) < 25 & PASS_ratio > 0.1 & ( ABHet > 0.25 | ABHet < 0 ) & MaxAAS > 4 )"  {input.vcf} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} 
    """

### smoove need unique output folder for each trio to have files like
# SC742298.histo
# output: {outdir}/{name}-smoove.vcf.gz
# using bam as input files
rule smoove:
    input: 
        # bams = lambda w: expand(output_dir + "/cram/{id}.cram", id= [person.id for person in families[w.fam]]),
        bams = get_bams_by_family,
        exclude_bed = EXCLUDE_BED,
        ref = genome
    output:
        vcf = output_dir + "/smoove/{fam}/{fam}-smoove.genotyped.vcf.gz",
        csi = output_dir + "/smoove/{fam}/{fam}-smoove.genotyped.vcf.gz.csi"
    threads: config["threads"]["smoove"]
    resources:
        runtime="24h",
        mem_mb=100000
    benchmark:
        "benchmarks/smoove/{fam}.tsv"
    params: 
        prefix = "{fam}",
        out_dir =  output_dir + "/smoove/{fam}"    
    container: "docker://brentp/smoove"
    shell: """
        smoove call -x --name {params.prefix} -o {params.out_dir} --exclude {input.exclude_bed} --fasta  {input.ref} -p {threads} --genotype {input.bams}
    """

### filter out dnSV using bcftools
# Note that samples are in the order: father, mother, proband
rule call_dnSV:
    input : 
        output_dir + "/smoove/{fam}/{fam}-smoove.genotyped.vcf.gz"
    output: 
        output_dir + "/dnSV/{fam}.smoove.dnSV.vcf"
    conda: "../envs/bcftools.yaml"
    shell: """
        bcftools view -i 'GT[2]="het" && GT[1]="RR" && GT[0]="RR" ' -O v  {input} -o {output}
    """

use rule call_dnSV as call_dnSV_gt2 with:
    input: 
        output_dir+"/graphtyper_filter/{fam}.gt2_filter.vcf.gz"
    output: 
        output_dir + "/dnSV/{fam}.gt2.dnSV.vcf"

use rule trio_manta_list as make_dnSV_list with:
    input:
        vcfs = expand(output_dir + "/dnSV/{{fam}}.{caller}.dnSV.vcf", caller=SV_CALLERS)
    output:     
        output_dir + "/joint_dnSV/{fam}.dnSV.lst"

rule joint_dnSV:
    input: 
        output_dir + "/joint_dnSV/{fam}.dnSV.lst"
    output:
        vcf = report(
            output_dir + "/joint_dnSV/{fam}.dnSV.vcf",
            caption="../report/dnSV.rst",
            category="De novo SVs/GraphTyper2+smoove",
            subcategory="Predictions",
            labels={
                "Family": "{fam}",
                "Desc": "dnSV",
                "File type": "VCF"
            }
        ),
        stats = output_dir + "/joint_dnSV/{fam}.dnSV.stats",
        others = multiext(output_dir+ "/joint_dnSV/{fam}.dnSV.stats", "_CHR", "support")
    conda: "../envs/survivor.yaml"
    shell: """
        SURVIVOR merge {input} 1000 2 0 0 0 50 {output.vcf}
        SURVIVOR stats {output.vcf} 50 -1 -1 {output.stats}
    """

rule dnSV_summary:
    input: expand(output_dir + "/joint_dnSV/{fam}.dnSV.vcf", fam=fam_ids)
    output:
        report(
            output_dir + "/dnSV_summary/dnSV_summary.txt",
            category="De novo SVs/GraphTyper2+smoove",
            subcategory="Summary",
            labels={
                "Desc": "Count of dnSVs",    
                "File type": "txt",
            }
        )
    shell: """
        (echo -e "TrioID\tdnSTR_Cnt" && ls {input} | parallel 'id=$(basename {{/}} .dnSV.vcf);cnt=$(zgrep -v -c "^#" {{}} ); echo -e "$id\t$cnt" ')  > {output}
    """

optional_output.append(expand(output_dir + "/joint_dnSV/{fam}.dnSV.vcf", fam=families))
optional_output.append(output_dir + "/dnSV_summary/dnSV_summary.txt")
