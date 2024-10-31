TriosCompass expects to use the parent folder of TriosCompass_v2 (the repo clone folder) as the working space, so as to separate the Snakemake workflow from the working space.  

Three configure files are required:
1. *Profile yaml file*: workflow/profiles/<PROFILE_NAME>/config.yaml
2. *Config yaml file* 
3. *Sample yaml file* 

The *config yaml file* can be specified by Snakemake command-line argument "\-\-configfile" or in the "*profile yaml file*, for example:

```yml
configfile: TriosCompass_v2/config/fullbam_config.yaml
snakefile: TriosCompass_v2/workflow/Snakefile
```

In turn, the *sample yaml file* is specified in the *config yaml file*, to [define sample input via PEPs](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#configuring-scientific-experiments-via-peps)

```yml
pepfile: "config/fullbam_pep.yaml"
pepschema: "../schemas/bam_schema.yaml"
```

+ Example of the PEP configure file for FASTQ input files
  + config/fastq_pep.yaml
    ```yml
    pep_version: 2.0.0
    sample_table: sample_fastq.csv
    
    # In manifest file, Sample_ID + Flowcell should be unique
    sample_modifiers:
      append:
        sample_name: "sn"
      derive:
        attributes: [sample_name]
        sources:
          sn: "{SAMPLE_ID}_{FLOWCELL}"
    ```

  + config/sample_fastq.csv 
    ```csv
    SAMPLE_ID,FLOWCELL,LANE,INDEX,R1,R2
    HG002,BH2JWTDSX5,1,CGGTTGTT-GTGGTATG,data/fq/HG002_NA24385_son_80X_R1.fq.gz,data/fq/HG002_NA24385_son_80X_R2.fq.gz
    HG003,BH2JWTDSX5,1,GCGTCATT-CAGACGTT,data/fq/HG003_NA24149_father_80X_R1.fq.gz,data/fq/HG003_NA24149_father_80X_R2.fq.gz
    HG004,BH2JWTDSX5,1,CTGTTGAC-ACCTCAGT,data/fq/HG004_NA24143_mother_80X_R1.fq.gz,data/fq/HG004_NA24143_mother_80X_R2.fq.gz
    ```
  + workflow/schemas/fastq_schema.yaml 
(schemas to validate config/fastq_pep.yaml) 
    ```yml
    workflow/schemas/fastq_schema.yaml 
    description: A example schema for a pipeline.
    imports:
      - http://schema.databio.org/pep/2.0.0.yaml
      # - TriosCompass_v2/workflow/schemas/2.0.0.yaml
      
    properties:
      samples:
        type: array
        items:
          type: object
          properties:
            SAMPLE_ID:
              type: string
              description: "sample id"
            FLOWCELL:
              type: string
              description: "Flowcell"
            INDEX:
              type: string
              description: "Library index"
            LANE:
              type: string
              description: "Lane number in flowcell"
              enum: ["1", "2"]
            R1:
              type: string
              description: "path to the R1 fastq file"
            R2:
              type: string
              description: "path to the R2 fastq file"
          required:
            - FLOWCELL
            - SAMPLE_ID
            - INDEX
            - R1
            - R2
    ```
+ Example of the PEP configure file for BAM input files
  + config/bam_pep.yaml
    ```yml
    pep_version: 2.0.0
    sample_table: sample_bam.csv
    
    
    sample_modifiers:
        append:
        sample_name: "sn"
        derive:
        attributes: [sample_name]
        sources:
            sn: "{SAMPLE_ID}"
    ```

  + config/sample_bam.csv 
    ```csv
    SAMPLE_ID,BAM
    HG002,sorted_bam/HG002_NA24385_son_80X.bam
    HG003,sorted_bam/HG003_NA24149_father_80X.bam
    HG004,sorted_bam/HG004_NA24143_mother_80X.bam
    ```
  + workflow/schemas/fastq_schema.yaml 
(schemas to validate config/bam_pep.yaml) 
    ```yml
    description: A example schema for a pipeline.
    imports:
      - http://schema.databio.org/pep/2.0.0.yaml
      # - TriosCompass_v2/workflow/schemas/2.0.0.yaml
      
    properties:
      samples:
        type: array
        items:
          type: object
          properties:
            SAMPLE_ID:
              type: string
              description: "sample id"
            BAM:
              type: string
              description: "path to the bam file"
          required:
            - SAMPLE_ID
            - BAM
    ``` 