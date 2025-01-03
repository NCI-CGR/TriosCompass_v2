import peppy
from itertools import compress
import os
import peds

# Process fastq sequences using fastp, fastqc, and fastq_screen

pepfile: os.path.join(workflow.basedir, "..", config["pepfile"])

### Disable pepschema, due to:
# Could not load remote file: http://schema.databio.org/pep/2.0.0.yaml. Original exception: URLError(TimeoutError(110, 'Connection timed out'))
# pepschema: config["pepschema"]

samples = pep.sample_table
subjs = list(set(samples["SAMPLE_ID"].tolist())) # 
sample_bam_dict=dict(zip(samples.SAMPLE_ID, samples.BAM))

# print (sample_bam_dict)

if config["bam_input"]["reset_RG"]:
    rule replace_rg:
        input: lambda w: sample_bam_dict[w.subj] 
        output:
            output_dir + "/fixed-rg/{subj}.bam"
        params:
            extra="--RGLB lib1 --RGPL illumina --RGPU {subj} --RGSM {subj} --CREATE_INDEX true "
        benchmark:
            output_dir +"/benchmark/replace_rg/{subj}.tsv"
        threads: config["threads"]["replace_rg"]
        wrapper:
            "v1.25.0/bio/picard/addorreplacereadgroups"
    
    def get_bam(subj):
        return output_dir +"/fixed-rg/{}.bam".format(subj) 


else:
    def get_bam(subj):
        return sample_bam_dict[subj] 

    # print(get_bam("SC501095"))

