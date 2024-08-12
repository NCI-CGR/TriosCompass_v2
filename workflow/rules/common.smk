import glob
from os import path

from snakemake.utils import validate


# validate(config, schema="../schemas/fastq_schema.yaml")

# location of genome fasta file with .fa or .fasta ext
genome = config["ref"]["sequence"]

refFile = os.path.basename(genome)
refDir = os.path.dirname(genome)
genome_prefix = os.path.splitext(refFile)[0]
genome_dict = refDir + '/' + genome_prefix + '.dict'
genome_fai = refDir + '/' + refFile + '.fai'

optional_output = list()
qc_output = list()


output_dir = config["output_dir"]

def get_bam_by_subj(wildcards):
    return get_bam(wildcards.subj)

def get_bams_by_family(wildcards):
    rv = [get_bam(person.id) for person in families[wildcards.fam]]
    return rv

def get_child_bam_by_family(wildcards):
    return get_bam(CHILD_DICT[wildcards.fam])
