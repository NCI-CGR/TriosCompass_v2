
MulitQC_output_dir = config["multiqc"]["output_dir"]

rule multiqc:
    input: qc_output
    output:
        report(MulitQC_output_dir + "/multiqc.html", caption="../report/multiqc.rst", category="Quality control"),
        directory(MulitQC_output_dir),
    wrapper:
        "v1.4.0/bio/multiqc"

# rule multiqc:
#     input: qc_output
#     output:
#          html=report(MulitQC_output_dir + "/multiqc.html", caption="../report/multiqc.rst", category="Quality control"),
#          dir=directory(MulitQC_output_dir),
#     conda: "../envs/multiqc.yaml"
#     shell: """
#         multiqc --title QC --filename $(basename {output.html}) --outdir {output.dir} {input}  --force
#     """

optional_output.append( MulitQC_output_dir + "/multiqc.html" ) 