
MulitQC_output_dir = config["multiqc"]["output_dir"]

rule multiqc:
    input: qc_output
    output:
        report(MulitQC_output_dir + "/multiqc.html", caption="../report/multiqc.rst", category="Quality control"),
        directory(MulitQC_output_dir),
    wrapper:
        "0.27.1/bio/multiqc"

optional_output.append( MulitQC_output_dir + "/multiqc.html" ) 