# ==================================
# @File    :   multiqc.snakefile
# @Time    :   2025/11/17 15:22:33
# @Author  :   Xin Dong
# @Contact :   xindong9511@gmail.com
# @License :   (C)Copyright 2020-2025, XinDong
# ==================================

# MultiQC report generation

def multiqc_targets(wildcards):
    ls = []
    ls.append(f"{RES_PATH}/multiqc/multiqc_report.zip")
    ls.append(f"{RES_PATH}/multiqc/multiqc_report.html")
    return ls

def get_multiqc_input(wildcards):
    ls = all_targets(wildcards)
    ls.remove(f"{RES_PATH}/multiqc/multiqc_report.zip")
    ls.remove(f"{RES_PATH}/multiqc/multiqc_report.html")
    inputs = ls
    return inputs


rule multiqc_report:
    input:
        get_multiqc_input
    output:
        f"{RES_PATH}/multiqc/multiqc_report.zip",
        f"{RES_PATH}/multiqc/multiqc_report.html"
    params:
        output_dir=f"{RES_PATH}/multiqc/",
    log:
        f"{RES_PATH}/logs/multiqc/multiqc_report.log"
    message: "MULTIQC: Generating MultiQC report"
    shell:
        """
        multiqc {RES_PATH} -o {params.output_dir} > {log} 2>&1
        """