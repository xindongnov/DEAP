# ==================================
# @File    :   rsem.snakefile
# @Time    :   2023/07/14 21:28:24
# @Author  :   Xin Dong @ref:RIMA
# @Contact :   xindong9511@gmail.com
# @License :   (C)Copyright 2020-2023, XinDong
# ==================================



import pandas as pd

# metadata = pd.read_csv(config["metasheet"], index_col=0, sep=',')

def merge_sep_inputs(inputs):
    inputs_format = ' -f '.join(str(i) for i in list(inputs)[0])
    return inputs_format   

def rsem_quantification_targets(wildcards):
    ls = []
    for sample in config["samples"]:
        pass
    #     ls.append("analysis/rsem/%s/%s.isoforms.results" % (sample, sample))
    #     ls.append("analysis/rsem/%s/%s.genes.results" % (sample, sample))
    # ls.append("analysis/rsem/rsem_tpm_iso_matrix.csv")
    # ls.append("analysis/rsem/rsem_tpm_gene_matrix.csv")
    return ls


rule rsem_quantification:
    input:
        "analysis/star/{sample}/{sample}.transcriptome.bam" #just to make sure STAR output is available before STAR_Fusion
    output:
        rsem_transcript_out = "analysis/rsem/{sample}/{sample}.isoforms.results",
        rsem_gene_out = "analysis/rsem/{sample}/{sample}.genes.results"
    log:
        "analysis/rsem/{sample}/{sample}.rsem.log"
    message: 
        "Running RSEM on {wildcards.sample}"
    benchmark:
        "analysis/rsem/{sample}/{sample}.rsem.benchmark.txt"
    conda: "../envs/rsem_env.yml"
    params:
        sample_name = lambda wildcards: "analysis/rsem/{sample}/{sample}".format(sample=wildcards.sample),
        # stranded = "--strand-specific" if config["stranded"] else "",
        # paired_end = "--paired-end" if len(config["mate"]) == 2 else "",
        gz_support=lambda wildcards: "--star-gzipped-read-file" if config["samples"][wildcards.sample][0][-3:] == '.gz' else ""
    threads:
        2
    shell:
        "rsem-calculate-expression -p {threads} {params.stranded} "
        "{params.paired_end} "
        "--estimate-rspd "
        "--bam --no-bam-output "
        "{input} {config[rsem_index]} "
        "{params.sample_name} > {log}"
        # "rsem-calculate-expression -p {threads} {params.stranded} "
  #       "{params.paired_end} "
  #       "--star {params.gz_support} "
  #       "--estimate-rspd --append-names "
  #       "{input} {config[rsem_index]} "
  #       "{params.sample_name} > {log}"

rule rsem_iso_matrix:
    input:
        rsem_iso_files = expand( "analysis/rsem/{sample}/{sample}.isoforms.results", sample=config["samples"] ),
        # metasheet = config['metasheet']
    output:
        rsem_iso_matrix = "analysis/rsem/rsem_tpm_iso_matrix.csv"
    message: 
        "Running RSEM matrix generation rule for isoforms"
    benchmark:
        "analysis/rsem/rsem_iso_matrix_benchmark.txt"
    params:
        args = lambda wildcards, input: merge_sep_inputs({input.rsem_iso_files})
    conda: "../envs/stat_perl_r.yml"
    shell:
        "perl src/raw_and_fpkm_count_matrix.pl --column 5 --metasheet {input.metasheet} --header -f {params.args} 1>{output.rsem_iso_matrix}"


rule rsem_gene_matrix:
    input:
        rsem_gene_files = expand( "analysis/rsem/{sample}/{sample}.genes.results", sample=config["samples"] ),
        # metasheet = config["metasheet"]
    output:
        rsem_gene_matrix = "analysis/rsem/rsem_tpm_gene_matrix.csv"
    message: 
        "Running RSEM matrix generation rule for genes"
    benchmark:
        "analysis/rsem/rsem_gene_matrix_benchmark.txt"
    params:
        args = lambda wildcards, input: merge_sep_inputs({input.rsem_gene_files})
    conda: "../envs/stat_perl_r.yml"
    shell:
        "perl src/raw_and_fpkm_count_matrix.pl --column 5 --metasheet {input.metasheet} --header -f {params.args} 1>{output.rsem_gene_matrix}"

# rule rsem_filter_gene_ct_matrix:
#   """filters the rsem gene count.
#   REPLACES: preprocess.snakefile- filter_Cuff_matrix"""
#   input:
#   #TODO: handle batch_effect correction
#       tpmFile = "analysis/" + config["token"] + "/rsem/rsem_gene_ct_matrix.csv",
#       annotFile=config['metasheet'],
#       force_run_upon_config_change = config['config_file']
#   output:
#       filtered_tpm = "analysis/" + config["token"] + "/rsem/rsem_gene_ct_matrix.filtered.csv"
#   params:
#       sample_names = " ".join(config["ordered_sample_list"])
#   message: "Generating Pre-processed RSEM TPM matrix file"
#   benchmark:
#       "benchmarks/" + config["token"] + "/rsem_filter_gene_ct_matrix.txt"
#   shell:
#       "Rscript viper/modules/scripts/rsem_filter_gene_ct_matrix.R "
#       "--tpm_file {input.tpmFile} "
#       "--min_samples {config[min_num_samples_expressing_at_threshold]} "
#       "--TPM_cutoff {config[TPM_threshold]} "
#       "--filter_miRNA {config[filter_mirna]} "
#       "--numgenes {config[numgenes_plots]} "
#       "--out_file {output.filtered_tpm} "
#       "--sample_names {params.sample_names} "

# rule batch_effect_removal_rsem:
#   input:
#   gene_cts = "analysis/" + config["token"] + "/rsem/rsem_gene_ct_matrix.csv",
#   annotFile = config["metasheet"]
#   output:
#   rsemcsvoutput="analysis/" + config["token"] + "/rsem/batch_corrected_rsem_gene_ct_matrix.csv",
#   rsempdfoutput="analysis/" + config["token"] + "/rsem/rsem_combat_qc.pdf"
#   params:
#   batch_column="batch",
#   datatype = "rsem"
#   message: "Removing batch effect from Gene Count matrix, if errors, check metasheet for batches, refer to README for specifics"
#   #priority: 2
#   benchmark:
#   "benchmarks/" + config["token"] + "/batch_effect_removal_rsem.txt"
#   shell:
#   "Rscript viper/modules/scripts/batch_effect_removal.R {input.gene_cts} {input.annotFile} {params.batch_column} "
#   "{params.datatype} {output.rsemcsvoutput} {output.rsempdfoutput} "
#   " && mv {input.gene_cts} analysis/{config[token]}/rsemlinks/without_batch_correction_rsem_gene_ct_matrix.csv "

# rule tpm_plot:
#   input:
#   #TODO: handle batch_effect correction
#   tpm_mat = "analysis/" + config['token'] + "/rsem/rsem_gene_ct_matrix.filtered.csv",
#   annotFile = config["metasheet"]
#   output:
#   tpm_png = "analysis/" + config["token"] + "/plots/gene_counts.tpm.png"
#   message: "Plot gene counts at various tpm cutoffs"
#   benchmark:
#   "benchmarks/" + config["token"] + "/rsem_tpm_plot.txt"
#   shell:
#   "Rscript viper/modules/scripts/fpkm_plot.R {input.tpm_mat} {output.tpm_png}"


