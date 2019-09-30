
# ================================
# @auther: Xin Dong
# @email: xindong9511@gmail.com
# @date: Sep 2019
# ================================

# {prefix}.1000.lisa_direct.csv

def lisa_target(wildcards):
    ls = []
    for run in config["MA_runs"]:
        ls.append("analysis/%s/expression/%s_exp_matrix.txt" % (run,run))
        ls.append("analysis/%s/expression/%s_PCA.png" % (run,run))
        for comp in config["MA_runs"][run]['compare']:
            ctrl = config["MA_runs"][run]['compare'][comp]['control']['name']
            treat = config["MA_runs"][run]['compare'][comp]['treat']['name']
            if len(config["MA_runs"][run]['compare'][comp]['control']['sample']) > 1 and len(config["MA_runs"][run]['compare'][comp]['treat']['sample']) > 1:
                ls.append("analysis/%s/expression/%s_results/%s_%s_limma_table.txt" % (run,run,treat,ctrl))
                ls.append("analysis/%s/expression/%s_results/%s_%s_DE.txt" % (run,run,treat,ctrl))
                ls.append("analysis/%s/expression/%s_results/%s_%s_volca_plot.png" % (run,run,treat,ctrl))
                ls.append("analysis/%s/expression/%s_results/%s_%s_upRegGenes.txt" % (run,run,treat,ctrl))
                ls.append("analysis/%s/expression/%s_results/%s_%s_downRegGenes.txt" % (run,run,treat,ctrl))
    for run in config["RS_runs"]:
        # print(run)
        # ls.append("analysis/%s/%s_design_matrix.txt" % (run,run))
        ls.append("analysis/%s/expression/%s_TPM_matrix.txt" % (run,run))
        ls.append("analysis/%s/expression/%s_Rawcount_matrix.txt" % (run,run))
        ls.append("analysis/%s/expression/%s_PCA.png" % (run,run))
        # ls.append("analysis/%s/expression/%s_results" %(run,run))
        # ls.append("analysis/%s/expression/%s_Compare_detail.txt" % (run,run))
        for comp in config["RS_runs"][run]['compare']:
            ctrl = config["RS_runs"][run]['compare'][comp]['control']['name']
            treat = config["RS_runs"][run]['compare'][comp]['treat']['name']
            if len(config["RS_runs"][run]['compare'][comp]['control']['sample']) > 1 and len(config["RS_runs"][run]['compare'][comp]['treat']['sample']) > 1:
                # ls.append("analysis/%s/expression/%s_Compare_detail.txt" % (run,run))
                ls.append("analysis/%s/expression/%s_results/%s_%s_DESeq_table.txt" % (run,run,treat,ctrl))
                ls.append("analysis/%s/expression/%s_results/%s_%s_DE.txt" % (run,run,treat,ctrl))
                ls.append("analysis/%s/expression/%s_results/%s_%s_volca_plot.png" % (run,run,treat,ctrl))
                ls.append("analysis/%s/expression/%s_results/%s_%s_upRegGenes.txt" % (run,run,treat,ctrl))
                ls.append("analysis/%s/expression/%s_results/%s_%s_downRegGenes.txt" % (run,run,treat,ctrl))
    return ls