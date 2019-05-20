# Differential expression

def DE_MA_A_targets(wildcards):
    ls = []
    for run in config["MA_runs"]:
        # print(run)
        if config["MA_runs"][run]['type'] == 'MA_A':
            ls.append("analysis/%s/expression/%s_exp_matrix_affy.txt" % (run,run))
            ls.append("analysis/%s/expression/%s_PCA_affy.png" % (run,run))
            for comp in config["MA_runs"][run]['compare']:
                ctrl = config["MA_runs"][run]['compare'][comp]['control']['name']
                treat = config["MA_runs"][run]['compare'][comp]['treat']['name']
                if len(config["MA_runs"][run]['compare'][comp]['control']['sample']) > 1 and len(config["MA_runs"][run]['compare'][comp]['treat']['sample']) > 1:
                    ls.append("analysis/%s/expression/%s_results/%s_%s_limma_table_affy.txt" % (run,run,treat,ctrl))
                    ls.append("analysis/%s/expression/%s_results/%s_%s_volca_plot_affy.png" % (run,run,treat,ctrl))
            # ls.append("" % (run,run))
        # print(ls)
    return ls

def getMAInput(wildcards):
    # print(wildcards.run)
    return config['MA_runs'][wildcards.run]['samples']


rule gather_affy_expresion:
    input:
        getMAInput
    output:
        "analysis/{run}/expression/{run}_exp_matrix_affy.txt"
    params:
        GPL=lambda wildcards: config['MA_runs'][wildcards.run]['GPL'],
        design=lambda wildcards: config['MA_runs'][wildcards.run]['matrix']
    shell:
        "Rscript ./DEAP/modules/scripts/get_affy_expression_profile.r {input} {params.GPL} {params.design} {output}"

rule Microarray_affy_PCA_plot:
    input:
        "analysis/{run}/expression/{run}_exp_matrix_affy.txt"
    output:
        "analysis/{run}/expression/{run}_PCA_affy.png"
    params:
        design_matrix=lambda wildcards: config['MA_runs'][wildcards.run]['matrix']
    shell:
        "Rscript ./DEAP/modules/scripts/draw_PCA_plot.r {input} Microarray {params.design_matrix} {output}"

rule Microarray_affy_DE:
    input:
        "analysis/{run}/expression/{run}_exp_matrix_affy.txt"
    output:
        table="analysis/{run}/expression/{run}_results/{treatment}_limma_table_affy.txt"
        # "analysis/{run}/expression/{run}_Compare_detail.txt"
    params:
        design_matrix=lambda wildcards: config['MA_runs'][wildcards.run]['matrix'],
        output_dir=lambda wildcards: "analysis/%s/expression/%s_results" % (wildcards.run,wildcards.run),
        foldchange=10,
        p=0.01
    shell:
        "Rscript ./DEAP/modules/scripts/Microarry_affy_olign_DEG_ananlysis.r {input} {params.design_matrix} {params.output_dir} {output.table} {params.foldchange} {params.p}"

rule Microarray_affy_plot:
    input:
        "analysis/{run}/expression/{run}_results/{treatment}_limma_table_affy.txt"
    output:
        "analysis/{run}/expression/{run}_results/{treatment}_volca_plot_affy.png"
    params:
        foldchange=10,
        pvalue=0.01
    shell:
        "Rscript ./DEAP/modules/scripts/draw_volca_plot.r {input} {output} {params.foldchange} {params.pvalue}"










