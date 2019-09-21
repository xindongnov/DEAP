# Differential expression

def DE_MA_O_targets(wildcards):
    ls = []
    for run in config["MA_runs"]:
        # print(run)
        if config["MA_runs"][run]['type'] == 'MA_O':
            ls.append("analysis/%s/expression/%s_exp_matrix_oligo.txt" % (run,run))
            ls.append("analysis/%s/expression/%s_PCA_oligo.png" % (run,run))
            for comp in config["MA_runs"][run]['compare']:
                ctrl = config["MA_runs"][run]['compare'][comp]['control']['name']
                treat = config["MA_runs"][run]['compare'][comp]['treat']['name']
                if len(config["MA_runs"][run]['compare'][comp]['control']['sample']) > 1 and len(config["MA_runs"][run]['compare'][comp]['treat']['sample']) > 1:
                    ls.append("analysis/%s/expression/%s_results/%s_%s_limma_table_oligo.txt" % (run,run,treat,ctrl))
                    ls.append("analysis/%s/expression/%s_results/%s_%s_volca_plot_oligo.png" % (run,run,treat,ctrl))
            # ls.append("" % (run,run))
        # print(ls)
    return ls

def getMA_O_Input(wildcards):
    MA_O_target = []
    if config['MA_runs'][wildcards.run]['type'] == 'MA_O':
        MA_O_target.extend(config['MA_runs'][wildcards.run]['samples'])
    # print(wildcards.run)
    return MA_O_target


rule gather_oligo_expression:
    input:
        getMA_O_Input
    output:
        "analysis/{run}/expression/{run}_exp_matrix_oligo.txt"
    params:
        GPL=lambda wildcards: config['MA_runs'][wildcards.run]['GPL'],
        design=lambda wildcards: config['MA_runs'][wildcards.run]['matrix'],
        species='Human' if config['assembly'] == 'hg38' else 'Mouse'
    shell:
        "Rscript ./DEAP/modules/scripts/get_olign_expression_profile.r {input} {params.GPL} {params.species} {params.design} {output}"

rule Microarray_oligo_PCA_plot:
    input:
        "analysis/{run}/expression/{run}_exp_matrix_oligo.txt"
    output:
        "analysis/{run}/expression/{run}_PCA_oligo.png"
    params:
        design_matrix=lambda wildcards: config['MA_runs'][wildcards.run]['matrix']
    shell:
        "Rscript ./DEAP/modules/scripts/draw_PCA_plot.r {input} Microarray {params.design_matrix} {output}"

rule Microarray_oligo_DE:
    input:
        "analysis/{run}/expression/{run}_exp_matrix_oligo.txt"
    output:
        table="analysis/{run}/expression/{run}_results/{treatment}_limma_table_oligo.txt"
        # "analysis/{run}/expression/{run}_Compare_detail.txt"
    params:
        design_matrix=lambda wildcards: config['MA_runs'][wildcards.run]['matrix'],
        output_dir=lambda wildcards: "analysis/%s/expression/%s_results" % (wildcards.run,wildcards.run),
        compare_detail=lambda wildcards: "analysis/%s/expression/%s_results/%s_compare_detail.txt" % (wildcards.run,wildcards.run,wildcards.treatment),
        foldchange=1,
        p=1,
        platform='oligo'
    shell:
        "Rscript ./DEAP/modules/scripts/Microarry_affy_olign_DEG_ananlysis.r {input} {params.design_matrix} {params.output_dir} {params.compare_detail} {params.foldchange} {params.p} {params.platform}"

rule Microarray_oligo_plot:
    input:
        "analysis/{run}/expression/{run}_results/{treatment}_limma_table_oligo.txt"
    output:
        "analysis/{run}/expression/{run}_results/{treatment}_volca_plot_oligo.png"
    params:
        foldchange=10,
        pvalue=0.01
    shell:
        "Rscript ./DEAP/modules/scripts/draw_volca_plot.r {input} {output} {params.foldchange} {params.pvalue}"










