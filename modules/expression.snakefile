# Differential expression

# ================================
# @auther: Xin Dong
# @email: xindong9511@gmail.com
# @date: Sep 2019
# ================================

def experssion_targets(wildcards):
    ls = []
    for run in config["runs"]:
        if config["runs"][run]["type"].startswith("MA_"):
            ls.append("analysis/%s/expression/%s_exp_matrix_probe.txt" % (run,run))
            ls.append("analysis/%s/expression/%s_exp_matrix_transcription.txt" % (run,run))
            ls.append("analysis/%s/expression/%s_exp_matrix_gene.txt" % (run,run))
            # ls.append("analysis/%s/expression/%s_PCA.png" % (run,run))
            for comp in config["runs"][run]['compare']:
                ctrl = config["runs"][run]['compare'][comp]['control']['name']
                treat = config["runs"][run]['compare'][comp]['treat']['name']
                if len(config["runs"][run]['compare'][comp]['control']['sample']) > 1 and len(config["runs"][run]['compare'][comp]['treat']['sample']) > 1:
                    ls.append("analysis/%s/expression/%s/%s_exp_matrix_probe.txt" % (run,comp,comp))
                    ls.append("analysis/%s/expression/%s/%s_exp_matrix_transcription.txt" % (run,comp,comp))
                    ls.append("analysis/%s/expression/%s/%s_exp_matrix_gene.txt" % (run,comp,comp))
                    ls.append("analysis/%s/expression/%s/%s_limma_table.txt" % (run,comp,comp))
            #         ls.append("analysis/%s/expression/%s_results/%s_%s_DE.txt" % (run,run,treat,ctrl))
            #         ls.append("analysis/%s/expression/%s_results/%s_%s_volca_plot.png" % (run,run,treat,ctrl))
            #         ls.append("analysis/%s/expression/%s_results/%s_%s.upRegGenes.txt" % (run,run,treat,ctrl))
            #         ls.append("analysis/%s/expression/%s_results/%s_%s.downRegGenes.txt" % (run,run,treat,ctrl))
        # elif config["runs"][run]["type"] == "RS":
            # ls.append("analysis/%s/expression/%s_TPM_matrix.txt" % (run,run))
            # ls.append("analysis/%s/expression/%s_Rawcount_matrix.txt" % (run,run))
            # ls.append("analysis/%s/expression/%s_PCA.png" % (run,run))
            # ls.append("analysis/%s/expression/%s_results" %(run,run))
            # ls.append("analysis/%s/expression/%s_Compare_detail.txt" % (run,run))
            # for comp in config["runs"][run]['compare']:
            #     ctrl = config["runs"][run]['compare'][comp]['control']['name']
            #     treat = config["runs"][run]['compare'][comp]['treat']['name']
            #     if len(config["runs"][run]['compare'][comp]['control']['sample']) > 1 and len(config["runs"][run]['compare'][comp]['treat']['sample']) > 1:
            #         # ls.append("analysis/%s/expression/%s_Compare_detail.txt" % (run,run))
            #         ls.append("analysis/%s/expression/%s_results/%s_%s_DESeq_table.txt" % (run,run,treat,ctrl))
            #         ls.append("analysis/%s/expression/%s_results/%s_%s_DE.txt" % (run,run,treat,ctrl))
            #         ls.append("analysis/%s/expression/%s_results/%s_%s_volca_plot.png" % (run,run,treat,ctrl))
            #         ls.append("analysis/%s/expression/%s_results/%s_%s.upRegGenes.txt" % (run,run,treat,ctrl))
            #         ls.append("analysis/%s/expression/%s_results/%s_%s.downRegGenes.txt" % (run,run,treat,ctrl))
    return ls

def get_MA_Input(wildcards):
    MA_target = []
    if config["runs"][wildcards.run]["type"].startswith("MA_"):
        MA_target.extend(config['runs'][wildcards.run]['samples'][s][0] for s in config['runs'][wildcards.run]['samples'])
    # print(MA_target)
    return MA_target

def get_MA_Compare_Input(wildcards):
    r = wildcards.run
    c = wildcards.compare
    MA_target = []
    if config["runs"][r]["type"].startswith("MA_"):
        for s in config['runs'][r]['compare'][c]["treat"]["sample"]:
            print(s)
            MA_target.extend(config['runs'][r]['samples'][s])
        for s in config['runs'][r]['compare'][c]["control"]["sample"]:
            MA_target.extend(config['runs'][r]['samples'][s])
    print(MA_target)
    return MA_target

def get_MA_Compare_Name(wildcards):
    r = wildcards.run
    c = wildcards.compare
    names = []
    if config["runs"][r]["type"].startswith("MA_"):
        names.extend(config['runs'][r]['compare'][c]["treat"]["sample"])
        names.extend(config['runs'][r]['compare'][c]["control"]["sample"])
    names_str = ",".join(names)
    return names_str

def expression_scripts_switcher(wildcards):
    if config['runs'][wildcards.run]['type'] == 'MA_A':
        return "expression_affy_exprs.r"
    elif config['runs'][wildcards.run]['type'] == 'MA_O':
        return "expression_oligo_exprs.r"
    else:
        return ""

def DE_scripts_switcher(wildcards):
    if config['runs'][wildcards.run]['type'] == 'MA_A':
        return "expression_affy_limma.r"
    elif config['runs'][wildcards.run]['type'] == 'MA_O':
        return "expression_oligo_limma.r"
    else:
        return ""

rule expression_MicroarrayGatherAllExpresion:
    input:
        get_MA_Input
    output:
        probe="analysis/{run}/expression/{run}_exp_matrix_probe.txt",
        trans="analysis/{run}/expression/{run}_exp_matrix_transcription.txt",
        gene="analysis/{run}/expression/{run}_exp_matrix_gene.txt",
    params:
        GPL=lambda wildcards: "%s/%s.txt" % (config['GPL'],config['runs'][wildcards.run]['platform']),
        scripts=lambda wildcards: expression_scripts_switcher(wildcards),
        inputs=lambda wildcards, input: ",".join(input),
        names=lambda wildcards: ",".join([s for s in config['runs'][wildcards.run]['samples']])
    shell:
        "Rscript ./DEAP/modules/scripts/{params.scripts} -f \"{params.inputs}\" -n {params.names} -l {params.GPL} -p {output.probe} -t {output.trans} -g {output.gene}"

rule expression_MicroarrayGatherExpression:
    input:
        get_MA_Compare_Input
    output:
        probe="analysis/{run}/expression/{compare}/{compare}_exp_matrix_probe.txt",
        trans="analysis/{run}/expression/{compare}/{compare}_exp_matrix_transcription.txt",
        gene="analysis/{run}/expression/{compare}/{compare}_exp_matrix_gene.txt",
    params:
        GPL=lambda wildcards: "%s/%s.txt" % (config['GPL'],config['runs'][wildcards.run]['platform']),
        scripts=lambda wildcards: expression_scripts_switcher(wildcards),
        inputs=lambda wildcards, input: ",".join(input),
        names=lambda wildcards: get_MA_Compare_Name(wildcards)
    shell:
        "Rscript ./DEAP/modules/scripts/{params.scripts} -f \"{params.inputs}\" -n {params.names} -l {params.GPL} -p {output.probe} -t {output.trans} -g {output.gene}"

rule expression_MicroarrayDifferentialExpression:
    # only for Microarray
    input:
        "analysis/{run}/expression/{compare}/{compare}_exp_matrix_probe.txt"
    output:
        table="analysis/{run}/expression/{compare}/{compare}_limma_table.txt"
        # "analysis/{run}/expression/{run}_Compare_detail.txt"
    params:
        scripts=lambda wildcards: DE_scripts_switcher(wildcards),
        GPL=lambda wildcards: "%s/%s.txt" % (config['GPL'],config['runs'][wildcards.run]['platform']),
        foldchange=0,
        p=1,
        treatsample=lambda wildcards: ",".join(config["runs"][wildcards.run]["compare"][wildcards.compare]["treat"]["sample"]), 
        contsample=lambda wildcards: ",".join(config["runs"][wildcards.run]["compare"][wildcards.compare]["control"]["sample"]),
        treatname=lambda wildcards: config["runs"][wildcards.run]["compare"][wildcards.compare]["treat"]["name"],
        contname=lambda wildcards: config["runs"][wildcards.run]["compare"][wildcards.compare]["control"]["name"],
    shell:
        "Rscript ./DEAP/modules/scripts/{params.scripts} -i {input} -l {params.GPL} -r {output.table} -f {params.foldchange} -q {params.p} -t {params.treatsample} "
        "-c {params.contsample} --treatname {params.treatname} --controlname {params.contname}"















# def get_quantsf(wildcards):
#     r = wildcards.run
#     sf =[]
#     if r in config['runs']:
#         if config['runs'][r]['samples']:
#             for s in config['runs'][r]['samples']:
#                 sf.append("analysis/%s/samples/%s/align/quant.sf" % (r,s))
#     # print(sf)
#     return sf


# def get_PCAplot_Input(wildcards):
#     if config["runs"][wildcards.run]["type"].startswith("MA_"):
#         return "analysis/%s/expression/%s_exp_matrix.txt" % (wildcards.run,wildcards.run)
#     elif config["runs"][wildcards.run]["type"] == "RS":
#         return ["analysis/%s/expression/%s_TPM_matrix.txt" % (wildcards.run,wildcards.run), "analysis/%s/%s_design_matrix.txt" % (wildcards.run,wildcards.run)]

# def get_format_input(wildcards):
#     if config["runs"][wildcards.run]["type"].startswith("MA_"):
#         return "analysis/%s/expression/%s_results/%s_limma_table.txt" % (wildcards.run,wildcards.run,wildcards.treatment)
#     elif config["runs"][wildcards.run]["type"] == "RS":
#         return "analysis/%s/expression/%s_results/%s_DESeq_table.txt" % (wildcards.run,wildcards.run,wildcards.treatment)


# rule expression_RNAseqGetSpecificDesign:
#     # only for RNA-seq
#     output:
#         temp("analysis/{run}/{run}_design_matrix.txt")
#     message: "Expression: Get specific design for {wildcards.run}"
#     run:
#         with open(str(output),'w') as op:
#             op.write(config["runs"][wildcards.run]['raw_design'].to_csv(sep=',',index=False))


# rule expression_RNAseqGatherTpmAndRawcount:
#     # only for RNA-seq
#     input:
#         get_quantsf
#     output:
#         TPM="analysis/{run}/expression/{run}_TPM_matrix.txt",
#         RawCount="analysis/{run}/expression/{run}_Rawcount_matrix.txt"
#     message: "Expression: gathering all TMP and raw count for {wildcards.run}"
#     params:
#         input_dir=lambda wildcards: 'analysis/%s' % wildcards.run,
#         species='Human' if config['assembly'] == 'hg38' else 'Mouse'
#     shell:
#         "Rscript ./DEAP/modules/scripts/RNASeq_get_TPM_Rawcount.r {params.species} {params.input_dir} {output.TPM} {output.RawCount}"

# rule expression_MicroarrayGatherExpresion:
#     # only for Microarray
#     input:
#         get_MA_Input
#     output:
#         "analysis/{run}/expression/{run}_exp_matrix.txt"
#     params:
#         GPL=lambda wildcards: config['runs'][wildcards.run]['platform'],
#         design=lambda wildcards: config['MA_runs'][wildcards.run]['matrix'],
#         scripts=lambda wildcards: expression_scripts_switcher(wildcards),
#         species= lambda wildcards: '' if config['MA_runs'][wildcards.run]['type'] == 'MA_A' else ("Mouse" if config["assembly"] == "mm10" else "Human")
#     shell:
#         "Rscript ./DEAP/modules/scripts/{params.scripts} {input} {params.GPL} {params.species} {params.design} {output}"

# rule expression_PCAplot:
#     # Both for RNA-seq and Microarray
#     input:
#         get_PCAplot_Input
#     output:
#         "analysis/{run}/expression/{run}_PCA.png"
#     params:
#         design_matrix=lambda wildcards, input: config['MA_runs'][wildcards.run]['matrix'] if wildcards.run in config['MA_runs'] else input[1],
#         protocal = lambda wildcards: "Microarray" if wildcards.run in config['MA_runs'] else "RNASeq"
#     shell:
#         "Rscript ./DEAP/modules/scripts/draw_PCA_plot.r {input[0]} {params.protocal} {params.design_matrix} {output}"

# rule expression_RNAseqDifferentialExpression:
#     # For RNA-seq
#     input:
#         rawcount="analysis/{run}/expression/{run}_Rawcount_matrix.txt",
#         design="analysis/{run}/{run}_design_matrix.txt"
#     output:
#         # detail="analysis/{run}/expression/{run}_Compare_detail.txt",
#         table="analysis/{run}/expression/{run}_results/{treatment}_DESeq_table.txt"
#     params:
#         output_dir=lambda wildcards: "analysis/%s/expression/%s_results" % (wildcards.run,wildcards.run),
#         detail=lambda wildcards: "analysis/%s/expression/%s_results/%s_compare_detail.txt" % (wildcards.run,wildcards.run,wildcards.treatment),
#     shell:
#         "Rscript ./DEAP/modules/scripts/RNASeq_DE_analysis.r {input.rawcount} {input.design} {params.output_dir} {params.detail}"

# rule expression_MicroarrayDifferentialExpression:
#     # only for Microarray
#     input:
#         "analysis/{run}/expression/{run}_exp_matrix.txt"
#     output:
#         table="analysis/{run}/expression/{run}_results/{treatment}_limma_table.txt"
#         # "analysis/{run}/expression/{run}_Compare_detail.txt"
#     params:
#         design_matrix=lambda wildcards: config['MA_runs'][wildcards.run]['matrix'],
#         output_dir=lambda wildcards: "analysis/%s/expression/%s_results" % (wildcards.run,wildcards.run),
#         compare_detail=lambda wildcards: "analysis/%s/expression/%s_results/%s_compare_detail.txt" % (wildcards.run,wildcards.run,wildcards.treatment),
#         foldchange=1,
#         p=1,
#         platform=lambda wildcards: 'affy' if config['MA_runs'][wildcards.run]['type'] == 'MA_A' else "oligo",
#     shell:
#         "Rscript ./DEAP/modules/scripts/Microarry_affy_oligo_DEG_ananlysis.r {input} {params.design_matrix} {params.output_dir} {params.compare_detail} {params.foldchange} {params.p} {params.platform}"


# rule expression_formatTable:
#     input:
#         get_format_input
#     output:
#         "analysis/{run}/expression/{run}_results/{treatment}_DE.txt"
#     params:
#         columns = lambda wildcards: r'{print $2"\t"$3"\t"$4"\t"$7"\t"$8}' if wildcards.run in config['RS_runs'] else r'{print $1"\t"$3"\t"$2"\t"$5"\t"$6}'
#     shell:
#         """awk -F'\\t' '{params.columns}' {input} > {output}"""


# rule expression_volcaPlot:
#     # both RNA-seq and microarray
#     input:
#         "analysis/{run}/expression/{run}_results/{treatment}_DE.txt"
#     output:
#         "analysis/{run}/expression/{run}_results/{treatment}_volca_plot.png"
#     params:
#         foldchange=100,
#         pvalue=0.01
#     shell:
#         "Rscript ./DEAP/modules/scripts/draw_volca_plot.r {input} {output} {params.foldchange} {params.pvalue}"

# rule expression_upRegGene:
#     input:
#         "analysis/{run}/expression/{run}_results/{treatment}_DE.txt"
#     output:
#         "analysis/{run}/expression/{run}_results/{treatment}.upRegGenes.txt"
#     # log: "analysis/{run}/log/{treatment}_upRegGenes.log"
#     params:
#         foldchange = 0,
#         pvalue = 0.05
#     shell:
#         """tail -n +2 {input} | awk -F '\\t' '$3>{params.foldchange} && $5<{params.pvalue} {{print}}'| sort -g -k5 | cut -f1 > {output}"""

# rule expression_downRegGene:
#     input:
#         "analysis/{run}/expression/{run}_results/{treatment}_DE.txt"
#     output:
#         "analysis/{run}/expression/{run}_results/{treatment}.downRegGenes.txt"
#     # log: "analysis/{run}/log/{treatment}_downRegGenes.log"
#     params:
#         foldchange = 0,
#         pvalue = 0.05
#     shell:
#         """tail -n +2 {input} | awk -F '\\t' '$3<{params.foldchange} && $5<{params.pvalue} {{print}}'| sort -g -k5 | cut -f1 > {output}"""









