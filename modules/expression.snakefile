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
            if config["check_compare"] == True and config["runs"][run]["samples"] != {}:
                ls.append("analysis/%s/expression/%s_exp_matrix_probe.txt" % (run,run))
                ls.append("analysis/%s/expression/%s_exp_matrix_transcription.txt" % (run,run))
                ls.append("analysis/%s/expression/%s_exp_matrix_gene.txt" % (run,run))
                ls.append("analysis/%s/expression/%s_PCA.png" % (run,run))
            elif config["check_compare"] == False:
                ls.append("analysis/%s/expression/%s_exp_matrix_probe.txt" % (run,run))
                ls.append("analysis/%s/expression/%s_exp_matrix_transcription.txt" % (run,run))
                ls.append("analysis/%s/expression/%s_exp_matrix_gene.txt" % (run,run))
                ls.append("analysis/%s/expression/%s_PCA.png" % (run,run))
            for comp in config["runs"][run]['compare']:
                ctrl = config["runs"][run]['compare'][comp]['control']['name']
                treat = config["runs"][run]['compare'][comp]['treat']['name']
                if len(config["runs"][run]['compare'][comp]['control']['sample']) > 1 and len(config["runs"][run]['compare'][comp]['treat']['sample']) > 1:
                    ls.append("analysis/%s/expression/%s/%s_%s_exp_matrix_probe.txt" % (run,comp,run,comp))
                    ls.append("analysis/%s/expression/%s/%s_%s_exp_matrix_transcription.txt" % (run,comp,run,comp))
                    ls.append("analysis/%s/expression/%s/%s_%s_exp_matrix_gene.txt" % (run,comp,run,comp))
                    ls.append("analysis/%s/expression/%s/%s_%s_limma_table.txt" % (run,comp,run,comp))
                    ls.append("analysis/%s/expression/%s/%s_%s_DE.txt" % (run,comp,run,comp))
                    ls.append("analysis/%s/expression/%s/%s_%s_volca_plot.png" % (run,comp,run,comp))
                    ls.append("analysis/%s/expression/%s/%s_%s.upRegGenes.txt" % (run,comp,run,comp))
                    ls.append("analysis/%s/expression/%s/%s_%s.downRegGenes.txt" % (run,comp,run,comp))
        elif config["runs"][run]["type"] == "RS":
            if config["check_compare"] == True and config["runs"][run]["samples"] != {}:
                ls.append("analysis/%s/expression/%s_TPM_matrix.txt" % (run,run))
                ls.append("analysis/%s/expression/%s_Rawcount_matrix.txt" % (run,run))
                ls.append("analysis/%s/expression/%s_PCA.png" % (run,run))
            elif config["check_compare"] == False:
                ls.append("analysis/%s/expression/%s_TPM_matrix.txt" % (run,run))
                ls.append("analysis/%s/expression/%s_Rawcount_matrix.txt" % (run,run))
                ls.append("analysis/%s/expression/%s_PCA.png" % (run,run))
            for comp in config["runs"][run]['compare']:
                ctrl = config["runs"][run]['compare'][comp]['control']['name']
                treat = config["runs"][run]['compare'][comp]['treat']['name']
                if len(config["runs"][run]['compare'][comp]['control']['sample']) > 1 and len(config["runs"][run]['compare'][comp]['treat']['sample']) > 1:
            #         # ls.append("analysis/%s/expression/%s_Compare_detail.txt" % (run,run))
                    ls.append("analysis/%s/expression/%s/%s_%s_DESeq_table.txt" % (run,comp,run,comp))
                    ls.append("analysis/%s/expression/%s/%s_%s_DE.txt" % (run,comp,run,comp))
                    ls.append("analysis/%s/expression/%s/%s_%s_volca_plot.png" % (run,comp,run,comp))
                    ls.append("analysis/%s/expression/%s/%s_%s.upRegGenes.txt" % (run,comp,run,comp))
                    ls.append("analysis/%s/expression/%s/%s_%s.downRegGenes.txt" % (run,comp,run,comp))
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
            # print(s)
            MA_target.extend(config['runs'][r]['samples'][s])
        for s in config['runs'][r]['compare'][c]["control"]["sample"]:
            MA_target.extend(config['runs'][r]['samples'][s])
    # print(MA_target)
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

def get_PCAplot_Input(wildcards):
    if config["runs"][wildcards.run]["type"].startswith("MA_"):
        return "analysis/%s/expression/%s_exp_matrix_gene.txt" % (wildcards.run,wildcards.run)
    elif config["runs"][wildcards.run]["type"] == "RS":
        return "analysis/%s/expression/%s_TPM_matrix.txt" % (wildcards.run,wildcards.run)

def get_quantsf(wildcards):
    r = wildcards.run
    sf =[]
    if r in config['runs']:
        if config['runs'][r]['samples']:
            for s in config['runs'][r]['samples']:
                sf.append("analysis/%s/samples/%s/align/quant.sf" % (r,s))
    # print(sf)
    return sf

def get_format_input(wildcards):
    if config["runs"][wildcards.run]["type"].startswith("MA_"):
        return "analysis/%s/expression/%s/%s_%s_limma_table.txt" % (wildcards.run,wildcards.compare,wildcards.run,wildcards.compare)
    elif config["runs"][wildcards.run]["type"] == "RS":
        return "analysis/%s/expression/%s/%s_%s_DESeq_table.txt" % (wildcards.run,wildcards.compare,wildcards.run,wildcards.compare)


rule expression_GetSpecificDesign:
    output:
        temp("analysis/{run}/expression/{run}_design_matrix.txt")
    message: "Expression: Get specific design for {wildcards.run}"
    run:
        with open(str(output),'w') as op:
            op.write(config["runs"][wildcards.run]['raw_design'].to_csv(sep=',',index=False))

rule expression_MicroarrayGatherAllExpresion:
    input:
        get_MA_Input
    output:
        probe="analysis/{run}/expression/{run}_exp_matrix_probe.txt",
        trans="analysis/{run}/expression/{run}_exp_matrix_transcription.txt",
        gene="analysis/{run}/expression/{run}_exp_matrix_gene.txt",
    log: "analysis/{run}/log/{run}_expression_MicroarrayGatherAllExpresion.log"
    message: "Expression: gathering all expression profile for microarray experiment {wildcards.run}"
    params:
        GPL=lambda wildcards: "%s/%s.txt" % (config['GPL'],config['runs'][wildcards.run]['platform']),
        scripts=lambda wildcards: expression_scripts_switcher(wildcards),
        inputs=lambda wildcards, input: ",".join(input),
        names=lambda wildcards: ",".join([s for s in config['runs'][wildcards.run]['samples']])
    shell:
        "Rscript ./DEAP/modules/scripts/{params.scripts} -f \"{params.inputs}\" -n {params.names} "
        "-l {params.GPL} -p {output.probe} -t {output.trans} -g {output.gene} > {log} 2>&1"

rule expression_MicroarrayGatherExpression:
    input:
        get_MA_Compare_Input
    output:
        probe="analysis/{run}/expression/{compare}/{run}_{compare}_exp_matrix_probe.txt",
        trans="analysis/{run}/expression/{compare}/{run}_{compare}_exp_matrix_transcription.txt",
        gene="analysis/{run}/expression/{compare}/{run}_{compare}_exp_matrix_gene.txt",
    log: "analysis/{run}/log/{run}_expression_MicroarrayGatherExpression.log"
    message: "Expression: gathering specific expression profile for {wildcards.run} in {wildcards.compare}"
    params:
        GPL=lambda wildcards: "%s/%s.txt" % (config['GPL'],config['runs'][wildcards.run]['platform']),
        scripts=lambda wildcards: expression_scripts_switcher(wildcards),
        inputs=lambda wildcards, input: ",".join(input),
        names=lambda wildcards: get_MA_Compare_Name(wildcards)
    shell:
        "Rscript ./DEAP/modules/scripts/{params.scripts} -f \"{params.inputs}\" -n {params.names} "
        "-l {params.GPL} -p {output.probe} -t {output.trans} -g {output.gene} > {log} 2>&1"


rule expression_PCAplot:
    # Both for RNA-seq and Microarray
    input:
        matrix=get_PCAplot_Input,
        design="analysis/{run}/expression/{run}_design_matrix.txt"
    output:
        "analysis/{run}/expression/{run}_PCA.png"
    log: "analysis/{run}/log/{run}_expression_PCAplot.log"
    message: "Expression: draw PCA plot for {wildcards.run}"
    params:
        protocal = lambda wildcards: "Microarray" if config['runs'][wildcards.run]["type"].startswith("MA_") else "RNASeq"
    shell:
        "Rscript ./DEAP/modules/scripts/draw_PCA_plot.r -i {input.matrix} -d {input.design} "
        "-r {output} -t {params.protocal} > {log} 2>&1"


rule expression_MicroarrayDifferentialExpression:
    # only for Microarray
    input:
        "analysis/{run}/expression/{compare}/{run}_{compare}_exp_matrix_probe.txt"
    output:
        table="analysis/{run}/expression/{compare}/{run}_{compare}_limma_table.txt"
        # "analysis/{run}/expression/{run}_Compare_detail.txt"
    log: "analysis/{run}/log/{run}_expression_MicroarrayDifferentialExpression.log"
    message: "Expression: call differential expression for microarray experiment {wildcards.run} in {wildcards.compare}"
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
        "Rscript ./DEAP/modules/scripts/{params.scripts} -i {input} -l {params.GPL} -r {output.table} "
        "-f {params.foldchange} -q {params.p} -t {params.treatsample} "
        "-c {params.contsample} --treatname {params.treatname} --controlname {params.contname} > {log} 2>&1"

rule expression_RNAseqGatherTpm:
    # only for RNA-seq
    input:
        get_quantsf
    output:
        TPM="analysis/{run}/expression/{run}_TPM_matrix.txt",
        # RawCount="analysis/{run}/expression/{run}_Rawcount_matrix.txt"
    message: "Expression: gathering all TPM for RNA-seq experiment {wildcards.run}"
    log: "analysis/{run}/log/{run}_expression_RNAseqGatherTpm.log"
    params:
        files=lambda wildcards, input: ",".join(input),
        species='Human' if config['assembly'].startswith("hg") else 'Mouse',
        names=lambda wildcards, input: ",".join([i.replace("analysis/%s/samples/" % wildcards.run, "").replace("/align/quant.sf", "") for i in input]),
    shell:
        "Rscript ./DEAP/modules/scripts/expression_RNASeq_get_TPM.r -s {params.species} "
        "-i {params.files} -o {output.TPM} -n {params.names} > {log} 2>&1"

rule expression_RNAseqGatherRawcount:
    # only for RNA-seq
    input:
        get_quantsf
    output:
        rawcount="analysis/{run}/expression/{run}_Rawcount_matrix.txt",
        # RawCount="analysis/{run}/expression/{run}_Rawcount_matrix.txt"
    message: "Expression: gathering all raw count for RNA-seq experiment {wildcards.run}"
    log: "analysis/{run}/log/{run}_expression_RNAseqGatherRawcount.log"
    params:
        files=lambda wildcards, input: ",".join(input),
        species='Human' if config['assembly'].startswith("hg") else 'Mouse',
        names=lambda wildcards, input: ",".join([i.replace("analysis/%s/samples/" % wildcards.run, "").replace("/align/quant.sf", "") for i in input]),
    shell:
        "Rscript ./DEAP/modules/scripts/expression_RNASeq_get_Rawcount.r -s {params.species} "
        "-i {params.files} -o {output.rawcount} -n {params.names} > {log} 2>&1"

rule expression_RNAseqDifferentialExpression:
    # For RNA-seq
    input:
        rawcount="analysis/{run}/expression/{run}_Rawcount_matrix.txt",
    output:
        # detail="analysis/{run}/expression/{run}_Compare_detail.txt",
        table="analysis/{run}/expression/{compare}/{run}_{compare}_DESeq_table.txt"
    log: "analysis/{run}/log/{run}_expression_RNAseqDifferentialExpression.log"
    message: "Expression: call differential expression for RNA-seq experiment {wildcards.run} in {wildcards.compare}"
    params:
        treat=lambda wildcards: ",".join(config["runs"][wildcards.run]["compare"][wildcards.compare]["treat"]["sample"]),
        treatname=lambda wildcards: config["runs"][wildcards.run]["compare"][wildcards.compare]["treat"]["name"],
        control=lambda wildcards: ",".join(config["runs"][wildcards.run]["compare"][wildcards.compare]["control"]["sample"]),
        controlname=lambda wildcards: config["runs"][wildcards.run]["compare"][wildcards.compare]["control"]["name"],
    shell:
        "Rscript ./DEAP/modules/scripts/expression_RNASeq_DE_analysis.r -i {input.rawcount} "
        "-r {output.table} -t {params.treat} -c {params.control} "
        "--treatname {params.treatname} --controlname {params.controlname} > {log} 2>&1"

rule expression_formatTable:
    input:
        get_format_input
    output:
        "analysis/{run}/expression/{compare}/{run}_{compare}_DE.txt"
    # log: "analysis/{run}/log/{run}_expression_formatTable.log"
    message: "Expression: format differential expression table for {wildcards.run} in {wildcards.compare}"
    params:
        columns = lambda wildcards: r'{print $2"\t"$3"\t"$4"\t"$7"\t"$8}' if config['runs'][wildcards.run]["type"]=="RS" else r'{print $1"\t"$3"\t"$2"\t"$5"\t"$6}'
    shell:
        """awk -F'\\t' '{params.columns}' {input} > {output}"""


rule expression_volcaPlot:
    # both RNA-seq and microarray
    input:
        "analysis/{run}/expression/{compare}/{run}_{compare}_DE.txt"
    output:
        "analysis/{run}/expression/{compare}/{run}_{compare}_volca_plot.png"
    log: "analysis/{run}/log/{run}_expression_volcaPlot.log"
    message: "Expression: draw volcano plot for {wildcards.run} in {wildcards.compare}"
    params:
        foldchange=100,
        pvalue=0.01
    shell:
        "Rscript ./DEAP/modules/scripts/draw_volca_plot.r -i {input} -r {output} "
        "-c {params.foldchange} -f {params.pvalue} > {log} 2>&1"

rule expression_upRegGene:
    input:
        "analysis/{run}/expression/{compare}/{run}_{compare}_DE.txt"
    output:
        "analysis/{run}/expression/{compare}/{run}_{compare}.upRegGenes.txt"
    # log: "analysis/{run}/log/{run}_expression_upRegGene.log"
    message: "Expression: extract up regulatory genes for {wildcards.run} in {wildcards.compare}"
    params:
        foldchange = 0,
        pvalue = 1
    shell:
        """tail -n +2 {input} | awk -F '\\t' '$3>{params.foldchange} && $5<{params.pvalue} {{print}}'| sort -g -k5 | cut -f1 > {output}"""

rule expression_downRegGene:
    input:
        "analysis/{run}/expression/{compare}/{run}_{compare}_DE.txt"
    output:
        "analysis/{run}/expression/{compare}/{run}_{compare}.downRegGenes.txt"
    # log: "analysis/{run}/log/{run}_expression_downRegGene.log"
    message: "Expression: extract down regulatory genes for {wildcards.run} in {wildcards.compare}"
    params:
        foldchange = 0,
        pvalue = 1
    shell:
        """tail -n +2 {input} | awk -F '\\t' '$3<{params.foldchange} && $5<{params.pvalue} {{print}}'| sort -g -k5 | cut -f1 > {output}"""



