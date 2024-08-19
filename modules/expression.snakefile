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
            ls.append(f"{RES_PATH}/expression/{run}/{run}_exp_matrix_probe.txt")
            ls.append(f"{RES_PATH}/expression/{run}/{run}_exp_matrix_transcription.txt")
            ls.append(f"{RES_PATH}/expression/{run}/{run}_exp_matrix_gene.txt")
            ls.append(f"{RES_PATH}/expression/{run}/{run}_PCA.png")
            for comp in config["runs"][run]['compare']:
                ctrl = config["runs"][run]['compare'][comp]['control']['name']
                treat = config["runs"][run]['compare'][comp]['treat']['name']
                if len(config["runs"][run]['compare'][comp]['control']['sample']) > 1 and len(config["runs"][run]['compare'][comp]['treat']['sample']) > 1:
                    ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}_exp_matrix_probe.txt")
                    ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}_exp_matrix_transcription.txt")
                    ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}_exp_matrix_gene.txt")
                    ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}_limma_table.txt")
                    ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}_DE.txt")
                    ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}_volca_plot_{config['lfc']}_{config['fdr']}.png")
                    ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}.upRegGenes.txt")
                    ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}.downRegGenes.txt")
        elif config["runs"][run]["type"] == "RS":
            ls.append(f"{RES_PATH}/expression/{run}/{run}_TPM_transcript_matrix.txt")
            ls.append(f"{RES_PATH}/expression/{run}/{run}_Rawcount_transcript_matrix.txt")
            ls.append(f"{RES_PATH}/expression/{run}/{run}_TPM_gene_matrix.txt")
            ls.append(f"{RES_PATH}/expression/{run}/{run}_Rawcount_gene_matrix.txt")
            ls.append(f"{RES_PATH}/expression/{run}/{run}_PCA.png")
            for comp in config["runs"][run]['compare']:
                ctrl = config["runs"][run]['compare'][comp]['control']['name']
                treat = config["runs"][run]['compare'][comp]['treat']['name']
                if len(config["runs"][run]['compare'][comp]['control']['sample']) > 1 and len(config["runs"][run]['compare'][comp]['treat']['sample']) > 1:
            # #         # ls.append("%s/expression/%s/%s_Compare_detail.txt" % (RES_PATH,run,run))
                    ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}_DESeq_table.txt")
                    ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}_DE.txt")
                    for lfc in config['lfc']:
                        for fdr in config['fdr']:
                            ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}_volca_plot_{lfc:.1f}_{fdr}.png")
            #         ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}_volca_plot_{config['lfc']}_{config['fdr']}.png")
                            ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}.upRegGenes_{lfc:.1f}_{fdr}.txt")
                            ls.append(f"{RES_PATH}/expression/{run}/{comp}/{run}_{comp}.downRegGenes_{lfc:.1f}_{fdr}.txt")
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
        return "%s/expression/%s/%s_exp_matrix_gene.txt" % (RES_PATH, wildcards.run, wildcards.run)
    elif config["runs"][wildcards.run]["type"] == "RS":
        return "%s/expression/%s/%s_TPM_gene_matrix.txt" % (RES_PATH, wildcards.run, wildcards.run)

def get_trans_quant(wildcards):
    r = wildcards.run
    lst =[]
    if r in config['runs']:
        if config['runs'][r]['samples']:
            for s in config['runs'][r]['samples']:
                if config['aligner'] == 'salmon':
                    lst.append("%s/salmon/%s/%s.quant.sf" % (RES_PATH,s,s))
                else:
                    lst.append("%s/STAR/%s/%s.ReadsPerTranscript.byfeatureCounts.txt" % (RES_PATH,s,s))
    return lst

def get_gene_quant(wildcards):
    r = wildcards.run
    lst =[]
    if r in config['runs']:
        if config['runs'][r]['samples']:
            for s in config['runs'][r]['samples']:
                if config['aligner'] == 'salmon':
                    lst.append("%s/salmon/%s/%s.quant.genes.sf" % (RES_PATH,s,s))
                else:
                    lst.append("%s/STAR/%s/%s.ReadsPerGene.byfeatureCounts.txt" % (RES_PATH,s,s))
    return lst

def get_format_input(wildcards):
    if config["runs"][wildcards.run]["type"].startswith("MA_"):
        return "%s/expression/%s/%s/%s_%s_limma_table.txt" % (RES_PATH,wildcards.run,wildcards.compare,wildcards.run,wildcards.compare)
    elif config["runs"][wildcards.run]["type"] == "RS":
        return "%s/expression/%s/%s/%s_%s_DESeq_table.txt" % (RES_PATH,wildcards.run,wildcards.compare,wildcards.run,wildcards.compare)

def get_species(wildcards):
    if config['assembly'].startswith("hg"):
        species = 'Human'
    elif config['assembly'].startswith("mm"):
        species = 'Mouse'
    elif config['assembly'].startswith("rn"):
        species = 'Rat'
    return species


rule expression_GetSpecificDesign:
    output:
        "%s/expression/{run}/{run}_design_matrix.txt" % RES_PATH
    message: "Expression: Get specific design for {wildcards.run}"
    run:
        with open(str(output),'w') as design_file:
            design_file.write(config["runs"][wildcards.run]['raw_design'].to_csv(sep=',',index=False))

rule expression_MicroarrayGatherAllExpresion:
    input:
        get_MA_Input
    output:
        probe="analysis/expression/{run}/{run}_exp_matrix_probe.txt",
        trans="analysis/expression/{run}/{run}_exp_matrix_transcription.txt",
        gene="analysis/expression/{run}/{run}_exp_matrix_gene.txt",
    log: "analysis/logs/expression/{run}/{run}_expression_MicroarrayGatherAllExpresion.log"
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
        probe="analysis/expression/{run}/{compare}/{run}_{compare}_exp_matrix_probe.txt",
        trans="analysis/expression/{run}/{compare}/{run}_{compare}_exp_matrix_transcription.txt",
        gene="analysis/expression/{run}/{compare}/{run}_{compare}_exp_matrix_gene.txt",
    log: "analysis/logs/expression/{run}/{run}_expression_MicroarrayGatherExpression_{compare}.log"
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
        design="analysis/expression/{run}/{run}_design_matrix.txt"
    output:
        "analysis/expression/{run}/{run}_PCA.png"
    log: "analysis/logs/expression/{run}/{run}_expression_PCAplot.log"
    message: "Expression: draw PCA plot for {wildcards.run}"
    # params:
    #     protocal = lambda wildcards: "Microarray" if config['runs'][wildcards.run]["type"].startswith("MA_") else "RNASeq"
    shell:
        "Rscript ./DEAP/modules/scripts/draw_PCA_plot.r -i {input.matrix} -d {input.design} -r {output} > {log} 2>&1"


rule expression_MicroarrayDifferentialExpression:
    # only for Microarray
    input:
        "analysis/expression/{run}/{compare}/{run}_{compare}_exp_matrix_probe.txt"
    output:
        table="analysis/expression/{run}/{compare}/{run}_{compare}_limma_table.txt"
        # "analysis/expression/{run}/{run}_Compare_detail.txt"
    log: "analysis/logs/expression/{run}/{run}_expression_MicroarrayDifferentialExpression_{compare}.log"
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

rule expression_RNAseqGatherTranscript:
    # only for RNA-seq
    input:
        get_trans_quant
    output:
        tpm="analysis/expression/{run}/{run}_TPM_transcript_matrix.txt",
        rawcount="analysis/expression/{run}/{run}_Rawcount_transcript_matrix.txt"
    message: "Expression: gathering all transcripts TPM and raw count for RNA-seq experiment {wildcards.run}"
    log: "analysis/logs/expression/{run}/{run}_expression_RNAseqGatherTpm.log"
    params:
        files=lambda wildcards, input: ",".join(input),
        script='expression_RNASeq_get_matrix_from_STAR.py' if config['aligner'] == 'STAR' else 'expression_RNASeq_get_matrix_from_salmon.py',
    shell:
        "python ./DEAP/modules/scripts/{params.script} -i {params.files} -t {output.tpm} -c {output.rawcount} > {log} 2>&1"

rule expression_RNAseqGatherGene:
    # only for RNA-seq
    input:
        get_gene_quant
    output:
        tpm="analysis/expression/{run}/{run}_TPM_gene_matrix.txt",
        rawcount="analysis/expression/{run}/{run}_Rawcount_gene_matrix.txt",
    message: "Expression: gathering all raw count for RNA-seq experiment {wildcards.run}"
    log: "analysis/logs/expression/{run}/{run}_expression_RNAseqGatherRawcount.log"
    params:
        files=lambda wildcards, input: ",".join(input),
        script='expression_RNASeq_get_matrix_from_STAR.py' if config['aligner'] == 'STAR' else 'expression_RNASeq_get_matrix_from_salmon.py'
    shell:
        "python ./DEAP/modules/scripts/{params.script} -i {params.files} -t {output.tpm} -c {output.rawcount} > {log} 2>&1"

# rule expression_RNAseqTranscriptToGene: # should use tximport, but rsem and salmon all have similar function
#     # only for RNA-seq
#     input:
#         transcript_matrix="analysis/expression/{run}/{run}_{type}_transcript_matrix.txt",
#     output:
#         gene_matrix="analysis/expression/{run}/{run}_{type}_gene_matrix.txt",
#     message: "Expression: Transfer gene name from transcript ID of {wildcards.run}"
#     log: "analysis/logs/expression/{run}/{run}_{type}_expression_RNAseqTranscriptToGene.log"
#     params:
#         ref = config['transformed_gtf']
#         # files=lambda wildcards, input: ",".join(input),
#         # species=get_species,
#         # names=lambda wildcards, input: ",".join([i.replace("analysis/%s/samples/" % wildcards.run, "").replace("/align/quant.sf", "") for i in input]),
#     shell:
#         "python ./DEAP/modules/scripts/expression_RNASeq_matrix_transcript_to_gene.py "
#         "-i {input.transcript_matrix} -o {output.gene_matrix} -r {params.ref} > {log} 2>&1"

rule expression_RNAseqDifferentialExpression:
    # For RNA-seq
    input:
        rawcount="analysis/expression/{run}/{run}_Rawcount_gene_matrix.txt",
        design="analysis/expression/{run}/{run}_design_matrix.txt"
    output:
        # detail="analysis/expression/{run}/{run}_Compare_detail.txt",
        table="analysis/expression/{run}/{compare}/{run}_{compare}_DESeq_table.txt"
    log: "analysis/logs/expression/{run}/{run}_expression_RNAseqDifferentialExpression_{compare}.log"
    message: "Expression: call differential expression for RNA-seq experiment {wildcards.run} in {wildcards.compare}"
    # params:
    #     treat=lambda wildcards: ",".join(config["runs"][wildcards.run]["compare"][wildcards.compare]["treat"]["sample"]),
    #     treatname=lambda wildcards: config["runs"][wildcards.run]["compare"][wildcards.compare]["treat"]["name"],
    #     control=lambda wildcards: ",".join(config["runs"][wildcards.run]["compare"][wildcards.compare]["control"]["sample"]),
    #     controlname=lambda wildcards: config["runs"][wildcards.run]["compare"][wildcards.compare]["control"]["name"],
    shell:
        "python DEAP/modules/scripts/expression_RNASeq_DE_analysis.py -c {input.rawcount} "
        "-o {output.table} -d {input.design} -cmp {wildcards.compare} > {log} 2>&1"

rule expression_formatTable:
    input:
        get_format_input
    output:
        "analysis/expression/{run}/{compare}/{run}_{compare}_DE.txt"
    # log: "analysis/logs/expression/{run}/{run}_expression_formatTable.log"
    message: "Expression: format differential expression table for {wildcards.run} in {wildcards.compare}"
    params:
        columns = lambda wildcards: r'{print $1"\t"$2"\t"$3"\t"$6"\t"$7}' if config['runs'][wildcards.run]["type"]=="RS" else r'{print $1"\t"$2"\t"$3"\t"$6"\t"$7}'
    shell:
        """awk -F'\\t' '{params.columns}' {input} > {output}"""


rule expression_volcaPlot:
    # both RNA-seq and microarray
    input:
        "analysis/expression/{run}/{compare}/{run}_{compare}_DE.txt"
    output:
        "analysis/expression/{run}/{compare}/{run}_{compare}_volca_plot_{lfc}_{fdr}.png"
    log: "analysis/logs/expression/{run}/{run}_expression_volcaPlot_{compare}_{lfc}_{fdr}.log"
    message: "Expression: draw volcano plot for {wildcards.run} in {wildcards.compare}"
    params:
        add_gene_name = config['add_gene_name'],
        n_gene = config['n_gene'],
        genename_handle = lambda wildcards: "-g" + ",".join(config['plot_gene_name']) if config['plot_gene_name'] != [] else ''
    shell:
        "python DEAP/modules/scripts/expression_draw_volca_plot.py "
        "-i {input} -o {output} "
        "-a {params.add_gene_name} -n {params.n_gene} {params.genename_handle} "
        "-l {wildcards.lfc} -f {wildcards.fdr} -t {wildcards.compare} > {log} 2>&1"


rule expression_upRegGene:
    input:
        "analysis/expression/{run}/{compare}/{run}_{compare}_DE.txt"
    output:
        "analysis/expression/{run}/{compare}/{run}_{compare}.upRegGenes_{lfc}_{fdr}.txt"
    # log: "analysis/logs/expression/{run}/{run}_expression_upRegGene.log"
    message: "Expression: extract up regulatory genes for {wildcards.run} in {wildcards.compare}"
    # params:
    #     foldchange = 0,
    #     pvalue = 1
    shell:
        """tail -n +2 {input} | awk -F '\\t' '$3>{wildcards.lfc} && $5<{wildcards.fdr} {{print}}'| sort -k3 -nr | cut -f1 > {output}"""

rule expression_downRegGene:
    input:
        "analysis/expression/{run}/{compare}/{run}_{compare}_DE.txt"
    output:
        "analysis/expression/{run}/{compare}/{run}_{compare}.downRegGenes_{lfc}_{fdr}.txt"
    # log: "analysis/logs/expression/{run}/{run}_expression_downRegGene.log"
    message: "Expression: extract down regulatory genes for {wildcards.run} in {wildcards.compare}"
    # params:
    #     foldchange = 0,
    #     pvalue = 1
    shell:
        """tail -n +2 {input} | awk -F '\\t' '$3<-{wildcards.lfc} && $5<{wildcards.fdr} {{print}}'| sort -k3 -n | cut -f1 > {output}"""



