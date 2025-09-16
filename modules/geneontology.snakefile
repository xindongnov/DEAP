
# ================================
# @auther: Xin Dong
# @email: xindong9511@gmail.com
# @date: Sep 2019
# ================================


def geneontology_targets(wildcards):
    # print(snakemake.snakemake)
    ls = []
    for run in config["runs"]:
        for comp in config["runs"][run]['compare']:
            ctrl = config["runs"][run]['compare'][comp]['control']['name']
            treat = config["runs"][run]['compare'][comp]['treat']['name']
            # if len(config["runs"][run]['compare'][comp]['control']['sample']) > 1 and len(config["runs"][run]['compare'][comp]['treat']['sample']) > 1:
            for geneType in ["upRegGenes","downRegGenes"]:
                for lfc in config['lfc']:
                    for fdr in config['fdr']:
                        for ont in ["BP", "MF", "CC"]:
                            ls.append(f"{RES_PATH}/geneontology/{run}/{comp}/{run}_{comp}_{lfc:.1f}_{fdr}_{geneType}_GO_{ont}.pdf")
                            ls.append(f"{RES_PATH}/geneontology/{run}/{comp}/{run}_{comp}_{lfc:.1f}_{fdr}_{geneType}_GO_{ont}.txt")
                        ls.append(f"{RES_PATH}/geneontology/{run}/{comp}/{run}_{comp}_{lfc:.1f}_{fdr}_{geneType}_kegg.pdf")
                        ls.append(f"{RES_PATH}/geneontology/{run}/{comp}/{run}_{comp}_{lfc:.1f}_{fdr}_{geneType}_kegg.txt")
    # print(ls)
    return ls


rule geneontology_run:
    input:
        "%s/expression/{run}/{compare}/{run}_{compare}.{geneType}_{lfc}_{fdr}.txt" % RES_PATH
    output:
        go_fig="%s/geneontology/{run}/{compare}/{run}_{compare}_{lfc}_{fdr}_{geneType}_GO_{ont}.pdf" % RES_PATH,
        go_txt="%s/geneontology/{run}/{compare}/{run}_{compare}_{lfc}_{fdr}_{geneType}_GO_{ont}.txt" % RES_PATH,
    message: "Gene Ontology: call GO for {wildcards.geneType} in {wildcards.run} {wildcards.compare} {wildcards.ont} {wildcards.lfc} {wildcards.fdr}"
    log: "%s/logs/geneontology/{run}/{run}_{compare}_{lfc}_{fdr}_{geneType}_GO_{ont}.log" % RES_PATH
    params:
        species='hs' if config["assembly"] in ['hg38', 'hg19'] else 'mm' if config["assembly"] in ['mm39', 'mm10', 'mm9'] else 'rn',
    shell:
        # "lisa oneshot {params.species} {input} > {output} 2>{log}"
        """
        if [ `wc -l {input} | cut -f 1 -d ' '` -lt 1 ]
        then
            touch {output.go_fig}
            touch {output.go_txt}
        else
            Rscript DEAP/modules/scripts/annotation_go.r -t {wildcards.compare}_{wildcards.geneType} \
            -l {input} -s {params.species} -w 30 \
            -o {wildcards.ont} --go_figure_path {output.go_fig} --go_table_path {output.go_txt} > {log} 2>&1
        fi
        """

rule kegg_run:
    input:
        "%s/expression/{run}/{compare}/{run}_{compare}.{geneType}_{lfc}_{fdr}.txt" % RES_PATH
    output:
        go_fig="%s/geneontology/{run}/{compare}/{run}_{compare}_{lfc}_{fdr}_{geneType}_kegg.pdf" % RES_PATH,
        go_txt="%s/geneontology/{run}/{compare}/{run}_{compare}_{lfc}_{fdr}_{geneType}_kegg.txt" % RES_PATH,
    message: "Gene Ontology: call kegg for {wildcards.geneType} in {wildcards.run} {wildcards.compare} {wildcards.lfc} {wildcards.fdr}"
    log: "%s/logs/geneontology/{run}/{run}_{compare}_{lfc}_{fdr}_{geneType}_kegg.log" % RES_PATH
    params:
        species='hs' if config["assembly"] in ['hg38', 'hg19'] else 'mm' if config["assembly"] in ['mm39', 'mm10', 'mm9'] else 'rn',
    shell:
        # "lisa oneshot {params.species} {input} > {output} 2>{log}"
        """
        if [ `wc -l {input} | cut -f 1 -d ' '` -lt 1 ]
        then
            touch {output.go_fig}
            touch {output.go_txt}
        else
            Rscript DEAP/modules/scripts/annotation_kegg.r -t {wildcards.compare}_{wildcards.geneType} \
            -l {input} -s {params.species} -w 30 \
            --kegg_figure_path {output.go_fig} --kegg_table_path {output.go_txt} > {log} 2>&1
        fi
        """



