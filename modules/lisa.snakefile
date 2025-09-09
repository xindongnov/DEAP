
# ================================
# @auther: Xin Dong
# @email: xindong9511@gmail.com
# @date: Sep 2019
# ================================


def lisa_targets(wildcards):
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
                        ls.append(f"{RES_PATH}/lisa/{run}/{comp}/{run}_{comp}_{lfc:.1f}_{fdr}_{geneType}.txt")
                        ls.append(f"{RES_PATH}/lisa/{run}/{comp}/{run}_{comp}_{lfc:.1f}_{fdr}_{geneType}.results.txt")
    # print(ls)
    return ls


rule lisa_copy_file:
    input:
        "analysis/expression/{run}/{compare}/{run}_{compare}.{geneType}_{lfc}_{fdr}.txt"
    output:
        "analysis/lisa/{run}/{compare}/{run}_{compare}_{lfc}_{fdr}_{geneType}.txt"
    shell:
        "head -n 500 {input} > {output}"

rule lisa_run:
    input:
        "analysis/lisa/{run}/{compare}/{run}_{compare}_{lfc}_{fdr}_{geneType}.txt"
    output:
        "analysis/lisa/{run}/{compare}/{run}_{compare}_{lfc}_{fdr}_{geneType}.results.txt"
    message: "LISA: call lisa for {wildcards.geneType} in {wildcards.run} {wildcards.compare}"
    log: "analysis/logs/lisa/{run}/{run}_lisa_{compare}_{geneType}_{lfc}_{fdr}.log"
    params:
        species='mm10' if config["assembly"].startswith('rn') else config["assembly"],
    # threads: 8
    shell:
        # "lisa oneshot {params.species} {input} > {output} 2>{log}"
        """
        if [ `wc -l {input} | cut -f 1 -d ' '` -lt 50 ]
        then
            touch {output}
        else
            lisa oneshot {params.species} {input} --seed=42 > {output} 2>{log}
        fi
        """





