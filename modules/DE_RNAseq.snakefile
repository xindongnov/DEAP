# Differential expression

def DE_RNAseq_targets(wildcards):
    ls = []
    for run in config["RS_runs"]:
        if config["RS_runs"][run]['compare']:
            # print(run)
            # ls.append("analysis/%s/%s_design_matrix.txt" % (run,run))
            ls.append("analysis/%s/expression/%s_TPM_matrix.txt" % (run,run))
            ls.append("analysis/%s/expression/%s_Rawcount_matrix.txt" % (run,run))
            ls.append("analysis/%s/expression/%s_PCA.png" % (run,run))
            ls.append("analysis/%s/expression/%s_Compare_detail.txt" % (run,run))
            # ls.append("analysis/%s/expression/%s_PCA.png" % (run,run))
        # print(ls)
    return ls

def get_quantsf(wildcards):
    r = wildcards.run
    sf =[]
    if config['RS_runs'][r]['samples']:
        for s in config['RS_runs'][r]['samples']:
            sf.append("analysis/%s/samples/%s/align/quant.sf" % (r,s))
    # print(sf)
    return sf


rule get_specific_design:
    output:
        temp("analysis/{run}/{run}_design_matrix.txt")
    run:
        with open(str(output),'w') as op:
            op.write(config["RS_runs"][wildcards.run]['raw_design'].to_csv(sep=',',index=False))

rule gather_TPM_and_Rawcount:
    input:
        get_quantsf
    output:
        TPM="analysis/{run}/expression/{run}_TPM_matrix.txt",
        RawCount="analysis/{run}/expression/{run}_Rawcount_matrix.txt"
    params:
        input_dir=lambda wildcards: 'analysis/%s' % wildcards.run,
        species='Human' if config['assembly'] == 'hg38' else 'Mouse'
    shell:
        "Rscript ./DEAP/modules/scripts/RNASeq_get_TPM_Rawcount.r {params.species} {params.input_dir} {output.TPM} {output.RawCount}"

rule RNASeq_PCA_plot:
    input:
        "analysis/{run}/expression/{run}_TPM_matrix.txt"
    output:
        "analysis/{run}/expression/{run}_PCA.png"
    params:
        expr_type='RNAseq',
        design_matrix=lambda wildcards: "analysis/%s/%s_design_matrix.txt" % (wildcards.run,wildcards.run)
    shell:
        "Rscript ./DEAP/modules/scripts/draw_PCA_plot.r {input} {params.expr_type} {params.design_matrix} {output}"

rule RNASeq_DE:
    input:
        rawcount="analysis/{run}/expression/{run}_Rawcount_matrix.txt",
        design="analysis/{run}/{run}_design_matrix.txt"
    output:
        "analysis/{run}/expression/{run}_Compare_detail.txt"
    params:
        output_dir=lambda wildcards: "analysis/%s/expression/%s_results" % (wildcards.run,wildcards.run)
    shell:
        "Rscript ./DEAP/modules/scripts/RNASeq_DE_analysis.r {input.rawcount} {input.design} {params.output_dir} {output}"

rule RNAseq_volca_plot:
    input:
        "analysis/{run}/expression/{run}_results/{treatment}_DESeq_table.txt"
    output:
        "analysis/{run}/expression/{run}_results/{treatment}_volca_plot.png"
    params:
        foldchange=100,
        pvalue=0.01
    shell:
        "Rscript ./DEAP/modules/scripts/draw_volca_plot.r {input} {output} {params.foldchange} {params.pvalue}"










