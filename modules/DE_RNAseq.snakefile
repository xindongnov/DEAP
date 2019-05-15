# Differential expression

def DE_RNAseq_targets(wildcards):
    ls = []
    for run in config["RS_runs"]:
        if not config["RS_runs"][run]['compare'].empty:
            print(run)
            ls.append("analysis/%s/%s_design_matrix.txt" % (run,run))
            ls.append("analysis/%s/DE/%s_TPM_matrix.txt" % (run,run))
            ls.append("analysis/%s/DE/%s_Rawcount_matrix.txt" % (run,run))
        # for sample in config["RS_runs"][run]['samples']:
        #     ls.append("analysis/%s/DE/%s_TPM_matrix.txt" % (run,run))
        print(ls)
    return ls

def get_quantsf(wildcards):
    r = wildcards.run
    sf =[]
    if config['RS_runs'][r]['samples']:
        for s in config['RS_runs'][r]['samples']:
            sf.append("analysis/%s/samples/%s/align/quant.sf" % (r,s))
    print(sf)
    return sf


rule get_specific_design:
    output:
        "analysis/{run}/{run}_design_matrix.txt"
    run:
        with open(str(output),'w') as op:
            op.write(config["RS_runs"][wildcards.run]['compare'].to_csv(sep=',',index=False))

rule gather_TPM:
    input:
        get_quantsf
    output:
        TPM="analysis/{run}/DE/{run}_TPM_matrix.txt",
        RawCount="analysis/{run}/DE/{run}_Rawcount_matrix.txt"
    params:
        input_dir=lambda wildcards: 'analysis/%s' % wildcards.run,
        species='Human' if config['assembly'] == 'hg38' else 'Mouse'
    shell:
        "Rscript ./DEAP/modules/scripts/RNASeq_get_TPM_Rawcount.r {params.species} {params.input_dir} {output.TPM} {output.RawCount}"












