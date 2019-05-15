# Differential expression

def DE_RNAseq_targets(wildcards):
    ls = []
    for run in config["RS_runs"]:
        ls.append("analysis/%s/%s_design_matrix.txt" % (run,run))
        # for sample in config["RS_runs"][run]['samples']:
        #     ls.append("analysis/%s/DE/%s_TPM_matrix.txt" % (run,run))
    return ls

def get_quantsf(wildcards):
    r = wildcards.run
    sf =[]
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
        "analysis/{run}/DE/{run}_TPM_matrix.txt"
    params:
        species=""
    shell:
        ""
