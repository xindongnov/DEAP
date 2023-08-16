# align salmon 

# ================================
# @auther: Xin Dong
# @email: xindong9511@gmail.com
# @date: Sep 2019
# ================================

salmon_threads = config['aligner_threads']
# global RES_PATH
# RES_PATH = config['res_path']

def getAlignFastq(wildcards):
    # r = wildcards.run
    s = wildcards.sample
    if config['trim'] == False:
        return config['samples'][s]
    else:
        tmp = []
        if len(config['samples'][s]) == 2:
            tmp.append('%s/trim/%s/%s_val_1.fq.gz' % (RES_PATH,s,s))
            tmp.append('%s/trim/%s/%s_val_2.fq.gz' % (RES_PATH,s,s))
        else:
            tmp.append('%s/trim/%s/%s_trimmed.fq.gz' % (RES_PATH,s,s))
        return tmp

def align_salmon_targets(wildcards):
    ls = []
    for run in config["runs"]:
        if config["runs"][run]["type"] == "RS" and config["runs"][run]['samples']:
            for sample in config["runs"][run]['samples']:
                ls.append('%s/salmon/%s/quant.sf' % (RES_PATH, sample))
    return ls

rule align_Salmon:
    input:
        getAlignFastq
    output:
        "%s/salmon/{sample}/quant.sf" % RES_PATH
    params:
        index=config["salmon_index"],
        _inputs=lambda wildcards,input: "-1 %s -2 %s" % (input[0], input[1]) if len(input) == 2 else '-r %s' % input[0],
        output_path=lambda wildcards: "%s/salmon/%s/" % (RES_PATH, wildcards.sample),
        bootstrap=100,
        gcbias=lambda wildcards: "--gcBias" if len(config['samples'][wildcards.sample]) == 2 else "",
    log: "%s/logs/salmon/{sample}_align_Salmon.log" % RES_PATH
    message: "ALIGN: Align {wildcards.sample} to the genome by salmon"
    threads: salmon_threads
    shell:
        "salmon quant -i {params.index} -l A {params._inputs} -o {params.output_path} "
        "--numBootstraps {params.bootstrap} -p {threads} {params.gcbias} --validateMappings > {log} 2>&1"



