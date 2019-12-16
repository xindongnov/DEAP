# align salmon 

# ================================
# @auther: Xin Dong
# @email: xindong9511@gmail.com
# @date: Sep 2019
# ================================

_threads = 8

def getAlignFastq(wildcards):
    r = wildcards.run
    s = wildcards.sample
    if config['trim'] == False:
        return config['runs'][r]['samples'][s]
    else:
        tmp = []
        if len(config['runs'][r]['samples'][s]) == 2:
            tmp.append('analysis/%s/samples/%s/trim/%s_val_1.fq.gz' % (r,s,s))
            tmp.append('analysis/%s/samples/%s/trim/%s_val_2.fq.gz' % (r,s,s))
        else:
            tmp.append('analysis/%s/samples/%s/trim/%s_trimmed.fq.gz' % (r,s,s))
        return tmp

def align_salmon_targets(wildcards):
    ls = []
    for run in config["runs"]:
        if config["runs"][run]["type"] == "RS" and config["runs"][run]['samples']:
            for sample in config["runs"][run]['samples']:
                ls.append('analysis/%s/samples/%s/align/quant.sf' % (run,sample))
    return ls

rule align_Salmon:
    input:
        getAlignFastq
    output:
        "analysis/{run}/samples/{sample}/align/quant.sf"
    params:
        index=config["salmon_index"],
        _inputs=lambda wildcards,input: "-1 %s -2 %s" % (input[0], input[1]) if len(input) == 2 else '-r %s' % input[0],
        output_path=lambda wildcards: "analysis/%s/samples/%s/align/" % (wildcards.run,wildcards.sample),
        bootstrap=100,
        gcbias=lambda wildcards: "--gcBias" if len(config['runs'][wildcards.run]['samples'][wildcards.sample]) == 2 else "",
    log: "analysis/{run}/log/{sample}_align_Salmon.log"
    message: "ALIGN: Align {wildcards.sample} to the genome by salmon"
    threads: _threads
    shell:
        "salmon quant -i {params.index} -l A {params._inputs} -o {params.output_path} "
        "--numBootstraps {params.bootstrap} -p {threads} {params.gcbias} --validateMappings > {log} 2>&1"



