# align salmon 
_threads = 8

def getAlignFastq(wildcards):
    r = wildcards.run
    s = wildcards.sample
    if config['trim'] == False:
        return config['RS_runs'][r]['samples'][s]
    else:
        tmp = []
        if len(config['RS_runs'][r]['samples'][s]) == 2:
            tmp.append('analysis/%s/samples/%s/trim/%s_R1_val_1.fq.gz' % (r,s,s))
            tmp.append('analysis/%s/samples/%s/trim/%s_R2_val_2.fq.gz' % (r,s,s))
        else:
            tmp.append('analysis/%s/samples/%s/trim/%s_trimmed.fq.gz' % (r,s,s))
        return tmp

def align_salmon_targets(wildcards):
    ls = []
    for run in config["RS_runs"]:
        if config["RS_runs"][run]['samples']:
            for sample in config["RS_runs"][run]['samples']:
                ls.append('analysis/%s/samples/%s/align/quant.sf' % (run,sample))
    return ls

rule align_salmon:
    input:
        getAlignFastq
    output:
        "analysis/{run}/samples/{sample}/align/quant.sf"
    params:
        index=config["salmon_index"],
        _inputs=lambda wildcards,input: "-1 %s -2 %s" % (input[0], input[1]) if len(input) == 2 else '-r %s' % input[0],
        output_path=lambda wildcards: "analysis/%s/samples/%s/align/" % (wildcards.run,wildcards.sample),
        bootstrap=100,
        gcbias=lambda wildcards: "--gcBias" if len(config['RS_runs'][wildcards.run]['samples'][wildcards.sample]) == 2 else "",
    # log: "analysis/{run}/log/align/{sample}.log"
    message: "ALIGN: Align {wildcards.sample} to the genome "
    threads: _threads
    shell:
        "salmon quant -i {params.index} -l A {params._inputs} -o {params.output_path} "
        "--numBootstraps {params.bootstrap} -p {threads} {params.gcbias} --validateMappings"



