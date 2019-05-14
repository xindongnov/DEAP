# trim adapter module
_threads=4

def getTrimFastq(wildcards):
    return config['RS_runs'][wildcards.run]['samples'][wildcards.sample]


def trim_targets(wildcards):
    ls = []
    for run in config["RS_runs"]:
        for sample in config["RS_runs"][run]['samples']:
            # print(run,sample)
            if len(config["RS_runs"][run]['samples'][sample]) == 2:
                ls.append('analysis/%s/%s/trim/%s_val_1_trimmed.fq.gz' % (run,sample,sample))
                ls.append('analysis/%s/%s/trim/%s_val_2_trimmed.fq.gz' % (run,sample,sample))
            else:
                ls.append('analysis/%s/%s/trim/%s_trimmed.fq.gz' % (run,sample,sample))
                ls.append('analysis/%s/%s/trim/%s_trimmed_fastqc.zip' % (run,sample,sample))
                ls.append('analysis/%s/%s/trim/%s_trimmed_fastqc.html' % (run,sample,sample))
    return ls


rule trim_paired_adapter:
    input:
        getTrimFastq
    output:
        'analysis/{run}/{sample}/trim/{sample}_val_1_trimmed.fq.gz',
        'analysis/{run}/{sample}/trim/{sample}_val_2_trimmed.fq.gz',
    params:
        quality=20,
        error_rate=0.1,
        stringency=8,
        length=20,
        output_dir=lambda wildcards: 'analysis/%s/%s/trim/' % (wildcards.run, wildcards.sample),
        # file=lambda wildcards, input: '--paired %s %s' % (input[0],input[1]) if len(input) == 2 else '%s' % input, 
        thread=_threads,
        basename=lambda wildcards: '%s' % wildcards.sample
    log: "analysis/log/{run}/trim/{sample}.log"
    message: "TRIM: Trim adaptor for {wildcards.sample} " 
    shell:
        "trim_galore -q {params.quality} -j {params.thread} --phred33 --stringency {params.stringency} "
        "-e {params.error_rate} --gzip --length {params.length} --fastqc "
        "--basename {params.basename} -o {params.output_dir} --paired {input} "

rule trim_single_adapter:
    input:
        getTrimFastq
    output:
        'analysis/{run}/{sample}/trim/{sample}_trimmed.fq.gz',
        'analysis/{run}/{sample}/trim/{sample}_trimmed_fastqc.zip',
        'analysis/{run}/{sample}/trim/{sample}_trimmed_fastqc.html',
    params:
        quality=20,
        error_rate=0.1,
        stringency=8,
        length=20,
        output_dir=lambda wildcards: 'analysis/%s/%s/trim/' % (wildcards.run, wildcards.sample),
        thread=_threads,
        basename=lambda wildcards: '%s' % wildcards.sample
    log: "analysis/log/{run}/trim/{sample}.log"
    message: "TRIM: Trim adaptor for {wildcards.sample} "
    shell:
        "trim_galore -q {params.quality} -j {params.thread} --phred33 --stringency {params.stringency} "
        "-e {params.error_rate} --gzip --length {params.length} --fastqc "
        "{input} --basename {params.basename} -o {params.output_dir} "


