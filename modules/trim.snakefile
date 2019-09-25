# trim adapter module
_threads=4

def getTrimFastq(wildcards):
    if config['RS_runs'][wildcards.run]['samples']:
        s = config['RS_runs'][wildcards.run]['samples'][wildcards.sample]
    return s


def trim_targets(wildcards):
    ls = []
    for run in config["RS_runs"]:
        if config["RS_runs"][run]['samples']:
            for sample in config["RS_runs"][run]['samples']:
                # print(run,sample)
                if len(config["RS_runs"][run]['samples'][sample]) == 2:
                    ls.append('analysis/%s/samples/%s/trim/%s_R1_val_1.fq.gz' % (run,sample,sample))
                    ls.append('analysis/%s/samples/%s/trim/%s_R2_val_2.fq.gz' % (run,sample,sample))
                else:
                    ls.append('analysis/%s/samples/%s/trim/%s_trimmed.fq.gz' % (run,sample,sample))
                    ls.append('analysis/%s/samples/%s/trim/%s_trimmed_fastqc.zip' % (run,sample,sample))
                    ls.append('analysis/%s/samples/%s/trim/%s_trimmed_fastqc.html' % (run,sample,sample))
    return ls


rule trim_paired_adapter:
    input:
        getTrimFastq
    output:
        'analysis/{run}/samples/{sample}/trim/{sample}_R1_val_1.fq.gz',
        'analysis/{run}/samples/{sample}/trim/{sample}_R2_val_2.fq.gz',
    params:
        quality=20,
        error_rate=0.1,
        stringency=8,
        length=20,
        output_dir=lambda wildcards: 'analysis/%s/samples/%s/trim/' % (wildcards.run, wildcards.sample),
        # file=lambda wildcards, input: '--paired %s %s' % (input[0],input[1]) if len(input) == 2 else '%s' % input, 
        basename=lambda wildcards: '%s' % wildcards.sample
    # log: "analysis/{run}/log/trim/{sample}.log"
    message: "TRIM: Trim adaptor for {wildcards.sample} " 
    threads: _threads
    shell:
        "trim_galore -q {params.quality} -j {threads} --phred33 --stringency {params.stringency} "
        "-e {params.error_rate} --gzip --length {params.length} --fastqc "
        "--basename {params.basename} -o {params.output_dir} --paired {input} "

rule trim_single_adapter:
    input:
        getTrimFastq
    output:
        'analysis/{run}/samples/{sample}/trim/{sample}_trimmed.fq.gz',
        'analysis/{run}/samples/{sample}/trim/{sample}_trimmed_fastqc.zip',
        'analysis/{run}/samples/{sample}/trim/{sample}_trimmed_fastqc.html',
    params:
        quality=20,
        error_rate=0.1,
        stringency=8,
        length=20,
        output_dir=lambda wildcards: 'analysis/%s/samples/%s/trim/' % (wildcards.run, wildcards.sample),
        basename=lambda wildcards: '%s' % wildcards.sample
    # log: "analysis/{run}/log/trim/{sample}.log"
    message: "TRIM: Trim adaptor for {wildcards.sample} "
    threads: _threads
    shell:
        "trim_galore -q {params.quality} -j {threads} --phred33 --stringency {params.stringency} "
        "-e {params.error_rate} --gzip --length {params.length} --fastqc "
        "{input} --basename {params.basename} -o {params.output_dir} "


