# trim adapter module
_threads=4

# ================================
# @auther: Xin Dong
# @email: xindong9511@gmail.com
# @date: Sep 2019
# ================================

def getTrimFastq(wildcards):
    if config['runs'][wildcards.run]['samples']:
        s = config['runs'][wildcards.run]['samples'][wildcards.sample]
    return s


def trim_targets(wildcards):
    ls = []
    for run in config["runs"]:
        if config["runs"][run]['type'] == "RS" and config["runs"][run]['samples']:
            for sample in config["runs"][run]['samples']:
                # print(run,sample)
                if len(config["runs"][run]['samples'][sample]) == 2:
                    ls.append('analysis/%s/samples/%s/trim/%s_val_1.fq.gz' % (run,sample,sample))
                    ls.append('analysis/%s/samples/%s/trim/%s_val_2.fq.gz' % (run,sample,sample))
                else:
                    ls.append('analysis/%s/samples/%s/trim/%s_trimmed.fq.gz' % (run,sample,sample))
                    ls.append('analysis/%s/samples/%s/trim/%s_trimmed_fastqc.zip' % (run,sample,sample))
                    ls.append('analysis/%s/samples/%s/trim/%s_trimmed_fastqc.html' % (run,sample,sample))
    return ls


rule trim_PairedEndAdapter:
    input:
        getTrimFastq
    output:
        'analysis/{run}/samples/{sample}/trim/{sample}_val_1.fq.gz',
        'analysis/{run}/samples/{sample}/trim/{sample}_val_2.fq.gz',
    params:
        quality=20,
        error_rate=0.1,
        stringency=8,
        length=20,
        output_dir=lambda wildcards: 'analysis/%s/samples/%s/trim/' % (wildcards.run, wildcards.sample),
        # file=lambda wildcards, input: '--paired %s %s' % (input[0],input[1]) if len(input) == 2 else '%s' % input, 
        basename=lambda wildcards: '%s' % wildcards.sample
    log: "analysis/{run}/log/{sample}_trim_PairedEndAdapter.log"
    message: "TRIM: Trim adaptor for paired-end sample - {wildcards.sample} " 
    threads: _threads
    shell:
        "trim_galore -q {params.quality} -j {threads} --phred33 --stringency {params.stringency} "
        "-e {params.error_rate} --gzip --length {params.length} --fastqc "
        "--basename {params.basename} -o {params.output_dir} --paired {input} > {log} 2>&1 "

rule trim_SingleEndAdapter:
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
    log: "analysis/{run}/log/{sample}_trim_SingleEndAdapter.log"
    message: "TRIM: Trim adaptor for paired-end sample - {wildcards.sample} "
    threads: _threads
    shell:
        "trim_galore -q {params.quality} -j {threads} --phred33 --stringency {params.stringency} "
        "-e {params.error_rate} --gzip --length {params.length} --fastqc "
        "{input} --basename {params.basename} -o {params.output_dir} > {log} 2>&1"


