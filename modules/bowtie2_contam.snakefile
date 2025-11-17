# ==================================
# @File    :   align_STAR.snakefile
# @Time    :   2023/07/14 21:28:11
# @Author  :   Xin Dong @ref:RIMA
# @Contact :   xindong9511@gmail.com
# @License :   (C)Copyright 2020-2023, XinDong
# ==================================



# alignment by bowtie2 with contamination panel

bowtie2_threads = config['aligner_threads']


def bowtie2_contam_targets(wildcards):
    ls = []
    for run in config["runs"]:
        if config["runs"][run]["type"] == "RS" and config["runs"][run]['samples']:
            for sample in config["runs"][run]['samples']:
                if config['contam'] == True:
                    for panel in config['contam_panel']:
                        if len(config['samples'][sample]) == 2:
                            ls.append('{res_path}/bowtie2_contam_panel/{sample}/{sample}.{panel}_PE.bam'.format(res_path=RES_PATH, panel=panel, sample=sample))
                            ls.append('{res_path}/bowtie2_contam_panel/{sample}/{sample}.{panel}.unmapped.fq.1.gz'.format(res_path=RES_PATH, panel=panel, sample=sample))
                            ls.append('{res_path}/bowtie2_contam_panel/{sample}/{sample}.{panel}.unmapped.fq.2.gz'.format(res_path=RES_PATH, panel=panel, sample=sample))
                            ls.append('{res_path}/bowtie2_contam_panel/{sample}/{sample}.{panel}_PE_flagstat.txt'.format(res_path=RES_PATH, panel=panel, sample=sample))
                        else:
                            ls.append('{res_path}/bowtie2_contam_panel/{sample}/{sample}.{panel}_SE.bam'.format(res_path=RES_PATH, panel=panel, sample=sample))
                            ls.append('{res_path}/bowtie2_contam_panel/{sample}/{sample}.{panel}.unmapped.fq.gz'.format(res_path=RES_PATH, panel=panel, sample=sample))
                            ls.append('{res_path}/bowtie2_contam_panel/{sample}/{sample}.{panel}_SE_flagstat.txt'.format(res_path=RES_PATH, panel=panel, sample=sample))
    return ls

def getAlignRawFastq(wildcards):
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


rule contamination_mapping_bowtie2_PE:
    input:
        getAlignRawFastq
    output:
        temp("%s/bowtie2_contam_panel/{sample}/{sample}.{panel}_PE.sam" % RES_PATH),
        "%s/bowtie2_contam_panel/{sample}/{sample}.{panel}.unmapped.fq.1.gz" % RES_PATH,
        "%s/bowtie2_contam_panel/{sample}/{sample}.{panel}.unmapped.fq.2.gz" % RES_PATH,
    params:
        index=lambda wildcards: config['contamination_panel'][wildcards.panel],
        res_path=RES_PATH
    threads: bowtie2_threads
    message: "ALIGN: Align {wildcards.sample} to contamination panel by bowtie2"
    log:
        "%s/logs/bowtie2_contam_panel/{sample}_contamination_mapping_bowtie2_PE_{panel}.log" % RES_PATH
    shell:
        """
        bowtie2 -x {params.index} -p {threads}  -1 {input[0]} -2 {input[1]} \
        --un-conc-gz {params.res_path}/bowtie2_contam_panel/{wildcards.sample}/{wildcards.sample}.{wildcards.panel}.unmapped.fq.gz \
        -S {params.res_path}/bowtie2_contam_panel/{wildcards.sample}/{wildcards.sample}.{wildcards.panel}_PE.sam > {log} 2>&1

        """

rule contamination_mapping_bowtie2_SE:
    input:
        getAlignRawFastq
    output:
        temp("%s/bowtie2_contam_panel/{sample}/{sample}.{panel}_SE.sam" % RES_PATH),
        "%s/bowtie2_contam_panel/{sample}/{sample}.{panel}.unmapped.fq.gz" % RES_PATH,
    params:
        index=lambda wildcards: config['contamination_panel'][wildcards.panel],
        res_path=RES_PATH
    threads: bowtie2_threads
    message: "ALIGN: Align {wildcards.sample} to contamination panel by bowtie2"
    log:
        "%s/logs/bowtie2_contam_panel/{sample}_contamination_mapping_bowtie2_SE_{panel}.log" % RES_PATH
    shell:
        """
        bowtie2 -x {params.index} -p {threads}  -U {input[0]} \
        --un-gz {params.res_path}/bowtie2_contam_panel/{wildcards.sample}/{wildcards.sample}.{wildcards.panel}.unmapped.fq.gz \
        -S {params.res_path}/bowtie2_contam_panel/{wildcards.sample}/{wildcards.sample}.{wildcards.panel}_SE.sam > {log} 2>&1
        """

rule contamination_sam2bam:
    input:
        "%s/bowtie2_contam_panel/{sample}/{sample}.{panel}_{type}.sam" % RES_PATH,
    output:
        "%s/bowtie2_contam_panel/{sample}/{sample}.{panel}_{type}.bam" % RES_PATH,
    threads: 16
    message: "SAM2BAM: Convert {wildcards.sample} bowtie2 contamination SAM to BAM"
    log:
        "%s/logs/bowtie2_contam_panel/{sample}_sam2bam_bowtie2_contam_{panel}_{type}.log" % RES_PATH
    shell:
        """
        samtools view -@ {threads} -bS {input} > {output} 2> {log}
        """

rule contamination_bam_flagstat:
    input:
        "%s/bowtie2_contam_panel/{sample}/{sample}.{panel}_{type}.bam" % RES_PATH,
    output:
        "%s/bowtie2_contam_panel/{sample}/{sample}.{panel}_{type}_flagstat.txt" % RES_PATH,
    threads: 4
    message: "FLAGSTAT: Get {wildcards.sample} bowtie2 contamination BAM flagstat"
    log:
        "%s/logs/bowtie2_contam_panel/{sample}_flagstat_bowtie2_contam_{panel}_{type}.log" % RES_PATH
    shell:
        """
        samtools flagstat -@ {threads} {input} > {output} 2> {log}
        """