# align salmon 

# ================================
# @auther: Xin Dong
# @email: xindong9511@gmail.com
# @date: Sep 2019
# ================================

salmon_threads = config['aligner_threads']
# global RES_PATH
# RES_PATH = config['res_path']

def align_salmon_targets(wildcards):
    ls = []
    for run in config["runs"]:
        if config["runs"][run]["type"] == "RS" and config["runs"][run]['samples']:
            for sample in config["runs"][run]['samples']:
                ls.append('%s/salmon/%s/%s.quant.sf' % (RES_PATH, sample, sample))
                ls.append('%s/salmon/%s/%s.quant.genes.sf' % (RES_PATH, sample, sample))
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

def getAlignrRNAFastq(wildcards):
    s = wildcards.sample
    if config['contam'] == True:
        for panel in config['contam_panel']:
            if len(config['samples'][s]) == 2:
                return ['%s/bowtie2_contam_panel/%s/%s.%s.unmapped.fq.1.gz' % (RES_PATH,s,s,panel),
                        '%s/bowtie2_contam_panel/%s/%s.%s.unmapped.fq.2.gz' % (RES_PATH,s,s,panel)]
            else:
                return ['%s/bowtie2_contam_panel/%s/%s.%s.unmapped.fq.gz' % (RES_PATH,s,s,panel)]
    else:
        tmp = getAlignRawFastq(wildcards)
        return tmp


def getAlignFastq(wildcards):
    s = wildcards.sample
    if config['aligner_rRNA'] == True:
        # first align to rRNA, then use the unmapped reads to align to genome
        tmp = []
        if len(config['samples'][s]) == 2:
            tmp.append('%s/STAR_rRNA/%s/%s.Unmapped.out.mate1' % (RES_PATH,s,s))
            tmp.append('%s/STAR_rRNA/%s/%s.Unmapped.out.mate2' % (RES_PATH,s,s))
        else:
            tmp.append('%s/STAR_rRNA/%s/%s.Unmapped.out.mate1' % (RES_PATH,s,s))
        return tmp
    else:
        if config['contam'] == True:
            for panel in config['contam_panel']:
                getAlignrRNAFastq(wildcards)
        else:
            tmp = getAlignRawFastq(wildcards)
            return tmp


rule align_rRNA_STAR:
    input:
        getAlignrRNAFastq
    output:
        "%s/STAR_rRNA/{sample}/{sample}.Unmapped.out.mate1" % RES_PATH,
        "%s/STAR_rRNA/{sample}/{sample}.Unmapped.out.mate2" % RES_PATH
    params:
        gz_support=lambda wildcards, input: "--readFilesCommand zcat" if str(input[0]).endswith('.gz') else "",
        prefix=lambda wildcards: "{res_path}/STAR_rRNA/{sample}/{sample}.".format(res_path=RES_PATH, sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample),
    threads: star_threads
    message: "ALIGN: Align rRNA reads of {wildcards.sample} to the rRNA genome by STAR"
    log:
        "%s/logs/STAR_rRNA/{sample}.align_rRNA_STAR.log" % RES_PATH
    shell:
        """
        STAR --runThreadN {threads} \
        --runMode alignReads \
        --genomeDir {config[STAR_rRNA_index]} \
        --readFilesIn {input} {params.gz_support} \
        --outFileNamePrefix {params.prefix} \
        --outReadsUnmapped Fastx \
        --outFilterMultimapNmax 10 \
        --winAnchorMultimapNmax 500 \
        --outFilterScoreMinOverLread 0 \
        --outFilterMatchNminOverLread 0 \
        --outSAMattributes All \
        --outSAMtype BAM Unsorted \
        --alignIntronMin 1 \
        --alignEndsType Local \
        --scoreGap 0 \
        --scoreGapNoncan 0 \
        --scoreGapGCAG 0 \
        --scoreGapATAC 0 \
        --scoreGenomicLengthLog2scale -1 \
        --chimFilter None \
        --chimOutType WithinBAM HardClip \
        --chimSegmentMin 5 \
        --chimJunctionOverhangMin 5 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreDropMax 80 \
        --chimNonchimScoreDropMin 20 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.05 \
        --limitOutSJcollapsed 10000000 \
        --limitBAMsortRAM 8000000000 > {log} 2>&1
        """


rule align_salmon:
    input:
        getAlignFastq
    output:
        "%s/salmon/{sample}/quant.sf" % RES_PATH,
        "%s/salmon/{sample}/quant.genes.sf" % RES_PATH
    params:
        index=config["salmon_index"],
        tx2genemap=lambda wildcards: "-g %s" % config['tx2genename'] if config['tx2genename'] else "",
        _inputs=lambda wildcards,input: "-1 %s -2 %s" % (input[0], input[1]) if len(input) == 2 else '-r %s' % input[0],
        output_path=lambda wildcards: "%s/salmon/%s/" % (RES_PATH, wildcards.sample),
        bootstrap=100,
        gcbias=lambda wildcards: "--gcBias" if len(config['samples'][wildcards.sample]) == 2 else "",
    log: "%s/logs/salmon/{sample}_align_Salmon.log" % RES_PATH
    message: "ALIGN: Align {wildcards.sample} to the genome by salmon"
    threads: salmon_threads
    shell:
        "salmon quant -i {params.index} {params.tx2genemap} "
        "-l A {params._inputs} -o {params.output_path} "
        "--numBootstraps {params.bootstrap} -p {threads} "
        "{params.gcbias} --validateMappings > {log} 2>&1"

rule align_salmonRenameTranscripts:
    input:
        "%s/salmon/{sample}/quant.sf" % RES_PATH
    output:
        "%s/salmon/{sample}/{sample}.quant.sf" % RES_PATH
    message: "ALIGN: rename salmon output of {wildcards.sample}"
    shell:
        "mv {input} {output}"

rule align_salmonRenameGenes:
    input:
        "%s/salmon/{sample}/quant.genes.sf" % RES_PATH
    output:
        "%s/salmon/{sample}/{sample}.quant.genes.sf" % RES_PATH
    message: "ALIGN: rename salmon output of {wildcards.sample}"
    shell:
        "mv {input} {output}"