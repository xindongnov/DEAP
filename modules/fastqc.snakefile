#MODULE: fastqc- sequence quality scores
#RULES:
#    sample_fastq: subsample 100k reads from each sample to perform fastqc analysis
_logfile="analysis/logs/fastqc.log"

def fastqc_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/fastqc/%s/%s_perSeqGC.txt" % (sample,sample))
        ls.append("analysis/fastqc/%s/%s_perSeqQual.txt" % (sample,sample))
        ls.append("analysis/fastqc/%s/%s_stats.csv" % (sample,sample))
        ls.append("analysis/fastqc/%s/%s_perSeqGC.png" % (sample,sample))
        ls.append("analysis/fastqc/%s/%s_perSeqGC_thumb.png" % (sample,sample))

    ls.append("analysis/fastqc/fastqc.csv")
    return ls

def getFastqcInput(wildcards):
    """Get the input file for fastqc. It's either the 100k.fastq sample or the 
     original bam"""
    s = wildcards.sample
    first_file = config["samples"][wildcards.sample][0]
    if first_file.endswith('.bam'):
        #CLEANER to check for .bam vs (.fastq, fq, fq.gz, etc)
        #ret = [first_file]
        #HACK to get fastqc to name things correctly.  IF we returned
        #the bam file, fastqc will name the output files accordingly
        ret = ["analysis/align/%s/%s_100k.bam" % (s,s)]
    else:
        #HACK: need to give the EVALUATED cannonical path
        ret = ["analysis/align/%s/%s_100k.fastq" % (s,s)]
    #print(ret)
    return ret

def getFastq3(wildcards):
    """Get associated fastqs for each sample.  
    NOTE: if PE, take the first pair"""
    s = config["samples"][wildcards.sample]
    return s[0]

rule fastqc_all:
    input:
        fastqc_targets

rule sample_fastq:
    """Subsample fastq"""
    input:
        getFastq3
    output:
        temp("analysis/align/{sample}/{sample}_100k.fastq")
    params:
        seed=11,
        #how many to sample
        size=100000
    message: "FASTQC: sample_fastq"
    log:_logfile
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "seqtk sample -s {params.seed} {input} {params.size} > {output} 2>>{log}"

rule sample_bam:
    """USED only when the sample-input is a bam file.  
    Subsample bam to 100k reads"""
    input:
        getFastqcBam
    params:
        n=100000
    output:
        temp("analysis/align/{sample}/{sample}_100k.bam")
    message: "FASTQC: sampling 100k reads from bam"
    log:_logfile
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "cidc_chips/modules/scripts/sampleBam.sh -i {input} -n {params.n} -o {output}"

rule convertBamToFastq:
    """USED only when the sample-input is a bam file."""
    input:
        "analysis/align/{sample}/{sample}_100k.bam"
    output:
        temp("analysis/align/{sample}/{sample}_100k.bam.fastq")
    message: "FASTQC: converting 100k.bam to 100k.fastq"
    log:_logfile
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "bamToFastq -i {input} -fq {output} 2>> {log}"
    
rule call_fastqc:
    """CALL FASTQC on each sub-sample"""
    input:
        getFastqcInput
    output:
        #MAKE temp
        "analysis/fastqc/{sample}/{sample}_100k_fastqc/fastqc_data.txt",
        temp("analysis/fastqc/{sample}/{sample}_100k_fastqc.html"),
        temp("analysis/fastqc/{sample}/{sample}_100k_fastqc.zip")
    #threads:
    params:
        sample = lambda wildcards: wildcards.sample
    message: "FASTQC: call fastqc"
    log:_logfile
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "fastqc {input} --extract -o analysis/fastqc/{params.sample} 2>>{log}"

rule get_PerSequenceQual:
    """extract per sequence quality from fastqc_data.txt"""
    input:
        "analysis/fastqc/{sample}/{sample}_100k_fastqc/fastqc_data.txt"
    output:
        "analysis/fastqc/{sample}/{sample}_perSeqQual.txt"
    params:
        #DON'T forget quotes
        section="'Per sequence quality'"
    message: "FASTQC: get_PerSequenceQual"
    log:_logfile
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "cidc_chips/modules/scripts/fastqc_data_extract.py -f {input} -s {params.section} > {output} 2>>{log}"

rule get_PerSequenceGC:
    """extract per sequence GC contentfrom fastqc_data.txt"""
    input:
        "analysis/fastqc/{sample}/{sample}_100k_fastqc/fastqc_data.txt"
    output:
        "analysis/fastqc/{sample}/{sample}_perSeqGC.txt"
    params:
        #DON'T forget quotes
        section="'Per sequence GC content'"
    message: "FASTQC: get_PerSequenceGC"
    log:_logfile
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "cidc_chips/modules/scripts/fastqc_data_extract.py -f {input} -s {params.section} > {output} 2>>{log}"

rule extract_FastQCStats:
    """extract per sequence GC content, and seq qual stats from fastqc run"""
    input:
        gc = "analysis/fastqc/{sample}/{sample}_perSeqGC.txt",
        qual = "analysis/fastqc/{sample}/{sample}_perSeqQual.txt"
    output:
        "analysis/fastqc/{sample}/{sample}_stats.csv"
    message: "FASTQC: extract_FastQCStats"
    log:_logfile
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "cidc_chips/modules/scripts/fastqc_stats.py -a {input.qual} -b {input.gc} > {output} 2>>{log}"

rule collect_fastQCStats:
    """Collect and parse out the fastqc stats for the ALL of the samples"""
    input:
        expand("analysis/fastqc/{sample}/{sample}_stats.csv", sample=sorted(config["samples"]))
    output:
        "analysis/fastqc/fastqc.csv"
    message: "FASTQC: collect and parse ALL mapping stats"
    log: _logfile
    #conda: "../envs/fastqc/fastqc.yaml"
    run:
        files = " -f ".join(input)
        shell("cidc_chips/modules/scripts/fastqc_getFastQCStats.py -f {files} > {output} 2>>{log}")

rule plot_fastQC_GC:
    """Plots the GC distribution of the sample according to data in
    perSeqGC.txt.  
    Generates a full-size image and *thumbnail image* (embedded into report)
    """
    input:
        gc = "analysis/fastqc/{sample}/{sample}_perSeqGC.txt",
    output:
        png="analysis/fastqc/{sample}/{sample}_perSeqGC.png",
        thumb="analysis/fastqc/{sample}/{sample}_perSeqGC_thumb.png",
    message:
        "FASTQC: generating GC content distrib. plots"
    log: _logfile
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "Rscript cidc_chips/modules/scripts/fastqc_plotGC.R {input.gc} {output.png} {output.thumb}"
