# ==================================
# @File    :   align_STAR.snakefile
# @Time    :   2023/07/14 21:28:11
# @Author  :   Xin Dong @ref:RIMA
# @Contact :   xindong9511@gmail.com
# @License :   (C)Copyright 2020-2023, XinDong
# ==================================



# alignment by STAR

_threads = 8


def align_STAR_targets(wildcards):
    ls = []
    for run in config["runs"]:
        if config["runs"][run]["type"] == "RS" and config["runs"][run]['samples']:
            for sample in config["runs"][run]['samples']:
                ls.append('analysis/STAR/{sample}/{sample}.unsorted.bam'.format(sample=sample))
                ls.append('analysis/STAR/{sample}/{sample}.sorted.bam'.format(sample=sample))
                ls.append('analysis/STAR/{sample}/{sample}.sorted.bam.bai'.format(sample=sample))
    return ls

def getAlignFastq(wildcards):
    s = wildcards.sample
    if config['trim'] == False:
        return config['samples'][s]
    else:
        tmp = []
        if len(config['samples'][s]) == 2:
            tmp.append('%s/trim/%s/%s_val_1.fq.gz' % (res_path,s,s))
            tmp.append('%s/trim/%s/%s_val_2.fq.gz' % (res_path,s,s))
        else:
            tmp.append('%s/trim/%s/%s_trimmed.fq.gz' % (res_path,s,s))
        return tmp

rule align_STAR:
    input:
        getAlignFastq
    output:
        unsortedBAM = "analysis/STAR/{sample}/{sample}.unsorted.bam",
        # sortedBAM = "analysis/STAR/{sample}/{sample}.sorted.bam",
        transcriptomeBAM = "analysis/STAR/{sample}/{sample}.transcriptome.bam",
        junction_file = "analysis/STAR/{sample}/{sample}.Chimeric.out.junction",
        #counts = "analysis/STAR/{sample}/{sample}.counts.tab",
        log_file = "analysis/STAR/{sample}/{sample}.Log.final.out"
    params:
        gz_support=lambda wildcards: "--readFilesCommand zcat" if config["samples"][wildcards.sample][0][-3:] == '.gz' else "",
        prefix=lambda wildcards: "analysis/STAR/{sample}".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample),
        # keepPairs = _keepPairs
    threads: _threads
    message: "ALIGN: Align {wildcards.sample} to the genome by STAR"
    log:
        "analysis/logs/STAR/{sample}.star_align.log"
    shell:
        "STAR --runThreadN {threads} "
        "--genomeDir {config[STAR_index]} "
        "--outReadsUnmapped Fastx "
        "--chimSegmentMin 12 "
        "--chimJunctionOverhangMin 12 "
        # "--chimOutJunctionFormat 1 "
        "--alignSJDBoverhangMin 10 "
        "--alignMatesGapMax 1000000 "
        "--alignIntronMax 1000000 "
        "--alignSJstitchMismatchNmax 5 -1 5 5 "
        "--outSAMstrandField intronMotif "
        "--outSAMunmapped Within "
        "--outSAMtype BAM Unsorted "
        "--readFilesIn {input} "
        # "--chimMultimapScoreRange 10 "
        # "--chimMultimapNmax 10 "
        # "--chimNonchimScoreDropMin 10 "
        # "--peOverlapNbasesMin 12 "
        # "--peOverlapMMp 0.1 "
        "--genomeLoad NoSharedMemory "
        "--outSAMheaderHD @HD VN:1.4 "
        "--twopassMode Basic "
        "{params.gz_support} "
        "--outFileNamePrefix {params.prefix} "
        "--quantMode TranscriptomeSAM"
        " && mv {params.prefix}Aligned.out.bam {output.unsortedBAM}"
        # " && samtools sort -T {params.prefix}TMP -o {output.sortedBAM} -@ {threads} {output.unsortedBAM}"
        " && mv {params.prefix}Aligned.toTranscriptome.out.bam {output.transcriptomeBAM}"
        #" && mv {params.prefix}ReadsPerGene.out.tab {output.counts}"
        " && mv {params.prefix}Chimeric.out.junction {output.junction_file}"
        " && mv {params.prefix}Log.final.out {output.log_file}"

rule align_STAR_sort_bam:
    input:
        "analysis/STAR/{sample}/{sample}.unsorted.bam",
    output:
        sortedBAM = "analysis/STAR/{sample}/{sample}.sorted.bam",
    message: "ALIGN: Sorting {wildcards.sample}.unsorted.bam"
    params:
        prefix=lambda wildcards: "analysis/STAR/{sample}".format(sample=wildcards.sample),
    shell:
        "samtools sort -T {params.prefix}TMP -o {output.sortedBAM} -@ {threads} {input}"


rule align_STAR_index_bam:
    # """INDEX the {sample}.sorted.bam file"""
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/STAR/{sample}/{sample}.sorted.bam.bai"
    message: "ALIGN: Indexing {wildcards.sample}.sorted.bam"
    shell:
        "samtools index {input}"

# rule generate_STAR_report:
#     input:
#         star_log_files=expand( "analysis/STAR/{sample}/{sample}.Log.final.out", sample=config["ordered_sample_list"] ),
#         star_gene_count_files=expand( "analysis/STAR/{sample}/{sample}.counts.tab", sample=config["ordered_sample_list"] ),
#         force_run_upon_meta_change = config['metasheet'],
#         force_run_upon_config_change = config['config_file']
#     output:
#         csv="analysis/" + config["token"] + "/STAR/STAR_Align_Report.csv",
#         png="analysis/" + config["token"] + "/STAR/STAR_Align_Report.png",
#         gene_counts="analysis/" + config["token"] + "/STAR/STAR_Gene_Counts.csv"
#     message: "ALIGN: Generating STAR report"
#     run:
#         log_files = " -l ".join( input.star_log_files )
#         count_files = " -f ".join( input.star_gene_count_files )
#         shell( "perl DEAP/modules/scripts/STAR_reports.pl -l {log_files} 1 > {output.csv}" )
#         shell( "Rscript DEAP/modules/scripts/map_stats.R {output.csv} {output.png}" )
#         shell( "perl DEAP/modules/scripts/raw_and_fpkm_count_matrix.pl -f {count_files} 1 > {output.gene_counts}" )

# rule batch_effect_removal_star:
#     input:
#         starmat = "analysis/" + config["token"] + "/STAR/STAR_Gene_Counts.csv",
#         annotFile = config["metasheet"]
#     output:
#         starcsvoutput="analysis/" + config["token"] + "/STAR/batch_corrected_STAR_Gene_Counts.csv",
#         starpdfoutput="analysis/" + config["token"] + "/STAR/star_combat_qc.pdf"
#     params:
#         batch_column="batch",
#         datatype = "star"
#     message: "Removing batch effect from STAR Gene Count matrix, if errors, check metasheet for batches, refer to README for specifics"
#     benchmark:
#         "benchmarks/" + config["token"] + "/batch_effect_removal_star.txt"
#     shell:
#         "Rscript viper/modules/scripts/batch_effect_removal.R {input.starmat} {input.annotFile} "
#         "{params.batch_column} {params.datatype} {output.starcsvoutput} {output.starpdfoutput} "
#         " && mv {input.starmat} analysis/{config[token]}/STAR/without_batch_correction_STAR_Gene_Counts.csv"


# rule run_STAR_fusion:
#     input:
#         bam="analysis/STAR/{sample}/{sample}.sorted.bam" #just to make sure STAR output is available before STAR_Fusion
#     output:
#         protected("analysis/STAR_Fusion/{sample}/{sample}.fusion_candidates.final"),
#         protected("analysis/STAR_Fusion/{sample}/{sample}.fusion_candidates.final.abridged")
#     log:
#         "analysis/STAR_Fusion/{sample}/{sample}.star_fusion.log"
#     message: "Running STAR fusion on {wildcards.sample}"
#     benchmark:
#         "benchmarks/{sample}/{sample}.run_STAR_fusion.txt"
#     shell:
#         "STAR-Fusion --chimeric_junction analysis/STAR/{wildcards.sample}/{wildcards.sample}.Chimeric.out.junction "
#         "--genome_lib_dir {config[genome_lib_dir]} --output_dir analysis/STAR_Fusion/{wildcards.sample} >& {log}"
#         " && mv analysis/STAR_Fusion/{wildcards.sample}/star-fusion.fusion_candidates.final {output[0]}"
#         " && mv analysis/STAR_Fusion/{wildcards.sample}/star-fusion.fusion_candidates.final.abridged {output[1]}"
#         " && touch {output[1]}" # For some sample, final.abridged is created but not .final file; temp hack before further investigate into this


# rule run_STAR_fusion_report:
#     input:
#         sf_list = expand("analysis/STAR_Fusion/{sample}/{sample}.fusion_candidates.final.abridged", sample=config["ordered_sample_list"]),
#         force_run_upon_meta_change = config['metasheet'],
#         force_run_upon_config_change = config['config_file']
#     output:
#         csv="analysis/" + config["token"] + "/STAR_Fusion/STAR_Fusion_Report.csv",
#         png="analysis/" + config["token"] + "/STAR_Fusion/STAR_Fusion_Report.png"
#     message: "Generating STAR fusion report"
#     benchmark:
#         "benchmarks/" + config["token"] + "/run_STAR_fusion_report.txt"
#     shell:
#         "python viper/modules/scripts/STAR_Fusion_report.py -f {input.sf_list} 1>{output.csv} "
#         "&& Rscript viper/modules/scripts/STAR_Fusion_report.R {output.csv} {output.png}"


# rule run_rRNA_STAR:
#     input:
#         getFastq
#     output:
#         bam=protected("analysis/STAR_rRNA/{sample}/{sample}.sorted.bam"),
#         log_file="analysis/STAR_rRNA/{sample}/{sample}.Log.final.out"
#     params:
#         stranded=rRNA_strand_command,
#         gz_support=gz_command,
#         prefix=lambda wildcards: "analysis/STAR_rRNA/{sample}/{sample}".format(sample=wildcards.sample),
#         readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample)
#     threads: 8
#     message: "Running rRNA STAR for {wildcards.sample}"
#     benchmark:
#         "benchmarks/{sample}/{sample}.run_rRNA_STAR.txt"
#     shell:
#         "STAR --runMode alignReads --runThreadN {threads}"
#         " --genomeDir {config[star_rRNA_index]}"
#         " --readFilesIn {input} {params.gz_support}"
#         " --outFileNamePrefix {params.prefix}."
#         " --outSAMmode Full --outSAMattributes All {params.stranded}"
#         " --outSAMattrRGline {params.readgroup}"
#         " --outSAMtype BAM SortedByCoordinate"
#         " --limitBAMsortRAM 45000000000"
#         " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"

# rule index_STAR_rRNA_bam:
#     """INDEX the STAR_rRNA/{sample}.sorted.bam file"""
#     input:
#         "analysis/STAR_rRNA/{sample}/{sample}.sorted.bam"
#     output:
#         "analysis/STAR_rRNA/{sample}/{sample}.sorted.bam.bai"
#     message: "Indexing STAR_rRNA {wildcards.sample}.sorted.bam"
#     benchmark:
#         "benchmarks/{sample}/{sample}.index_STAR_rRNA_bam.txt"
#     shell:
#         "samtools index {input}"

# rule generate_rRNA_STAR_report:
#     input:
#         star_log_files=expand( "analysis/STAR_rRNA/{sample}/{sample}.Log.final.out", sample=config["ordered_sample_list"] ),
#         force_run_upon_meta_change = config['metasheet'],
#         force_run_upon_config_change = config['config_file']
#     output:
#         csv="analysis/" + config["token"] + "/STAR_rRNA/STAR_rRNA_Align_Report.csv",
#         png="analysis/" + config["token"] + "/STAR_rRNA/STAR_rRNA_Align_Report.png"
#     message: "Generating STAR rRNA report"
#     benchmark:
#         "benchmarks/" + config["token"] + "/run_rRNA_STAR_report.txt"
#     run:
#         log_files = " -l ".join( input.star_log_files )
#         shell( "perl viper/modules/scripts/STAR_reports.pl -l {log_files} 1>{output.csv}" )
#         shell( "Rscript viper/modules/scripts/map_stats_rRNA.R {output.csv} {output.png}" )

# rule align_SJtab2JunctionsBed:
#     """Convert STAR's SJ.out.tab to (tophat) junctions.bed BED12 format"""
#     input:
#         "analysis/STAR/{sample}/{sample}.SJ.out.tab"
#     output:
#         "analysis/STAR/{sample}/{sample}.junctions.bed"
#     benchmark:
#         "benchmarks/{sample}/{sample}.align_SJtab2JunctionsBed.txt"
#     shell:
#         "viper/modules/scripts/STAR_SJtab2JunctionsBed.py -f {input} > {output}"
