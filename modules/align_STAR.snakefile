# ==================================
# @File    :   align_STAR.snakefile
# @Time    :   2023/07/14 21:28:11
# @Author  :   Xin Dong @ref:RIMA
# @Contact :   xindong9511@gmail.com
# @License :   (C)Copyright 2020-2023, XinDong
# ==================================



# alignment by STAR

star_threads = config['aligner_threads']


def align_STAR_targets(wildcards):
    ls = []
    for run in config["runs"]:
        if config["runs"][run]["type"] == "RS" and config["runs"][run]['samples']:
            for sample in config["runs"][run]['samples']:
                # if config['contam'] == True:
                #     for panel in config['contam_panel']:
                #         ls.append('{res_path}/bowtie2_contam_panel/{sample}/{panel}.sam'.format(res_path=RES_PATH, panel=panel, sample=sample))
                #         if len(config['samples'][sample]) == 2:
                #             ls.append('{res_path}/bowtie2_contam_panel/{sample}/{panel}.unmapped.1.fq'.format(res_path=RES_PATH, panel=panel, sample=sample))
                #             ls.append('{res_path}/bowtie2_contam_panel/{sample}/{panel}.unmapped.2.fq'.format(res_path=RES_PATH, panel=panel, sample=sample))
                #         else:
                #             ls.append('{res_path}/bowtie2_contam_panel/{sample}/{panel}.unmapped.fq'.format(res_path=RES_PATH, panel=panel, sample=sample))
                if config['aligner_rRNA'] == True:
                    ls.append('{res_path}/STAR_rRNA/{sample}/{sample}.Unmapped.out.mate1'.format(res_path=RES_PATH, sample=sample))
                    if len(config['samples'][sample]) == 2:
                        ls.append('{res_path}/STAR_rRNA/{sample}/{sample}.Unmapped.out.mate2'.format(res_path=RES_PATH, sample=sample))
                ls.append('{res_path}/STAR/{sample}/{sample}.Aligned.out.bam'.format(res_path=RES_PATH, sample=sample))
                ls.append('{res_path}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam'.format(res_path=RES_PATH, sample=sample))
                ls.append('{res_path}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bw'.format(res_path=RES_PATH, sample=sample))
                ls.append('{res_path}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai'.format(res_path=RES_PATH, sample=sample))
                ls.append('{res_path}/STAR/{sample}/{sample}.Aligned.toTranscriptome.out.bam'.format(res_path=RES_PATH, sample=sample))
                # ls.append('{res_path}/STAR/{sample}/{sample}.Aligned.toTranscriptome.out.bam.bai'.format(res_path=RES_PATH, sample=sample))
                ls.append('{res_path}/STAR/{sample}/{sample}.ReadsPerTranscript.byfeatureCounts.txt'.format(res_path=RES_PATH, sample=sample))
                ls.append('{res_path}/STAR/{sample}/{sample}.ReadsPerGeneID.byfeatureCounts.txt'.format(res_path=RES_PATH, sample=sample))
                # ls.append('{res_path}/STAR/{sample}/{sample}.ReadsPerGeneName.byfeatureCounts.txt'.format(res_path=RES_PATH, sample=sample))
                ls.append('{res_path}/STAR/{sample}/{sample}.rsem.genes.results'.format(res_path=RES_PATH, sample=sample))
                ls.append('{res_path}/STAR/{sample}/{sample}.rsem.isoforms.results'.format(res_path=RES_PATH, sample=sample))
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


# rule contamination_mapping_bowtie2_PE:
#     input:
#         getAlignRawFastq
#     output:
#         "%s/bowtie2_contam_panel/{sample}/{panel}.sam" % RES_PATH,
#         "%s/bowtie2_contam_panel/{sample}/{panel}.unmapped.1.fq" % RES_PATH,
#         "%s/bowtie2_contam_panel/{sample}/{panel}.unmapped.2.fq" % RES_PATH,
#     params:
#         index=lambda wildcards: config['contamination_panel'][wildcards.panel],
#         # gz_support=lambda wildcards, input: "--un-gz" if str(input[0]).endswith('.gz') else "",
#         res_path=RES_PATH
#     threads: star_threads
#     message: "ALIGN: Align {wildcards.sample} to contamination panel by bowtie2"
#     log:
#         "%s/logs/bowtie2_contam_panel/{sample}_align_bowtie2_contam_{panel}.log" % RES_PATH
#     shell:
#         """
#         bowtie2 -x {params.index} -p {threads}  -1 {input[0]} -2 {input[1]} \
#         --un-conc {params.res_path}/bowtie2_contam_panel/{wildcards.sample}/{wildcards.panel}.unmapped.fq \
#         -S {params.res_path}/bowtie2_contam_panel/{wildcards.sample}/{wildcards.panel}.sam > {log} 2>&1
#         """


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



rule align_STAR:
    input:
        getAlignFastq
    output:
        "%s/STAR/{sample}/{sample}.Aligned.out.bam" % RES_PATH, 
        "%s/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam" % RES_PATH,
        "%s/STAR/{sample}/{sample}.Aligned.toTranscriptome.out.bam" % RES_PATH,
        "%s/STAR/{sample}/{sample}.ReadsPerGene.out.tab" % RES_PATH,
        "%s/STAR/{sample}/{sample}.Chimeric.out.junction" % RES_PATH,
        # "%s/STAR/{sample}/{sample}.Chimeric.out.sam" % RES_PATH,
        "%s/STAR/{sample}/{sample}.SJ.out.tab" % RES_PATH,
        "%s/STAR/{sample}/{sample}.Log.final.out" % RES_PATH
    params:
        gz_support=lambda wildcards, input: "--readFilesCommand zcat" if str(input[0]).endswith('.gz') else "",
        prefix=lambda wildcards: "{res_path}/STAR/{sample}/{sample}.".format(res_path=RES_PATH, sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample),
        # keepPairs = _keepPairs
    threads: star_threads
    message: "ALIGN: Align {wildcards.sample} to the genome by STAR"
    log:
        "%s/logs/STAR/{sample}.align_STAR.log" % RES_PATH
    shell:
        "STAR --runThreadN {threads} "
        "--genomeDir {config[STAR_index]} "
        "--genomeLoad NoSharedMemory "
        "--readFilesIn {input} "
        "{params.gz_support} "
        "--outSAMtype BAM Unsorted SortedByCoordinate "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMheaderHD @HD VN:1.4 "
        "--outSAMstrandField intronMotif "
        "--outSAMunmapped Within "
        "--outReadsUnmapped Fastx "
        "--outWigType bedGraph read1_5p "
        "--alignSJDBoverhangMin 10 "
        "--alignMatesGapMax 1000000 "
        "--alignIntronMax 1000000 "
        "--alignSJstitchMismatchNmax 5 -1 5 5 "
        "--chimSegmentMin 12 "
        "--chimJunctionOverhangMin 12 "
        # "--chimOutJunctionFormat 1 "
        # "--chimMultimapScoreRange 10 "
        # "--chimMultimapNmax 10 "
        # "--chimNonchimScoreDropMin 10 "
        # "--peOverlapNbasesMin 12 "
        # "--peOverlapMMp 0.1 "
        "--twopassMode Basic "
        "--quantMode TranscriptomeSAM GeneCounts > {log} 2>&1"


rule align_STARIndexSortedBam:
    # """INDEX the {sample}.sorted.bam file"""
    input:
        "%s/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam" % RES_PATH
    output:
        "%s/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai" % RES_PATH
    message: "ALIGN: Indexing {wildcards.sample}.Aligned.sortedByCoord.out.bam"
    shell:
        "samtools index {input}"

# rule align_STARIndexTranscripBam:
#     # """INDEX the {sample}.sorted.bam file"""
#     input:
#         "%s/STAR/{sample}/{sample}.Aligned.toTranscriptome.out.bam" % RES_PATH
#     output:
#         "%s/STAR/{sample}/{sample}.Aligned.toTranscriptome.out.bam.bai" % RES_PATH
#     message: "ALIGN: Indexing {wildcards.sample}.Aligned.toTranscriptome.out.bam"
#     shell:
#         "samtools index {input}"

rule align_STARBamtoBw:
    input:
        bam="%s/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam" % RES_PATH,
        bai="%s/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai" % RES_PATH
    output:
        "%s/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bw" % RES_PATH
    params:
        effective_genome_size=config['effective_genome_size']
    threads: 8
    log:
        "%s/logs/STAR/{sample}.align_STARBamtoBw.log" % RES_PATH
    shell:
        "bamCoverage -b {input.bam} -o {output} -p {threads} "
        "--binSize 10 --normalizeUsing CPM --effectiveGenomeSize {params.effective_genome_size} > {log} 2>&1"


rule align_STARfeatureCountTranscript:
    input:
        bam="%s/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam" % RES_PATH,
        bai="%s/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai" % RES_PATH
    output:
        "%s/STAR/{sample}/{sample}.ReadsPerTranscript.byfeatureCounts.txt" % RES_PATH
    params:
        gtf=config['gtf'],
        paired=lambda wildcards: '-p --countReadPairs' if len(config['samples'][wildcards.sample]) == 2 else ''
    threads: 4
    log:
        "%s/logs/STAR/{sample}.align_STARfeatureCountTranscript.log" % RES_PATH
    shell:
        "featureCounts -T {threads} -t exon -g transcript_id {params.paired} -o {output} -a {params.gtf} {input.bam} > {log} 2>&1"

rule align_STARfeatureCountGeneID:
    input:
        bam="%s/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam" % RES_PATH,
        bai="%s/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai" % RES_PATH
    output:
        "%s/STAR/{sample}/{sample}.ReadsPerGeneID.byfeatureCounts.txt" % RES_PATH
    params:
        gtf=config['gtf'],
        paired=lambda wildcards: '-p --countReadPairs' if len(config['samples'][wildcards.sample]) == 2 else ''
    threads: 4
    log:
        "%s/logs/STAR/{sample}.align_STARfeatureCountGene.log" % RES_PATH
    shell:
        "featureCounts -T {threads} -t exon -g gene_id {params.paired} -o {output} -a {params.gtf} {input.bam} > {log} 2>&1"

# rule align_STARfeatureCountGeneName:
#     input:
#         bam="%s/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam" % RES_PATH,
#         bai="%s/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai" % RES_PATH
#     output:
#         "%s/STAR/{sample}/{sample}.ReadsPerGeneName.byfeatureCounts.txt" % RES_PATH
#     params:
#         gtf=config['gtf'],
#         paired=lambda wildcards: '-p --countReadPairs' if len(config['samples'][wildcards.sample]) == 2 else ''
#     threads: 4
#     log:
#         "%s/logs/STAR/{sample}.align_STARfeatureCountGene.log" % RES_PATH
#     shell:
#         "featureCounts -T {threads} -t gene -g gene_name {params.paired} -o {output} -a {params.gtf} {input.bam} > {log} 2>&1"

rule align_STARrsemCalTPM:
    input:
        bam="%s/STAR/{sample}/{sample}.Aligned.toTranscriptome.out.bam" % RES_PATH,
    output:
        "%s/STAR/{sample}/{sample}.rsem.genes.results" % RES_PATH,
        "%s/STAR/{sample}/{sample}.rsem.isoforms.results" % RES_PATH
    params:
        rsem_index=config['rsem_index'],
        paired=lambda wildcards: '--paired-end' if len(config['samples'][wildcards.sample]) == 2 else '',
        prefix=lambda wildcards: "{res_path}/STAR/{sample}/{sample}.rsem".format(res_path=RES_PATH, sample=wildcards.sample)
    message: "ALIGN: Calculating TPM of {wildcards.sample} by RSEM"
    threads: 16
    log:
        "%s/logs/STAR/{sample}.align_STARrsemCalTPM.log" % RES_PATH
    shell:
        "rsem-calculate-expression --alignments -p {threads} {params.paired} "
        "{input.bam} {params.rsem_index}rsem {params.prefix} > {log} 2>&1"



# rule align_STAR_sort_bam:
#     input:
#         "%s/STAR/{sample}/{sample}.unsorted.bam" % RES_PATH,
#     output:
#         sortedBAM = "%s/STAR/{sample}/{sample}.sorted.bam" % RES_PATH,
#     message: "ALIGN: Sorting {wildcards.sample}.unsorted.bam"
#     threads: 4
#     params:
#         prefix=lambda wildcards: "{res_path}/STAR/{sample}/".format(res_path=RES_PATH, sample=wildcards.sample),
#     shell:
#         "samtools sort -T {params.prefix}TMP -o {output.sortedBAM} -@ {threads} {input}"


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
