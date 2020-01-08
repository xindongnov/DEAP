#!/usr/bin/env Rscript

# command: Rscript RNASeq_DE_analysis.R -i "/Users/shixiaoying/analysis/GSE110708/expression/GSE110708_Rawcount_matrix.txt" -r "/Users/shixiaoying/work2019/Test_DEAP/test_result/" -t "IFNy" -c "NT" --controlname "GSM3014871,GSM3014872,GSM3014873" --treatname "GSM3014868,GSM3014869,GSM3014870"

required_Packages = c("DESeq2", "optparse", "sva")
if(!all(required_Packages %in% installed.packages())){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install(setdiff(required_Packages, installed.packages()))
}

require(DESeq2)
require(optparse)
require(sva)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default=FALSE,
              help="Input matrix path"),
  make_option(c("-r", "--result"), type="character", default=FALSE,
              help="Input differential expression result path"),
  make_option(c("-t", "--treat"), type="character", default=TRUE,
              help="Input treatment"),
  make_option(c("-c", "--control"), type="character", default=TRUE,
              help="Input control"),
  make_option(c("--treatname"), type="character", default=TRUE,
              help="Input treat sample"),
  make_option(c("--controlname"), type="character", default=TRUE,
              help="Input control sample")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

profile_path=opt$input 
result_path=opt$result
treatname = opt$treatname
controlname = opt$controlname
treatsample = strsplit(opt$treat,',')[[1]]
controlsample = strsplit(opt$control,',')[[1]]

rawcount=as.matrix(read.table(profile_path,header =T,sep="\t"))
tran_gene=rawcount[,c(1,2)]
rawcount=rawcount[,c(-1,-2)]
rownames(rawcount) <- tran_gene[,1]
mode(rawcount)<-"integer"
design_matrix <- as.data.frame(c(treatsample,controlsample))
colnames(design_matrix)="sample"
design_matrix$condition <- c(rep(controlname,nrow(design_matrix)))
design_matrix$condition[match(treatsample, as.character(design_matrix$sample))] <- treatname
design_matrix$label_c <- c(rep(1,nrow(design_matrix)))
design_matrix$label_c[match(treatsample,design_matrix$sample)] <- 2
# print(design_matrix)
# print(match(design_matrix$sample,colnames(rawcount)),])
# design_matrix <- design_matrix[match(design_matrix$sample,colnames(rawcount)),]

if(length(treatsample) >= 2 & length(controlsample) >= 2){
  rawcount = rawcount[,match(design_matrix$sample,colnames(rawcount))]
  sampleTable <- data.frame(sampleName = design_matrix$sample,  condition = design_matrix$condition)
  sampleTable$condition <- factor(sampleTable$condition)
  sampleTable$condition <- relevel(sampleTable$condition, ref = controlname)

  ddsFullCountTable <- DESeqDataSetFromMatrix(countData = rawcount, colData = sampleTable, design= ~ condition)
  dds <- DESeq(ddsFullCountTable)
  # dds <- dds[rowSums(counts(dds)) > 0,]
  dat <- counts(dds, normalized=TRUE)
  idx <- rowMeans(dat) > 1 
  dat <- dat[idx,]
  mod <- model.matrix(~ design_matrix$label_c, colData(dds))
  mod0 <- model.matrix(~ 1, colData(dds))
  # n.sv <- num.sv(dat,mod,method="leek")
  svseq <- svaseq(dat, mod, mod0, n.sv=1)
  ddssva <- dds
  ddssva$SV1 <- svseq$sv[,1]
  design(ddssva) <- ~ SV1 + condition
  ddssva <- DESeq(ddssva)
  des.re <- results(ddssva)

  df <- data.frame(matrix(unlist(des.re@listData), ncol=6))
  colnames(df)=names(des.re@listData)
  df <- cbind(tran_gene,df)
  colnames(df)=c('refseq','symbol','AveExpr','logFC','lfcSE','stat','P.Value','adj.P.Val')
  df_dropna <- na.omit(df)
  # aggregate(df_dropna$AveExpr, by = list(df_dropna$symbol), max)
  df_dropna_dedup <- df_dropna[!duplicated(df_dropna$symbol),]
  write.table(df_dropna_dedup, result_path, quote=F, row.names = F, sep="\t")
}

