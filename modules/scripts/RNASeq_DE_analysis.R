# command: Rscript RNASeq_DE_analysis.R -i "/Users/shixiaoying/analysis/GSE110708/expression/GSE110708_Rawcount_matrix.txt" -r "/Users/shixiaoying/work2019/Test_DEAP/test_result/" -t "IFNy" -c "NT" --controlname "GSM3014871,GSM3014872,GSM3014873" --treatname "GSM3014868,GSM3014869,GSM3014870"

require(DESeq2)
require(optparse)

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
              help="Input treatname sample"),
  make_option(c("--controlname"), type="character", default=TRUE,
              help="Input control sample")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

profile_path=opt$input 
result_path=opt$result
treat = opt$treat
control = opt$control
treatsample = opt$treatname
controlsample = opt$controlname
treatsample = strsplit(treatsample,',')[[1]]
controlsample = strsplit(controlsample,',')[[1]]

rawcount=as.matrix(read.table(profile_path,header =T,sep="\t"))
tran_gene=rawcount[,c(1,2)]
rawcount=rawcount[,c(-1,-2)]
mode(rawcount)<-"integer"

design_matrix <- as.data.frame(c(treatsample,controlsample))
colnames(design_matrix)="sample"
design_matrix$condition <- c(rep(control,nrow(design_matrix)))
design_matrix$condition[match(treatsample, as.character(design_matrix$sample))] <- treat
design_matrix$label_c <- c(rep(1,nrow(design_matrix)))
design_matrix$label_c[match(treatsample,design_matrix$sample)] <- 2
design_matrix <- design_matrix[match(design_matrix$sample,colnames(rawcount)),]

##function used to performance specific treatment DSE
# get_DSE<-function(x){
#   con_sample=as.vector(unlist(design_matrix$sample[which(design_matrix[,]=="1")]))
#   case_sample=as.vector(unlist(design_matrix$sample[which(design_matrix[,x]=="2")]))
# compa_samp=c(con_sample,case_sample)
#   
if(length(treatsample) >= 2 & length(controlsample) >= 2){
  # comp_cond=design_matrix$treatment[match(compa_samp,design_matrix$sample)]
  # control=as.character(unlist(design_matrix$treatment[which(design_matrix[,3]=="1")][1]))
  # case=as.character(unlist(design_matrix$treatment[which(design_matrix[,3]=="2")][1]))
  # temp_rawcount=rawcount[,match(compa_samp,colnames(rawcount))] 
  #
  rawcount = rawcount[,match(design_matrix$sample,colnames(rawcount))]
  sampleTable <- data.frame(sampleName = design_matrix$sample,  condition = design_matrix$condition)
  sampleTable$condition <- factor(sampleTable$condition)
  sampleTable$condition <- relevel(sampleTable$condition, ref = control)
  ddsFullCountTable <- DESeqDataSetFromMatrix(countData = rawcount,colData = sampleTable,  design= ~ condition)
  dds <- DESeq(ddsFullCountTable) 
  des.re=results(dds)
  df <- data.frame(matrix(unlist(des.re@listData), ncol=6))
  colnames(df)=names(des.re@listData)
  df=cbind(tran_gene,df)
  colnames(df)=c('refseq','symbol','AveExpr','logFC','lfcSE','stat','P.Value','adj.P.Val')
  write.table(df,result_path,quote=F,row.names = F,sep="\t")
}

