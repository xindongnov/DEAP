
# ================================
# @auther: Ya Han
# @date: Sep 2019
# ================================

required_Packages = c("DESeq2","ggplot2","ggrepel")

if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}

require(DESeq2)
require(ggplot2)
require(ggrepel)
rm(list = ls())

args = commandArgs(T)
file_path = args[1] #data file path
design_path = args[2] #design matrix path
result_path = args[3] #result path
species = args[4] #cdf names
GPL_path=arg[5]#GPL information path
FC_setting = args[6] #setting fold change
pval_setting = args[7] #setting adj.p.val
FC_setting <- as.numeric(FC_setting)
pval_setting <- as.numeric(pval_setting)
logFC_setting <- log2(FC_setting)

#DESeq 是直接

library(readr)
design_matrix=read.table(file_path,header = F,sep="\t")
colnames(design_matrix)=c("sample","condition")
rawcount=read.table(file_path,header = T,row.names = 1,sep="\t")
rawcount=rawcount[,na.omit(match(design_matrix$sample,colnames(rawcount)))]
design_matrix <- data.frame(lapply(design_matrix, function(x) factor(x, levels=unique(x))))
library(DESeq2)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = rawcount,colData = design_matrix,  design= ~ condition)
dds <- DESeq(ddsFullCountTable) 

#draw a TPM plot and export the TPM expression profile
getTPM <- function(countMat, idType = "Ensembl"){
  ensembl = read.csv("/Users/yahan/Documents/SPARCS/SPARCS_TCGA/gdac.broadinstitute.org_SKCM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/ensembl_hg38_reference.txt",check.names = FALSE)
  ensembl$Length <- abs(ensembl$`Gene end (bp)` - ensembl$`Gene start (bp)`)
  if(toupper(idType) == "ENSEMBL"){
    len <- ensembl[match(rownames(countMat),ensembl$`Gene stable ID`), "Length"]
    rownames(countMat) = ensembl[match(rownames(countMat),ensembl$`Gene stable ID`), "Gene name"]
  }
  else if(toupper(idType) == "SYMBOL")
    len <- ensembl[match(rownames(countMat),ensembl$`Gene name`), "Length"]
  else
    stop("Please input right type of gene name, such as Ensembl or gene Symbol ...")

  na_idx = which(is.na(len))
  if(length(na_idx)>0){
    warning(paste0("Omit ", length(na_idx), " genes of which length is not available !"))
    countMat = countMat[!is.na(len),]
    len = len[!is.na(len)]
  }
  tmp <- countMat / len
  TPM <- 1e6 * t(t(tmp) / colSums(tmp))
  TPM = TPM[!duplicated(rownames(TPM)),]
  return(TPM)
}
tpm=getTPM(rawcount,"SYMBOL")
pca_raw <- prcomp(t(tpm), center = TRUE, scale. = F)
edata_pc_df <- as.data.frame(pca_raw$x)
edata_pc_df$class=design_matrix$condition
pca_plot=ggplot(edata_pc_df, aes(x = PC1, y = PC2, color = class)) +
    geom_point()
ggsave(paste0(result_path,"/PCA_plot.png"), pca_plot, width = 5, height = 5)
rm(tpm)

#calculate the different expression gene
condition_matrix=matrix(as.vector(unlist(t(combn(unique(design_matrix$condition),2)))),ncol=2)
condition_DEG<-function(x){
  contrastV <- c("TY", x[1], x[2])
  rescond <- results(dds,  contrast=contrastV)
  write.table(rescond,file=paste0(result_path,"/",x[1],"_",x[2],"DESeq_table.txt"),quote=F,sep="\t")
  
  result=matrix(unlist(rescond@listData),ncol=6)
  colnames(result)=names(rescond@listData)
  rownames(result)=rescond@rownames
  result=data.frame(result)
  
  #select significant up-regulated and down-regulated gene to draw Volcano plot
  sigout <- subset(result, abs(log2FoldChange) > logFC_setting & padj < pval_setting)
  upR <- sigout[which(sigout$log2FoldChange > 0), 0]
  downR <- sigout[which(sigout$log2FoldChange < 0), 0]

volca <- ggplot(result,aes(log2FoldChange,-log10(padj))) + 
  geom_point(aes(color=ifelse((log2FoldChange > logFC_setting | log2FoldChange < -logFC_setting) & padj < pval_setting, "sign","non-sign")), cex = 1, alpha = 0.3) + 
  scale_colour_manual("Significance", values = c("steelblue","red")) +
  labs(title = 'Volcano plot', x = "log2 Fold Change", y = "-log10 adj.P.Val")+
  geom_text_repel(data = sigout, aes(label = rownames(sigout)), cex = 3,box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))
ggsave(paste0(result_path,'/',x[1],"_",x[2],"Volcano_plot.png"), volca, width = 7, height = 7)
}
apply(condition_matrix,1,condition_DEG)
