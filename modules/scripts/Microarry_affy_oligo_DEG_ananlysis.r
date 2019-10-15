
# ================================
# @auther: Ya Han
# @email: xindong9511@gmail.com
# @date: Sep 2019
# ================================

# required_Packages = c("limma")
# if(!all(required_Packages %in% installed.packages())){
#   source("https://bioconductor.org/biocLite.R")
#   biocLite(setdiff(required_Packages, installed.packages()))
# }
require(limma)
args = commandArgs(T)
file_path = args[1] #data file path
design_path = args[2] #design matrix path
result_path=args[3]
result_path_com=args[4]
FC_setting = args[5] #setting fold change
pval_setting = args[6] #setting adj.p.val
platform = args[7]
FC_setting <- as.numeric(FC_setting)
pval_setting <- as.numeric(pval_setting)
logFC_setting <- log2(FC_setting)

#logFC_setting=log2(2)
#pval_setting=0.01
#make sure the same sample order in design_matrix and expression profi
expres_matr=read.table(file_path,header=T,sep="\t",fill=T,row.names = 1,quote = "")
#design_path <- "/Users/yahan/Downloads/GSE12056_RAW/gse12056_design_matrix.txt"
design_mat = read.table(design_path, sep=",", header=T)
design_mat$treatment=gsub("-"," ",design_mat$treatment)#- is a sepcific signal in limma contrast_matrix
design_mat$treatment=gsub(" ",".",design_mat$treatment)

calcu_DEG<-function(x){
  con_sam=as.vector(unlist(design_mat$sample[which(design_mat[,x]=="1")]))
  cas_sam=as.vector(unlist(design_mat$sample[which(design_mat[,x]=="2")]))
  if(length(con_sam) >=2 & length(cas_sam) >= 2){
  control=as.vector(unlist(design_mat$treatment[which(design_mat[,x]=="1")]))[1]
  case=as.vector(unlist(design_mat$treatment[which(design_mat[,x]=="2")]))[1]
  sam=c(con_sam,cas_sam)
  tem_mat=expres_matr[,sam]
  grouplist = as.vector(design_mat[which(!is.na(design_mat[,x])),2])
  # read design
  design <- model.matrix(~0+factor(grouplist))
  colnames(design) = levels(factor(grouplist))
  rownames(design) = colnames(tem_mat)
  fit <- lmFit(tem_mat,design)
  cont.matrix <- makeContrasts(contrasts = paste0(case,"-",control),levels = colnames(fit$coefficients))
  fit <- contrasts.fit(fit, cont.matrix)
  fit <- eBayes(fit)
  output <- topTable(fit, number = Inf, lfc = 0, p.value = 1)
  output<-cbind(rownames(output),output)
  colnames(output) <-c("SYMBOL","log2FoldChange","AveExpr","t","P.Value","padj","B")
  # output file
  write.table(output,paste0(result_path,'/',case,"_",control,'_limma_table.txt'),sep = '\t',quote = F,col.names = T,row.names = F)
  return(c(case,control))}
}
result=apply(matrix(3:dim(design_mat)[2]),1,calcu_DEG)
if(length(unlist(result)) >= 2){
  result=matrix(unlist(result),ncol=2)
colnames(result)=c("case","control")
write.table(result,result_path_com,quote=F,row.names = F,sep="\t")}
