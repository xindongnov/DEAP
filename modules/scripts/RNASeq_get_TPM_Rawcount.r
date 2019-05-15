required_Packages = c("DESeq2","ggplot2","ggrepel")

if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}

require(DESeq2)
require(ggplot2)
require(ggrepel)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
rm(list = ls())

args = commandArgs(T)
file_path = args[1] #data file path
design_path = args[2] #design matrix path
result_path = args[3] #result path
species = args[4] #cdf names

#DESeq calculate different expression gene on raw count

#design_path
design_matrix=read.table(design_path,header = T,sep="\t")
#now the colnames of design matrix is sample,compare1,compare2,compare3
#file_path="/Volumes/Macintosh HD/Users/yahan/Downloads/Microarray_data"
#file_path folder which is including all sample data in one dataset 
sample_path=list.dirs(file_path,recursive=T)#each sample folder alseo including many folder, only one is the result of alignment
sample_path=sample_path[grep("align",sample_path)]
tpm_matrix=data.frame(Name=NA)
rawcount=data.frame(Name=NA)
get_gene<-function(x){
  temp_path=paste0(x,"/","quant.sf")
  temp_matrix=read.table(temp_path,header=T,sep = "\t")
  temp_matrix_tpm=temp_matrix[,c("Name","TPM")]
  temp_matrix_count=temp_matrix[,c("Name","NumReads")]
  tpm_matrix<<-merge(tpm_matrix,temp_matrix_tpm,by="Name",all=TRUE)
  rawcount<<-merge(rawcount,temp_matrix_count,by="Name",all=TRUE)
}
result=apply(matrix(sample_path),1,get_gene)

sample=unlist(lapply(strsplit(sample_path,"/"),function(x) tail(x,2)[1]))
tran_tpm=unlist(lapply(strsplit(tpm_matrix[,1],"\\."),function(x) x[1]))
tran_raw=unlist(lapply(strsplit(rawcount[,1],"\\."),function(x) x[1]))
#transform the ENSG to GENE ID and symbol



if(species="Mouse"){
  tran_gene_tpm=select(org.Mm.eg.db,keys=tran_tpm,columns = c("SYMBOL"), keytype="ENSEMBLTRANS")
}

if(species="Human"){
  tran_gene_raw=select(org.Hs.eg.db,keys=tran_raw,columns = c("SYMBOL"), keytype="ENSEMBLTRANS")#
}
tran_gene_tpm=tran_gene_tpm[match(tran_tpm,tran_gene_tpm$ENSEMBLTRANS),]
tpm_matrix=cbind(tran_gene_tpm,tpm_matrix[-1])
tran_gene_raw=tran_gene_raw[match(tran_raw,tran_gene_raw$ENSEMBLTRANS),]
rawcount=cbind(tran_gene_raw,rawcount[-1])
tpm_matrix=tpm_matrix[which(!is.na(tpm_matrix[,1])),]
rawcount=rawcount[which(!is.na(rawcount[,1])),]
colnames(rawcount)=c("ENSEMBLTRANS","SYMBOL",sample)
colnames(tpm_matrix)=c("ENSEMBLTRANS","SYMBOL",sample)
write.table(tpm_matrix,paste0(result_path,"/","TPM_matrix.txt"),col.names=T,row.names=F,sep="\t",quote = F)
write.table(rawcount,paste0(result_path,"/","RawCount_matrix.txt"),col.names=T,row.names=F,sep="\t",quote = F)