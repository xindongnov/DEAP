
required_Packages = c("DESeq2")

if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}

library(DESeq2)
args = commandArgs(T)
profile_path=args[1]
design_path=args[2]
result_path=args[3]
compare_output_path=args[4]

#profile_path="/Users/yahan/Downloads/Microarray_data/RawCount_matrix.txt"
#design_path="/Users/yahan/Desktop/design_matrix.txt"

rawcount=as.matrix(read.table(profile_path,header =T,sep="\t"))
tran_gene=rawcount[,c(1,2)]
rownames(rawcount)=rawcount[,1]
rawcount=rawcount[,c(-1,-2)]
sample=colnames(rawcount)

rawcount=matrix(as.numeric(as.vector(unlist(rawcount))),ncol=length(sample))
colnames(rawcount)=sample
mode(rawcount)<-"integer"
design_matrix=read.table(design_path,header = T,sep=",")

num_com=dim(design_matrix)[2]

##function used to performance specific treatment DSE
get_DSE<-function(x){
  con_sample=as.vector(unlist(design_matrix$sample[which(design_matrix[,x]=="1")]))
  case_sample=as.vector(unlist(design_matrix$sample[which(design_matrix[,x]=="2")]))
  compa_samp=c(con_sample,case_sample)
  
  if(length(case_sample) >= 2 & length(con_sample) >= 2){
    comp_cond=design_matrix$treatment[match(compa_samp,design_matrix$sample)]
    contral=as.character(unlist(design_matrix$treatment[which(design_matrix[,x]=="1")][1]))
    case=as.character(unlist(design_matrix$treatment[which(design_matrix[,x]=="2")][1]))
    temp_rawcount=rawcount[,match(compa_samp,colnames(rawcount))] 
    
    sampleTable <- data.frame(sampleName = compa_samp,  condition = comp_cond)
    sampleTable$condition <- factor(sampleTable$condition)
    sampleTable$condition <- relevel(sampleTable$condition, ref = contral)
    ddsFullCountTable <- DESeqDataSetFromMatrix(countData = temp_rawcount,colData = sampleTable,  design= ~ condition)
    dds <- DESeq(ddsFullCountTable) 
    des.re=results(dds)
    df <- data.frame(matrix(unlist(des.re@listData), ncol=6))
    colnames(df)=names(des.re@listData)
    df=cbind(tran_gene,df)
    write.table(df,file=paste0(result_path,"/",case,"_",contral,"_","DESeq_table.txt"),quote=F,row.names = F,sep="\t")
  return(c(case,contral))}
}
result=apply(matrix(seq(3:num_com)),1,get_DSE)
write.table(result,compare_output_path,quote=F,,col.names=c("case","contral"),row.names = F,sep="\t")

