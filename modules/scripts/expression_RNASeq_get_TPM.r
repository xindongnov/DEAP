#command: Rscript /Users/shixiaoying/work2019/Test_DEAP/DEAP_scripts/RNASeq_get_TPM.R -i "/Users/shixiaoying/analysis/GSE110708" -r "/Users/shixiaoying/work2019/Test_DEAP/test_result/t.txt" -s "Mouse"

# ================================
# @auther: Ya Han
# @date: Sep 2019
# ================================

required_Packages = c("optparse")
if(!all(required_Packages %in% installed.packages())){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install(setdiff(required_Packages, installed.packages()))
}

require(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default=FALSE, help="Input quant.sf files, use comma to split"),
  make_option(c("-n", "--name"), type = "character", default=FALSE, help="Input samples name"),
  make_option(c("-o", "--output"), type="character", default=FALSE, help="Output of TPM result path"),
  make_option(c("-s", "--species"), type="character", default=FALSE, help="Input the species: Mouse or Human")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

species = opt$species
files = strsplit(opt$input,",")[[1]]
result_path_TPM = opt$output
names = strsplit(opt$name,",")[[1]]

# sample_path=list.dirs(file_path,recursive=T)
# sample_path=sample_path[grep("align$",sample_path)]
tpm_matrix=data.frame(Name=NA)
get_gene<-function(x){
  # x=paste0(x,"/","quant.sf")
  temp_matrix=read.table(x,header=T,sep = "\t")
  temp_matrix_tpm=temp_matrix[,c("Name","TPM")]
  tpm_matrix<<-merge(tpm_matrix,temp_matrix_tpm,by="Name",all=TRUE)
}
result=apply(matrix(files),1,get_gene)

# sample=unlist(lapply(strsplit(sample_path,"/"),function(x) tail(x,2)[1]))
tran_tpm=unlist(lapply(strsplit(tpm_matrix[,1],"\\."),function(x) x[1]))

#transform the ENSG to GENE ID and symbol

if(species == "Mouse"){
  library(org.Mm.eg.db)
  tran_gene_tpm=select(org.Mm.eg.db,keys=tran_tpm,columns = c("SYMBOL"), keytype="REFSEQ")
}

if(species == "Human"){
  library(org.Hs.eg.db)
  tran_gene_tpm=select(org.Hs.eg.db,keys=tran_tpm,columns = c("SYMBOL"), keytype="REFSEQ")

}

tran_gene_tpm=tran_gene_tpm[match(tran_tpm,tran_gene_tpm$REFSEQ),]
tpm_matrix=cbind(tran_gene_tpm,tpm_matrix[-1])
tpm_matrix=tpm_matrix[which(!is.na(tpm_matrix[,1])),]
colnames(tpm_matrix)=c("refseq","symbol",names)
write.table(tpm_matrix,result_path_TPM,col.names=T,row.names=F,sep="\t",quote = F)
