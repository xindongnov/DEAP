# command :Rscript /Users/shixiaoying/work2019/Test_DEAP/DEAP_scripts/RNASeq_get_Rawcount.R -i "/Users/shixiaoying/analysis/GSE110708" -r "/Users/shixiaoying/work2019/Test_DEAP/test_result/c.txt" -s "Mouse"

# ================================
# @auther: Ya Han
# @date: Sep 2019
# ================================

required_Packages = c("optparse","org.Mm.eg.db","org.Hs.eg.db")
if(!all(required_Packages %in% installed.packages())){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install(setdiff(required_Packages, installed.packages()))
}

require(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default=FALSE,
              help="Input sample path folder"),
  make_option(c("-n", "--name"), type = "character", default=FALSE,
              help="Input samples name"),
  make_option(c("-o", "--output"), type="character", default=FALSE,
              help="Input count result path"),
  make_option(c("-s", "--species"), type="character", default=FALSE,
              help="Input the species: Mouse or Human")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

species = opt$species
files = strsplit(opt$input,",")[[1]]
result_path_COUNT = opt$output
names = strsplit(opt$name,",")[[1]]

#DESeq calculate different expression gene on raw count

#design_path
#design_matrix=read.table(design_path,header = T,sep="\t")
#now the colnames of design matrix is sample,compare1,compare2,compare3
#file_path="/Users/shixiaoying/analysis/GSE110708"
#file_path folder which is including all sample data in one dataset 

# sample_path=list.dirs(file_path,recursive=T)#each sample folder also including many folder, only one is the result of alignment
# sample_path=sample_path[grep("align$",sample_path)]
rawcount=data.frame(Name=NA)
get_gene<-function(x){
  # x=paste0(x,"/","quant.sf")
  temp_matrix=read.table(x,header=T,sep = "\t")
  temp_matrix_count=temp_matrix[,c("Name","NumReads")]
  rawcount<<-merge(rawcount,temp_matrix_count,by="Name",all=TRUE)
}
result=apply(matrix(files),1,get_gene)

# sample=unlist(lapply(strsplit(sample_path,"/"),function(x) tail(x,2)[1]))
tran_raw=unlist(lapply(strsplit(rawcount[,1],"\\."),function(x) x[1]))

#transform the ENSG to GENE ID and symbol

if(species == "Mouse"){
  library(org.Mm.eg.db)
  tran_gene_raw=select(org.Mm.eg.db,keys=tran_raw,columns = c("SYMBOL"), keytype="REFSEQ")# 
  }

if(species == "Human"){
  library(org.Hs.eg.db)
  tran_gene_raw=select(org.Hs.eg.db,keys=tran_raw,columns = c("SYMBOL"), keytype="REFSEQ")#
  
}
tran_gene_raw=tran_gene_raw[match(tran_raw,tran_gene_raw$REFSEQ),]
rawcount=cbind(tran_gene_raw,rawcount[-1])
rawcount=rawcount[which(!is.na(rawcount[,1])),]
colnames(rawcount)=c("refseq","symbol",names)

write.table(rawcount,result_path_COUNT,col.names=T,row.names=F,sep="\t",quote = F)
