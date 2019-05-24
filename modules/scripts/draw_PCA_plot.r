
required_Packages = c("ggplot2","ggrepel")

if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}

library(ggplot2)
library(ggrepel)
args = commandArgs(T)
profile_path=args[1] #which is used to save the tpm matrix or expression matrix
express_type=args[2] # tell the funtion which is RNAseq or Microarray
design_path=args[3] #describe the subclass of each sample used for PCA

result_path=args[4] #save the figure

expre_matr=as.matrix(read.table(profile_path,header =T,sep="\t"))
if(express_type == "RNASeq"){
  expre_matr=expre_matr[,c(-1,-2)]
} else {
  expre_matr=expre_matr[,-1]
}
num_sam=dim(expre_matr)[2]
expre_matr=matrix(as.numeric(as.vector(unlist(expre_matr))),ncol=num_sam)
expre_matr=log2(expre_matr)
pca_raw <- prcomp(t(expre_matr), center = TRUE, scale. = F)
edata_pc_df <- as.data.frame(pca_raw$x)
design_matrix=read.table(design_path,header = T,sep=",")
edata_pc_df$class=design_matrix$treatment
pca_plot=ggplot(edata_pc_df, aes(x = PC1, y = PC2, color = class)) +
  geom_point()+theme_bw()+theme(axis.text.x = element_text(size = 10))+theme(axis.text.y = element_text(size = 10))
ggsave(result_path, pca_plot, width = 6, height = 5)
