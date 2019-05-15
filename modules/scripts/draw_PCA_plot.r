library(ggplot2)
library(ggrepel)

profile_path #which is used to save the tpm matrix or expression matrix
express_type # tell the funtion which is RNAseq or Microarray
design_path #describe the subclass of each sample used for PCA
result_path #save the figure

expre_matr=as.matrix(read.table(profile_path,header =T,sep="\t"))
if(express_type == "RNASeq"){
  expre_matr=expre_matr[,c(-1,-2)]
} else {
  expre_matr=expre_matr[,-1]
}
num_sam=dim(expre_matr)[2]
expre_matr=matrix(as.numeric(as.vector(unlist(expre_matr))),ncol=num_sam)
pca_raw <- prcomp(t(expre_matr), center = TRUE, scale. = F)
edata_pc_df <- as.data.frame(pca_raw$x)
design_matrix=read.table(design_path,header = T,sep="\t")
edata_pc_df$class=design_matrix$condition
pca_plot=ggplot(edata_pc_df, aes(x = PC1, y = PC2, color = class)) +
  geom_point()
ggsave(paste0(result_path,"/PCA_plot.png"), pca_plot, width = 5, height = 5)
