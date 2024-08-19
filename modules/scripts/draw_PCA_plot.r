#!/usr/bin/env Rscript

# usage: Rscript draw_PCA_plot.r -i <input matrix path> -d <input control sample> -r <input plot result path>
# I recommend using the R package ggplot2 and ggrepel to draw the PCA plot

# required_Packages = c("optparse", "ggplot2", "ggrepel")
# if(!all(required_Packages %in% installed.packages())){
#   if (!requireNamespace("BiocManager", quietly = TRUE)){
#     install.packages("BiocManager")
#   }
#   BiocManager::install(setdiff(required_Packages, installed.packages()))
# }

require(optparse)
require(ggplot2)
require(ggrepel)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default=FALSE,
              help = "Input matrix path"),
  # make_option(c("-t", "--type"), type = "character", default=FALSE,
  #             help = "Whether the input data is RNAseq or Microarray"),
  make_option(c("-d", "--design"), type = "character", default=FALSE,
              help = "Input control sample"),
  make_option(c("-r", "--result"), type = "character", default=FALSE,
              help = "Input plot result path")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


profile_path = opt$input # using tpm matrix or expression matrix
# express_type = opt$type # tell the funtion which is RNAseq or Microarray
design_path = opt$design # describe the subclass of each sample used for PCA
result_path = opt$result # save the figure

expre_matr <- as.matrix(read.table(profile_path, header=T, sep="\t", row.names=1))
num_sam <- dim(expre_matr)[2]
# if(express_type == "RNASeq"){
#   expre_matr=expre_matr[,3:num_sam]
#   expre_matr=matrix(as.numeric(as.vector(unlist(expre_matr))),ncol=num_sam-2)
# } else {
#   expre_matr=expre_matr
#   expre_matr=matrix(as.numeric(as.vector(unlist(expre_matr))),ncol=num_sam)
# }
expre_matr <- log2(expre_matr + 1)
pca_raw <- prcomp(t(expre_matr), center = TRUE, scale. = FALSE)
edata_pc_df <- as.data.frame(pca_raw$x)
design_matrix <- read.table(design_path, header = TRUE, sep = ",")
na <- rep(TRUE, nrow(design_matrix))
# for (i in 3:ncol(design_matrix)){
#     na=na & is.na(design_matrix[i])
# }
# design_matrix <- design_matrix[,!na]
edata_pc_df$class <- design_matrix$condition
pca_plot <- ggplot(edata_pc_df, aes(x = PC1, y = PC2, color = class)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10))
ggsave(result_path, pca_plot, width = 6, height = 5)





# profile_path = opt$input
# express_type = opt$type
# result_path = opt$result
# treat = opt$treatname
# control = opt$controlname
# treatsample = strsplit(opt$treat,',')[[1]]
# controlsample = strsplit(opt$control,',')[[1]]

# # profile_path  = "/Users/shixiaoying/Downloads/ll.txt" 
# expre_matr=as.matrix(read.table(profile_path, sep = "\t", header = TRUE,
#                                 stringsAsFactors = FALSE,
#                                 quote = "", fill = TRUE))
# # colnames(expre_matr) <- sapply(strsplit(colnames(expre_matr),"_"),function(x) x[1])

# design_matrix <- as.data.frame(c(treatsample,controlsample))
# colnames(design_matrix) <- 'sample'
# design_matrix$condition <- c(rep(control,nrow(design_matrix)))
# design_matrix$condition[match(treatsample, as.character(design_matrix$sample))] <- treat
# design_matrix$label_c <- c(rep(1,nrow(design_matrix)))
# design_matrix$label_c[match(treatsample,design_matrix$sample)] <- 2
# design_matrix <- design_matrix[match(colnames(expre_matr),as.character(design_matrix$sample)),]

# if(express_type == "RNASeq"){
#   expre_matr=expre_matr[,c(-1,-2)]  
# } else {
#   expre_matr=expre_matr[,(ncol(expre_matr)-nrow(design_matrix)+1):ncol(expre_matr)]
# }
# num_sam=dim(expre_matr)[2]
# expre_matr=matrix(as.numeric(as.vector(unlist(expre_matr))),ncol=num_sam)
# expre_matr=log2(expre_matr+1)
# pca_raw <- prcomp(t(expre_matr), center = TRUE, scale. = F)
# edata_pc_df <- as.data.frame(pca_raw$x)
# edata_pc_df$class=design_matrix$condition
# pca_plot=ggplot(edata_pc_df, aes(x = PC1, y = PC2, color = class)) +
#   geom_point()+theme_bw()+theme(axis.text.x = element_text(size = 10))+theme(axis.text.y = element_text(size = 10))
# ggsave(result_path, pca_plot, width = 6, height = 5)
