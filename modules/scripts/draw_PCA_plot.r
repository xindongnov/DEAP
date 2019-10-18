#command:Rscript /Users/shixiaoying/work2019/Test_DEAP/DEAP_scripts/draw_PCA_plot.r -i "/Users/shixiaoying/Downloads/ll.txt"  -k "Microarray" -c 'control' -t treat --controlname  "GSM1017442,GSM1017443,GSM1017444,GSM1017445" --treatname 'GSM1017446,GSM1017447' -r /Users/shixiaoying/work2019/Test_DEAP/test_result/pca.pdf

require(optparse)
require(ggplot2)
require(ggrepel)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default=FALSE,
              help="Input matrix path"),
  make_option(c("-k", "--type"), type = "character", default=FALSE,
              help="Whether the input data is RNAseq or Microarray"),
  make_option(c("-t", "--treat"), type="character", default=TRUE,
              help="Input treatment"),
  make_option(c("-c", "--control"), type="character", default=TRUE,
              help="Input control"),
  make_option(c("--treatname"), type="character", default=TRUE,
              help="Input treatname sample"),
  make_option(c("--controlname"), type="character", default=TRUE,
              help="Input control sample"),
  make_option(c("-r", "--result"), type="character", default=FALSE,
              help="Input plot result path")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

profile_path = opt$input
express_type = opt$type
result_path = opt$result
treat = opt$treat
control = opt$control
treatsample = opt$treatname
controlsample = opt$controlname
treatsample = strsplit(treatsample,',')[[1]]  #"GSM1017442,GSM1017443,GSM1017444,GSM1017445"
controlsample = strsplit(controlsample,',')[[1]] #'GSM1017446,GSM1017447' 


# profile_path  = "/Users/shixiaoying/Downloads/ll.txt" 
expre_matr=as.matrix(read.table(profile_path, sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE,
                                quote = "", fill = TRUE))
# colnames(expre_matr) <- sapply(strsplit(colnames(expre_matr),"_"),function(x) x[1])

design_matrix <- as.data.frame(c(treatsample,controlsample))
colnames(design_matrix) <- 'sample'
design_matrix$treatment <- c(rep(control,nrow(design_matrix)))
design_matrix$treatment[match(treatsample, as.character(design_matrix$sample))] <- treat
design_matrix$label_c <- c(rep(1,nrow(design_matrix)))
design_matrix$label_c[match(treatsample,design_matrix$sample)] <- 2
design_matrix <- design_matrix[match(colnames(expre_matr),as.character(design_matrix$sample)),]

if(express_type == "RNASeq"){
  expre_matr=expre_matr[,c(-1,-2)]  
} else {
  expre_matr=expre_matr[,(ncol(expre_matr)-nrow(design_matrix)+1):ncol(expre_matr)]
}
num_sam=dim(expre_matr)[2]
expre_matr=matrix(as.numeric(as.vector(unlist(expre_matr))),ncol=num_sam)
expre_matr=log2(expre_matr+1)
pca_raw <- prcomp(t(expre_matr), center = TRUE, scale. = F)
edata_pc_df <- as.data.frame(pca_raw$x)
edata_pc_df$class=design_matrix$treatment
pca_plot=ggplot(edata_pc_df, aes(x = PC1, y = PC2, color = class)) +
  geom_point()+theme_bw()+theme(axis.text.x = element_text(size = 10))+theme(axis.text.y = element_text(size = 10))
ggsave(result_path, pca_plot, width = 6, height = 5)
