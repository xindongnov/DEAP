required_Packages = c("optparse","ggplot2","ggrepel")
if(!all(required_Packages %in% installed.packages())){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install(setdiff(required_Packages, installed.packages()))
}


require(optparse)
require(ggplot2)
require(ggrepel)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default=FALSE,
              help="Input matrix path"),
  make_option(c("-r", "--result"), type="character", default=FALSE,
              help="Input plot result path"),
  make_option(c("-c", "--foldchange"), type="character", default='0',
              help="Set fold change value"),
  make_option(c("-f", "--fdr"), type="character", default='1',
              help="Set adj.p.val")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input_path = opt$input
result_path = opt$result
FC_setting = opt$foldchange
fdr_setting = opt$fdr
FC_setting <- as.numeric(FC_setting)
fdr_setting <- as.numeric(fdr_setting)
logFC_setting <- log2(FC_setting)

# file_path="/Users/shixiaoying/Downloads/result_path.txt"
file=read.table(input_path,sep = "\t", header = TRUE,
                stringsAsFactors = FALSE,
                quote = "", fill = TRUE)
sigout <- subset(file, abs(logFC) > logFC_setting & adj.P.Val < fdr_setting)
index=order(sigout$logFC)
if(length(index) > 20){
sigout=sigout[c(index[1:10],rev(index)[1:10]),]}

upR <- as.vector(unlist(sigout$symbol[which(sigout$logFC > 0)]))
downR <- as.vector(unlist(sigout$symbol[which(sigout$logFC < 0)]))

volca <- ggplot(file,aes(logFC,-log10(adj.P.Val))) + 
  geom_point(aes(color=ifelse((logFC > logFC_setting | logFC < -logFC_setting) & adj.P.Val < fdr_setting, "sign","non-sign")), cex = 1, alpha = 0.3) + 
  scale_colour_manual("Significance", values = c("grey","red")) +
  labs(title = 'Volcano plot', x = "log2 Fold Change", y = "-log10 adj.P.Val")+
  geom_text_repel(data = sigout, aes(label = sigout$symbol), cex = 3,box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))+
  theme_bw()+geom_hline(aes(yintercept=-log10(fdr_setting)),linetype="dashed",color="grey",size=0.5)+
  geom_vline(aes(xintercept=logFC_setting),linetype="dashed",color="grey",size=0.5)+
  geom_vline(aes(xintercept=-logFC_setting),linetype="dashed",color="grey",size=0.5)+theme(axis.text.x = element_text(size = 10))+theme(axis.text.y = element_text(size = 10))

ggsave(result_path, volca, width = 7, height = 7)
