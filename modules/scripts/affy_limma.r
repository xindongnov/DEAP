required_Packages = c("affy", "limma","sva","ggplot2","ggrepel")

if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}

require(affy)
require(limma)
require(sva)
require(ggplot2)
require(ggrepel)
rm(list = ls())

args = commandArgs(T)
file_path = args[1] #data file path
design_path = args[2] #design matrix path
result_path = args[3] #result path
species = args[4] #decide species,such as human ,mouse
GPL_path=arg[5]#GPL information path
FC_setting = args[6] #setting fold change
pval_setting = args[7] #setting adj.p.val
FC_setting <- as.numeric(FC_setting)
pval_setting <- as.numeric(pval_setting)
logFC_setting <- log2(FC_setting)

dir.create(result_path)
# read design
design_mat = read.table(design_path, sep="\t", header=F)
gsm = as.vector(design_mat[,1])
grouplist = as.vector(design_mat[,2])
grouplist=gsub("-"," ",grouplist)#- is a sepcific signal in limma contrast_matrix
grouplist=gsub(" ",".",grouplist)

# read files in and process data
celFiles <- list.celfiles(file_path, full.names = T)
data.affy <- ReadAffy(filenames = celFiles)
data.rma <- rma(data.affy)
data.expr <- exprs(data.rma)

#Get gene symbol
gpl_inf=read.table(GPL_path,header=T,sep="\t",fill=T,quote = "")
symbol=gpl_inf$Gene.Symbol[match(rownames(data.expr),gpl_inf$ID)]
symbol[which(symbol=="")]="NULL"
symbol=sapply(strsplit(as.vector(unlist(symbol))," //"),function(x) x[1])
data.expr=aggregate(data.expr,by=list(symbol),mean)
rownames(data.expr)=data.expr[,1]
data.expr=data.expr[,-1]

#get exprSet  
#make sure the same sample order in design_matrix and expression profile
newcolname <- colnames(data.expr)[sapply(gsm,grep,colnames(data.expr),USE.NAMES = F)]
data.expr <- data.expr[,newcolname]
colnames(data.expr) <- gsm
write.table(data.expr,paste0(result_path,'/Expression_profile.txt'),sep = '\t',quote = F)
exprSet <- data.expr[,gsm]
# do pca analysis
pca_raw <- prcomp(t(data.expr), center = TRUE, scale. = F)
edata_pc_df <- as.data.frame(pca_raw$x)
edata_pc_df$class=grouplist
pca_plot=ggplot(edata_pc_df, aes(x = PC1, y = PC2, color = class)) +
    geom_point()
ggsave(paste0(result_path,"/PCA_plot.png"), pca_plot, width = 5, height = 5)

#do limma
design <- model.matrix(~0+factor(grouplist))
colnames(design) = levels(factor(grouplist))
rownames(design) = colnames(exprSet)
fit <- lmFit(exprSet,design)

condition_matrix=t(combn(unique(grouplist),2))
cond_matr=matrix((unlist(apply(condition_matrix,1,function(x) paste0(x,collapse = "-")))),nrow=1)
contrast.matrix <- makeContrasts(contrasts=cond_matr, levels =colnames(fit$coefficients))
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

each_condi_DEG<-function(x){
output <- topTable(fit,coef=x, number = Inf, lfc = 0, p.value = 1)
colnames(output) <- c("log2FoldChange","AveExpr","t","P.Value","padj","B")
# output file
write.table(output,paste0(result_path,'/',x,'Different_Expression.txt'),sep = '\t',quote = F)

sigout <- subset(output, abs(log2FoldChange) > logFC_setting & padj < pval_setting)
upR <- sigout[which(sigout$log2FoldChange > 0), 0]
downR <- sigout[which(sigout$log2FoldChange < 0), 0]

volca <- ggplot(output,aes(log2FoldChange,-log10(padj))) + 
  geom_point(aes(color=ifelse((log2FoldChange > logFC_setting | log2FoldChange < -logFC_setting) & padj < pval_setting, "sign","non-sign")), cex = 1, alpha = 0.3) + 
  scale_colour_manual("Significance", values = c("steelblue","red")) +
  labs(title = 'Volcano plot', x = "log2 Fold Change", y = "-log10 adj.P.Val")+
  geom_text_repel(data = sigout, aes(label = rownames(sigout)), cex = 3,box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))
ggsave(paste0(result_path,'/',x,"Volcano_plot.png"), volca, width = 7, height = 7)
}
apply(matrix(colnames(fit$coefficients)),1,each_condi_DEG)
