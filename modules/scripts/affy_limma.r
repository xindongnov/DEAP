required_Packages = c("affy", "limma","sva","org.Hs.eg.db","ggplot2")

if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}

require(affy)
require(limma)
# require(hgu133plus2hsrefseqcdf)
require(org.Hs.eg.db)
require(sva)
require(ggplot2)

rm(list = ls())

# need GSEnum
# GSEnum <- 'GSE12056'
# designname <- 'A2_GSE12056_CREB1_KO.txt'
# #need design_path
# design_path <- paste0('microarray/design_matrix/',designname)
# #need data_path
# file_path <- paste0('microarray/data/',GSEnum,'/')
# #need result_path
# result_path <- paste0('microarray/result/GSE12056/')
# logFC_setting = log2(1) #setting log2 fold change
# pval_setting = 0.1 #setting adj.p.val
# cdf_name = 'hgu133plus2hsrefseqcdf'

args = commandArgs(T)
file_path = args[1] #data file path
design_path = args[2] #design matrix path
result_path = args[3] #result path
cdf_name = args[4] #cdf names
FC_setting = args[5] #setting fold change
pval_setting = args[6] #setting adj.p.val
FC_setting <- as.numeric(FC_setting)
pval_setting <- as.numeric(pval_setting)
logFC_setting <- log2(FC_setting)
require(cdf_name, character.only = T)

dir.create(result_path)

# read design
design_mat = read.table(design_path, sep="\t", header=F)
gsm = as.vector(design_mat[,1])
grouplist = as.vector(design_mat[,2])

# read files in and process data
celFiles <- list.celfiles(path = file_path, full.names = T)
data.affy <- ReadAffy(filenames = celFiles)
data.affy@cdfName = cdf_name
data.rma <- rma(data.affy)
data.expr <- exprs(data.rma)
# batcheffect
data.affy@protocolData@data
# groupName <- colnames(data.expr)
# batchIndex <- c(1,2,2,3,3,1)
# Condition <- c(1,1,1,2,2,2)
# batchInfo <- data.frame(groupName, batchIndex, Condition)
# cbmod <- model.matrix(~Condition, data = batchInfo)
# data.expr <- ComBat(dat = data.expr, batch = batchIndex, mean.only = TRUE, mod = cbmod)

#get exprSet
newcolname <- colnames(data.expr)[sapply(gsm,grep,colnames(data.expr),USE.NAMES = F)]
data.expr <- data.expr[,newcolname]
colnames(data.expr) <- gsm
#colnames(data.expr) <- sapply(sapply(colnames(data.expr),strsplit,'_',USE.NAMES = F),head,1,USE.NAMES = F)
exprSet <- data.expr[,gsm]

#do limma
design <- model.matrix(~0+factor(grouplist))
colnames(design) = levels(factor(grouplist))
rownames(design) = colnames(exprSet)
fit <- lmFit(exprSet,design)
contrast.matrix <- makeContrasts(paste0(rev(unique(grouplist)),collapse = '-'), levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
tempoutput <- topTable(fit, number = Inf, lfc = 0, p.value = 1)

#Get gene symbol
rownames(tempoutput) <- sapply(sapply(rownames(tempoutput),strsplit,'[.]',USE.NAMES = F),head,1,USE.NAMES = F)
tempoutput$SYMBOL <- select(org.Hs.eg.db, keys=rownames(tempoutput), columns="SYMBOL", keytype="REFSEQ")[,2]
output <- na.omit(tempoutput[!duplicated(tempoutput$SYMBOL),])
rownames(output) <- output$SYMBOL
output <- output[,-7]
colnames(output) <- c("log2FoldChange","AveExpr","t","P.Value","padj","B")
# output file
write.table(output,paste0(result_path,'/Different_Expression.txt'),sep = '\t',quote = F)

sigout <- subset(output, abs(log2FoldChange) > logFC_setting & padj < pval_setting)
upR <- sigout[which(sigout$log2FoldChange > 0), 0]
downR <- sigout[which(sigout$log2FoldChange < 0), 0]
write.table(upR, paste0(result_path,'/Upregulated_gene.txt'),sep = '\t',quote = F,col.names = F)
write.table(downR, paste0(result_path,'/Downregulated_gene.txt'),sep = '\t',quote = F,col.names = F)

volca <- ggplot(output,aes(log2FoldChange,-log10(padj))) + 
  geom_point(aes(color=ifelse((log2FoldChange > logFC_setting | log2FoldChange < -logFC_setting) & padj < pval_setting, "sign","non-sign")), cex = 1, alpha = 0.3) + 
  scale_colour_manual("Significance", values = c("steelblue","red")) +
  labs(title = 'Volcano plot', x = "log2 Fold Change", y = "-log10 adj.P.Val")+
  geom_text(data = sigout, aes(x = log2FoldChange, y=-log10(padj), label = rownames(sigout)), cex = 1, vjust = "inward", hjust = "inward")
ggsave(paste0(result_path,"/Volcano_plot.png"), volca, width = 7, height = 7)

