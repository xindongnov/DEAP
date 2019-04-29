required_Packages = c("oligo", "limma","sva")
if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}

require(oligo)
require(limma)
require(sva)
require(ggplot2)

args = commandArgs(T)
file_path = args[1] #data file path
design_path = args[2] #design matrix path
result_path = args[3] #result path
gpl_path = args[4]
FC_setting = args[5] #setting fold change
pval_setting = args[6] #setting adj.p.val
FC_setting <- as.numeric(FC_setting)
pval_setting <- as.numeric(pval_setting)
logFC_setting <- log2(FC_setting)
dir.create(result_path)

#design_path <- "/Users/shixiaoying/rongbin/gsm.txt"
#file_path <- "/Users/shixiaoying/rongbin/GSE41867_RAW"
#gpl_path  <- "/Users/shixiaoying/Downloads/GSE41867_family.soft"

ff <- gpl_path
nn <- grep("^[^#!^]", readLines(ff))[1] - 1
tmp <- system(paste0("grep -n '!Sample_organism_ch1' ",ff),intern = TRUE)[1]
if(gsub(".*= ","",tmp)=="Mus musculus"){
  db ="org.Mm.eg.db"
  if(db %in% installed.packages()){
    library(org.Mm.eg.db)
  }else{
    source("https://bioconductor.org/biocLite.R")
    biocLite("org.Mm.eg.db")
  }
}
if(gsub(".*= ","",tmp)=="Homo sapiens"){
  db="org.Hs.eg.db"
  if(db %in% installed.packages()){
    library(org.Hs.eg.db)
  }else{
    source("https://bioconductor.org/biocLite.R")
    biocLite("org.Hs.eg.db")
    library(org.Mm.eg.db)
  }
}

# Read design
design_mat = read.table(design_path, sep="\t", header=F)
gsm = as.vector(design_mat[,1])
grouplist = as.vector(design_mat[,2])

setwd(file_path)
# read files in and process data
celFiles <- list.celfiles(path = file_path)
data.raw <- read.celfiles(filenames = file.path('.', celFiles))
#pdf("Pre-normalization MA plot",width = 3.5,height = 4)
#MAplot(data=data.raw)
data.rma <- rma(data.raw)
data.expr <- exprs(data.rma)
time_batch <- substr(as.character(protocolData(data.raw)$dates),1,10)
dates <- as.character(protocolData(data.raw)$dates)
tmp <- sample(1:50, length(unique(time_batch)), replace=F)
data.batch <- time_batch
for(i in 1:length(unique(time_batch))){
  data.batch[which(time_batch==unique(time_batch)[i])] <- tmp[i]
}
# tmp <- sample(1:50, length(unique(grouplist)), replace=F)
# for(i in 1:length(unique(grouplist))){
#   grouplist[which(grouplist==unique(grouplist)[i])] <- tmp[i]
# }

sampleInfo <- data.frame(groupName=colnames(data.expr), date=dates, Condition=grouplist, batchIndex=time_batch)
{
  if(length(unique(data.batch)) == 1)
    data.combat <- data.expr else
      data.combat <- ComBat(dat = data.expr, batch = as.numeric(data.batch), mod = NULL)
}
colnames(data.expr) <- gsm

#Differential expression analysis
#convert probe ids to gene symbols
pfinfo <- read.table(ff, sep = "\t", quote = "", header = TRUE, skip = nn, fill = TRUE)
ID_idx <- which(colnames(pfinfo)=='ID')
gene_idx <- which(colnames(pfinfo)=='gene_assignment')
#gpl <- getGEO(getGEO(gsub("./|_.*",'',celFiles[1]))@header$platform_id)
#ID_idx <- which(colnames(Table(gpl))=='ID')
#gene_idx <- which(colnames(Table(gpl))=='gene_assignment')
#write.csv(Table(gpl)[,c(ID_idx,gene_idx)],"GPL16570.csv")
ids <- pfinfo[,c(ID_idx,gene_idx)]
ids_ok <- ids[grepl('^NM', ids$gene_assignment),] #select the lines with gene information
data.exprs_ann = merge(data.expr, ids_ok, by.x = 0, by.y = 'ID')
dim(data.exprs_ann)
data.exprs_ann$gene_assignment <- gsub(' /.*','',data.exprs_ann$gene_assignment)
data.exprs_ann_rd <- data.exprs_ann[!duplicated(data.exprs_ann$gene_assignment),]
dim(data.exprs_ann_rd)
rownames(data.exprs_ann_rd) <- data.exprs_ann_rd$gene_assignment
data.exprs_ann_rd <- data.exprs_ann_rd[,c(-1,-ncol(data.exprs_ann_rd))]
data.exprs_ann_rd1 <- data.exprs_ann_rd
{
  if(db=="org.Mm.eg.db")
    data.exprs_ann_rd1$SYMBOL <- select(org.Mm.eg.db, keys=rownames(data.exprs_ann_rd), columns="SYMBOL", keytype="REFSEQ")[,2]
  if(db=="org.Hs.eg.db")
    data.exprs_ann_rd1$SYMBOL <- select(org.Hs.eg.db, keys=rownames(data.exprs_ann_rd), columns="SYMBOL", keytype="REFSEQ")[,2]
}
data.exprs_ann_rd_sym <- na.omit(data.exprs_ann_rd1[!duplicated(data.exprs_ann_rd1$SYMBOL),])
rownames(data.exprs_ann_rd_sym) <- data.exprs_ann_rd_sym$SYMBOL
data.exprs_ann_rd_sym <- data.exprs_ann_rd_sym[,-ncol(data.exprs_ann_rd_sym)]
write.table(data.exprs_ann_rd_sym, paste0('.','/Full_data.txt'),sep = '\t',quote = F)


#do limma
design <- model.matrix(~0+factor(grouplist))
colnames(design) = levels(factor(grouplist))
rownames(design) = colnames(data.exprs_ann_rd_sym)
fit <- lmFit(data.exprs_ann_rd_sym,design)
contrast.matrix <- makeContrasts(paste0(rev(unique(grouplist)),collapse = '-'), levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
tempoutput <- topTable(fit, number = Inf, lfc = 0, p.value = 1)
colnames(tempoutput) <- c("log2FoldChange","AveExpr","t","P.Value","padj","B")
# output file
write.table(tempoutput,paste0('.','/Different_Expression.txt'),sep = '\t',quote = F)

sigout <- subset(tempoutput, abs(log2FoldChange) > logFC_setting & padj < pval_setting)
upR <- sigout[which(sigout$log2FoldChange > 0), 0]
downR <- sigout[which(sigout$log2FoldChange < 0), 0]
write.table(upR, paste0('.','/Upregulated_gene.txt'),sep = '\t',quote = F,col.names = F)
write.table(downR, paste0('.','/Downregulated_gene.txt'),sep = '\t',quote = F,col.names = F)

volca <- ggplot(tempoutput,aes(log2FoldChange,-log10(padj))) + 
  geom_point(aes(color=ifelse((log2FoldChange > logFC_setting | log2FoldChange < -logFC_setting) & padj < pval_setting, "sign","non-sign")), cex = 1, alpha = 0.3) + 
  scale_colour_manual("Significance", values = c("blue","red")) +
  labs(title = 'Volcano plot', x = "log2 Fold Change", y = "-log10 adj.P.Val")+
  geom_text(data = sigout, aes(x = log2FoldChange, y=-log10(padj), label = rownames(sigout)), cex = 1, vjust = "inward", hjust = "inward")
ggsave(paste0('.',"/Volcano_plot.png"), volca, width = 7, height = 7)
