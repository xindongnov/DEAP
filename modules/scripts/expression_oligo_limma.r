required_Packages = c("limma","optparse")
if(!all(required_Packages %in% installed.packages())){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install(setdiff(required_Packages, installed.packages()))
}

require(limma)
require(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default=TRUE,dest='input',
              help="Input matrix path"),
  make_option(c("-l", "--gpl"), type="character", default=FALSE,
              help="Input GPL file path"),
  make_option(c("-r", "--result"), type="character", default=FALSE,
              help="Input differential expression result path"),
  make_option(c("-f", "--foldchange"), type="character", default='0',
              help="Set fold change value"),
  make_option(c("-q", "--fdr"), type="character", default="1",
              help="Set adj.p.val"),
  make_option(c("-t", "--treat"), type="character", default=TRUE,
              help="Input treatment"),
  make_option(c("-c", "--control"), type="character", default=TRUE,
              help="Input control"),
  make_option(c("--treatname"), type="character", default=TRUE,
              help="Input treatname sample"),
  make_option(c("--controlname"), type="character", default=TRUE,
              help="Input control sample")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input_path = opt$input
GPL_path = opt$gpl
result_path = opt$result
FC_setting = opt$foldchange
fdr_setting = opt$fdr
treatsample = opt$treat
controlsample = opt$control
treatname = opt$treatname
controlname = opt$controlname

treatsample = strsplit(treatsample,',')[[1]]
controlsample = strsplit(controlsample,',')[[1]]

FC_setting <- as.numeric(FC_setting)
fdr_setting <- as.numeric(fdr_setting)
logFC_setting <- log2(FC_setting)


input_matrix <- read.table(input_path, sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE,
                           quote = "", fill = TRUE)
dm <- as.data.frame(c(treatsample,controlsample))
colnames(dm) <- 'sample'
dm$treatment <- c(rep(controlname,nrow(dm)))
dm$treatment[match(treatsample, as.character(dm$sample))] <- treatname
dm$compare <- c(rep(1,nrow(dm)))
dm$compare[match(treatsample,dm$sample)] <- 2
dm <- dm[match(colnames(input_matrix),dm$sample),]

design <- model.matrix(~0+factor(dm[,3]))
rownames(design) <- dm$sample
colnames(design) <- c("Control", "Treatment")

fit <- lmFit(input_matrix,design)
contrast.mat <- makeContrasts(Treatment-Control, levels=design)
fit2 <- contrasts.fit(fit, contrast.mat)
fit2 <- eBayes(fit2)
limma_res <- topTable(fit2, number=nrow(input_matrix))

# Annotation
gpl_inf = read.table(GPL_path, sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE,
                     quote = "", fill = TRUE)

symbol=gpl_inf$gene_assignment[match(rownames(limma_res),gpl_inf$ID)]
symbol[which(symbol=="")]="NULL"
symbol=sapply(strsplit(as.vector(unlist(symbol))," //"),function(x) x[2])
limma_res_agg=aggregate(limma_res,by=list(symbol),max)
colnames(limma_res_agg)[1] <- 'symbol'

limma_res.new = limma_res_agg[order(limma_res_agg$adj.P.Val),]
limma_res.new1 = limma_res.new[!duplicated(limma_res.new$symbol),]
colnames(limma_res.new1) <- c("symbol","logFC","AveExpr","t","P.Value","adj.P.Val","B")

# Cut off the result
limma_res.new2 = limma_res.new1[(abs(limma_res.new1$logFC)>logFC_setting)&(limma_res.new1$adj.P.Val<fdr_setting),]
write.table(limma_res.new2, result_path, row.names = F, col.names = T, sep = '\t', quote = F)
