required_Packages = c("oligo", "limma","sva")

if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}

library(oligo)
library(limma)
#library(hgu133plus2.db)
#library(sva)

#need GSEnum
GSEnum <- 'GSE38956'
designname <- 'A1_GSE38956_SIN3A_KD.txt'
#need design_path
design_path <- paste0('microarray/design_matrix/',designname)
#need data_path
file_path <- paste0('microarray/data/',GSEnum,'/')
#need result_path
result_path <- paste0('microarray/result/',designname)

# Read design
design_mat = read.table(design_path, sep="\t", header=F)
gsm = as.vector(design_mat[,1])
grouplist = as.vector(design_mat[,2])

# read file in and process data
celFiles <- list.celfiles(path = file_path, full.names = TRUE)
rawData <- read.celfiles(celFiles)
data.rma <- rma(rawData)
data.expr <- exprs(data.rma)

# batcheffect
rawData@protocolData@data
# groupName <- colnames(data.expr)
# batchIndex <- c()
# Condition <- c()
# batchInfo <- data.frame(groupName, batchIndex, Condition)
# cbmod <- model.matrix(~Condition, data = batchInfo)
# data.expr <- ComBat(dat = data.expr, batch = batchIndex, mean.only = TRUE, mod = cbmod)

#Get gene symbol
#Annot <- data.frame(SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "))
#cpltedata <- merge(Annot, data.expr, by.x=0, by.y=0, all = TRUE)

#get exprSet
dataindex <- sapply(gsm,paste0,'.CEL',USE.NAMES = F)
exprSet <- data.expr[,dataindex]

# do limma
design <- model.matrix(~0+factor(grouplist))
colnames(design) = levels(factor(grouplist))
rownames(design) = colnames(exprSet)
fit <- lmFit(exprSet,design)
contrast.matrix <- makeContrasts(paste0(unique(grouplist),collapse = '-'), levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
output <- topTable(fit, number = Inf, lfc = log2(1.5), p.value = 0.05)
write.table(output,result_path,sep = '\t',quote = F)


