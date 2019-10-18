
# ================================
# @auther: Xin Dong
# @email: xindong9511@gmail.com
# @date: Sep 2019
# ================================

required_Packages = c("lumi", "limma")

if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}

require(lumi)
require(limma)

#need GSEnum
GSEnum <- 'GSE30622'
#need design_path
design_path <- ''
file_path <- ''

fileName <- file_path
x.lumi <- lumiR.batch(fileName)
lumi.N.Q <- lumiExpresso(x.lumi[[1]])

design_mat = read.table(design_path, sep="\t", header=F)
gsm = as.vector(design_mat[,1])
grouplist = as.vector(design_mat[,2])

exprSet <- exprs(lumi.N.Q)[,gsm]
ID <- IlluminaID2nuID(rownames(exprSet),species = "Human")
rownames(exprSet) <- ID[,4]

design <- model.matrix(~0+factor(grouplist))
colnames(design) = levels(factor(grouplist))
rownames(design) = colnames(exprSet)
fit <- lmFit(exprSet,design)
contrast.matrix <- makeContrasts(paste0(unique(grouplist),collapse = '-'), levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
output <- topTable(fit, number = Inf, lfc = log2(1.5), p.value = 0.05)
write.table(output,'Result_limma.txt',sep = '\t',quote = F)