
# ================================
# @auther: Ya Han, Xin Dong
# @email: xindong9511@gmail.com
# @date: Sep 2019
# ================================

required_Packages = c("affy", "sva")
if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}
require(affy)
require(sva)

args = commandArgs(T)
file_path = args[1] #data file path
GPL_path=args[2]#GPL information path
design_path = args[3] #design matrix path
result_path = args[4] #result path

#file_path <- "/Users/yahan/Downloads/GSE12056_RAW"
#GPL_path  <- "/Users/yahan/Downloads/GSE12056_RAW/GPL570-55999.txt"
#result_path="/Users/yahan/Downloads/GSE12056_RAW/exp.txt"
#design_path <- "/Users/yahan/Downloads/GSE12056_RAW/gse12056_design_matrix.txt"
design_mat = read.table(design_path, sep=",", header=T)

# read files in and process data
celFiles <- list.celfiles(file_path, full.names = T)
data.affy <- ReadAffy(filenames = celFiles)
data.rma <- rma(data.affy)
data.expr <- exprs(data.rma)

#Get gene symbol
gpl_inf=read.table(GPL_path,header=T,sep="\t",fill=T,quote = "")
symbol=gpl_inf$Gene.Symbol[match(rownames(data.expr),gpl_inf$ID)]
symbol[which(symbol=="")]=NA
symbol=sapply(strsplit(as.vector(unlist(symbol))," //"),function(x) x[1])
data.expr=aggregate(data.expr,by=list(symbol),median)

gsm = as.vector(design_mat$sample)
newcolname <- c("Group.1",colnames(data.expr)[sapply(gsm,grep,colnames(data.expr),USE.NAMES = F)])
data.expr <- data.expr[,newcolname]
colnames(data.expr) <- c("SYMBOL",gsm)
write.table(data.expr,result_path,sep = '\t',quote = F,row.names = F,col.names = T)