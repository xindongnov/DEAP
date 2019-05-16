
required_Packages = c("oligo","org.Mm.eg.db","org.Hs.eg.db")
if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}

require(oligo)
require(org.Mm.eg.db)
require(org.Hs.eg.db)
args = commandArgs(T)
file_path = args[1] #data file path
GPL_path=arg[2]#GPL information path
species=args[3]
design_path = args[4] #design matrix path
result_path = args[5] #result path

#file_path="/Users/yahan/Downloads/GSE41867_RAW"
#design_path="/Users/yahan/Downloads/GSE41867_RAW/GSE41867_designmatrix.txt"
#GPL_path="/Users/yahan/Downloads/GPL6246-18741.txt"
# Read design
design_mat = read.table(design_path, sep=",", header=T)
gsm = as.vector(design_mat[,1])

# read files in and process data
celFiles <- list.celfiles(file_path, listGzipped=TRUE,full.names=TRUE)
data.raw <- read.celfiles(celFiles)
data.rma <- rma(data.raw)
data.expr <- exprs(data.rma)

#Get gene symbol
gpl_inf=read.table(GPL_path,header=T,sep="\t",fill=T,quote = "")
symbol=gpl_inf$gene_assignment[match(rownames(data.expr),gpl_inf$ID)]
symbol=sapply(strsplit(as.vector(unlist(symbol))," //"),function(x) x[1])

if(species=="Mouse")
    new_symbol <- select(org.Mm.eg.db, keys=symbol, columns="SYMBOL", keytype="REFSEQ")[,2]
if(species=="Human")
    new_symbol <- select(org.Hs.eg.db, keys=symbol, columns="SYMBOL", keytype="REFSEQ")[,2]
new_symbol=toupper(new_symbol)
data.expr=aggregate(data.expr,by=list(new_symbol),median)

gsm = as.vector(design_mat$sample)
newcolname <-c("Group.1",colnames(data.expr)[sapply(gsm,grep,colnames(data.expr),USE.NAMES = F)])
data.expr <- data.expr[,newcolname]
colnames(data.expr) <- c("SYMBOL",gsm)
write.table(data.expr,result_path,sep = '\t',col.names = T,row.names = F,quote = F)



