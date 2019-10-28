# Microarray
# -*- coding: utf-8 -*-
# Created on Mon Oct 14 14:46:37 2019


required_Packages = c("affy","affyPLM")
if(!all(required_Packages %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(required_Packages, installed.packages()))
}

require(affyPLM)
require(affy)
require(optparse)

option_list <- list(
  make_option(c("-f", "--file"), type = "character", default=TRUE,
              help="Input data file path"),
  make_option(c("-l", "--gpl"), type="character", default=FALSE,
              help="Input GPL file path"),
  make_option(c("-p", "--probe"), type="character", default=FALSE,
              help="Input probe result file path"),
  make_option(c("-t", "--transcript"), type="character", default=FALSE,
              help="Input transcript result file path"),
  make_option(c("-g", "--gene"), type="character", default=FALSE,
              help="Input gene result file path"),
  make_option(c("-n", "--name"), type = "character", default=TRUE,
              help="Column name of output table header")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

names = opt$name
file_name = opt$file
GPL_path = opt$gpl
probe_path = opt$probe
transcript_path = opt$transcript
gene_path = opt$gene


data.affy <- ReadAffy(filenames = unlist(strsplit(file_name,',')))
eset = rma(data.affy)
data.expr = data.frame(exprs(eset))
data.expr = subset(data.expr, !grepl('AFFX-', rownames(data.expr))) # remove affy non-gene probes
colnames(data.expr) <- unlist(strsplit(names,','))
if(!is.na(probe_path)){
  write.table(data.expr,probe_path,col.names = T, sep = '\t',quote = F)
}

gpl_inf = read.table(GPL_path, sep = "\t", header = TRUE,
                 stringsAsFactors = FALSE,
                 quote = "", fill = TRUE)
#expression matrix: gene symbol
if(!is.na(gene_path)){
  symbol=gpl_inf$Gene.Symbol[match(rownames(data.expr),gpl_inf$ID)]
  symbol[which(symbol=="")]="NULL"
  symbol=sapply(strsplit(as.vector(unlist(symbol))," //"),function(x) x[1])
  data.expr1=aggregate(data.expr,by=list(symbol),max)
  rownames(data.expr1)=data.expr1[,1]
  data.expr1=data.expr1[,-1]
  write.table(data.expr1,gene_path,col.names = T, sep = '\t',quote = F)
}

#expression matrix: transcript symbol
if(!is.na(transcript_path)){
  transcript=gpl_inf$RefSeq.Transcript.ID[match(rownames(data.expr),gpl_inf$ID)]
  transcript[which(transcript=="")]="NULL"
  transcript=sapply(strsplit(as.vector(unlist(transcript))," //"),function(x) x[1])
  data.expr2=aggregate(data.expr,by=list(transcript),max)
  rownames(data.expr2)=data.expr2[,1]
  data.expr2=data.expr2[,-1]
  write.table(data.expr2,transcript_path, col.names = T, sep = '\t',quote = F)
}

