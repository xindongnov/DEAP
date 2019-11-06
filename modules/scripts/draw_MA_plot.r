required_Packages = c("affy")
if(!all(required_Packages %in% installed.packages())){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install(setdiff(required_Packages, installed.packages()))
}

require(optparse)
require(ggplot2)
require(ggrepel)
require(affy)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default=FALSE,
              help="Input matrix path"),
  make_option(c("-r", "--result"), type="character", default=FALSE,
              help="Input plot result path")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input_path = opt$input
result_path = opt$result

# input_path  = "/Users/shixiaoying/Downloads/ll.txt" 
expr_matr=read.table(input_path, sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE,
                                quote = "", fill = TRUE)

pdf(result_path)
mva.pairs(data.frame(expr_matr),cex=0.5,main="MA plot")
dev.off()
