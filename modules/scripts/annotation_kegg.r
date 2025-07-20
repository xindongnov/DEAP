library(clusterProfiler, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(optparse)
library(stringr)


option_list <- list(
  make_option(c("-t", "--title"),
    type = "character", default = FALSE,
    help = "Figure title"
  ),
  make_option(c("-l", "--list"),
    type = "character", default = FALSE,
    help = "gene list"
  ),
  make_option(c("-s", "--species"),
    type = "character", default = "hs",
    help = "species, hs, mm or rn"
  ),
  make_option(c("--kegg_figure_path"),
    type = "character", default = "kegg.png",
    help = "figure path"
  ),
  make_option(c("--kegg_table_path"),
    type = "character", default = "kegg.csv",
    help = "table path"
  ),
  make_option(c("-w", "--text_width"),
    type = "integer", default = 30,
    help = "table path"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

title <- opt$title
gene_list <- opt$list
species <- opt$species
ont <- opt$ont
ekegg_figure_path <- opt$kegg_figure_path
ekegg_table_path <- opt$kegg_table_path
text_width <- opt$text_width

gene_list <- read.csv(gene_list, header = FALSE)$V1

if (species == "hs") {
  orgdb_string <- "org.Hs.eg.db"
  orgdb <- org.Hs.eg.db
  organism <- "hsa"
} else if (species == "mm") {
  orgdb_string <- "org.Mm.eg.db"
  orgdb <- org.Mm.eg.db
  organism <- "mmu"
} else if (species == "rn") {
  orgdb_string <- "org.Rn.eg.db"
  orgdb <- org.Rn.eg.db
  organism <- "rno"
} else {
  stop("Unsupported species. Please use 'hs', 'mm', or 'rn'.")
}


if (length(gene_list) == 0) {
  file.create(ekegg_figure_path)
  file.create(ekegg_table_path)
  print("no gene in gene list.")
  quit(save = "no", status = 0)
}

tryCatch(
  gene_entriz <- bitr(gene_list,
    fromType = "SYMBOL", toType = c("ENTREZID"),
    OrgDb = orgdb_string
  ),
  error = function(e) {
    file.create(ekegg_figure_path)
    file.create(ekegg_table_path)
    print("No enough genes transfered with ENTREZID.")
    quit(save = "no", status = 0)
  }
)

if (nrow(gene_entriz) == 0) {
  file.create(ekegg_figure_path)
  file.create(ekegg_table_path)
  print("No enough genes transfered with ENTREZID.")
  quit(save = "no", status = 0)
}

ekegg <- enrichKEGG(
  gene = gene_entriz$ENTREZID,
  organism = organism,
  keyType = "kegg",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1
)

if (is.null(ekegg) || nrow(ekegg) == 0) {
  file.create(ekegg_figure_path)
  file.create(ekegg_table_path)
  print("No KEGG enrichment results found.")
  quit(save = "no", status = 0)
} else {
  p <- dotplot(ekegg, showCategory = 15, title = title) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = text_width))
  ggsave(ekegg_figure_path, p, width = 4, height = 5, scale = 1.5)

  ekegg_res_table <- as.data.frame(ekegg)
  write.table(ekegg_res_table, ekegg_table_path,
    quote = FALSE, row.names = FALSE, sep = "\t"
  )
}