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
  make_option(c("-o", "--ont"),
    type = "character", default = "BP",
    help = "species, CC or BP or MF"
  ),
  make_option(c("--go_figure_path"),
    type = "character", default = "GO.png",
    help = "figure path"
  ),
  make_option(c("--go_table_path"),
    type = "character", default = "GO.csv",
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
ego_figure_path <- opt$go_figure_path
ego_table_path <- opt$go_table_path
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
  file.create(ego_figure_path)
  file.create(ego_table_path)
  print("no gene in gene list.")
  quit(save = "no", status = 0)
}

tryCatch(
  gene_entriz <- bitr(gene_list,
    fromType = "SYMBOL", toType = c("ENTREZID"),
    OrgDb = orgdb_string
  ),
  error = function(e) {
    file.create(ego_figure_path)
    file.create(ego_table_path)
    print("No enough genes transfered with ENTREZID.")
    quit(save = "no", status = 0)
  }
)

if (nrow(gene_entriz) == 0) {
  file.create(ego_figure_path)
  file.create(ego_table_path)
  print("No enough genes transfered with ENTREZID.")
  quit(save = "no", status = 0)
}


ego <- enrichGO(gene_entriz$ENTREZID,
  OrgDb = orgdb,
  ont = ont,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  keyType = "ENTREZID"
)

if (is.null(ego)) {
  file.create(ego_figure_path)
  file.create(ego_table_path)
  print("No enough genes for GO enrichment.")
  quit(save = "no", status = 0)
} else {
  p <- dotplot(ego, showCategory = 15, title = title) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = text_width))

  ggsave(ego_figure_path, p, width = 4, height = 5, scale = 1.5)

  ego_res_table <- as.data.frame(ego)
  write.table(ego_res_table, ego_table_path,
    quote = FALSE, row.names = FALSE, sep = "\t"
  )
}