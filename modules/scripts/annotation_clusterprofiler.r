library(clusterProfiler, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
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
    type = "character", default = "hg",
    help = "species, hg or mm"
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
  make_option(c("-k", "--kegg"),
    type = "character", default = "y",
    help = "whether do KEGG, y or n"
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
dokegg <- opt$kegg
ego_figure_path <- opt$go_figure_path
ekegg_figure_path <- opt$kegg_figure_path
ego_table_path <- opt$go_table_path
ekegg_table_path <- opt$kegg_table_path
text_width <- opt$text_width

gene_list <- strsplit(gene_list, ",")[[1]]

if (species == "hs") {
  gene_entriz <- bitr(gene_list,
                      fromType = "SYMBOL", toType = c("ENTREZID"),
                      OrgDb = "org.Hs.eg.db")

  ego <- enrichGO(gene_entriz$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    keyType = "ENTREZID"
  )
  if (dokegg == "y") {
    ekegg <- enrichKEGG(
      gene = gene_entriz$ENTREZID,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      qvalueCutoff = 1
    )
  }
}

if (species == "mm") {
  gene_entriz <- bitr(gene_list,
                      fromType = "SYMBOL", toType = c("ENTREZID"),
                      OrgDb = "org.Mm.eg.db")

  ego <- enrichGO(gene_entriz$ENTREZID,
    OrgDb = org.Mm.eg.db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    keyType = "ENTREZID"
  )
  if (dokegg == "y") {
    ekegg <- enrichKEGG(
      gene = gene_entriz$ENTREZID,
      organism = "mmu",
      keyType = "kegg",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      qvalueCutoff = 1
    )
  }
}

p <- dotplot(ego, showCategory = 15, title = title) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = text_width))

ggsave(ego_figure_path, p, width = 4, height = 5, scale = 1.5)

ego_res_table <- as.data.frame(ego)
write.table(ego_res_table, ego_table_path,
            quote = FALSE, row.names = FALSE, sep = "\t")

if (dokegg == "y") {
  p <- dotplot(ekegg, showCategory = 15, title = title) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = text_width))
  ggsave(ekegg_figure_path, p, width = 4, height = 5, scale = 1.5)
  ekegg_res_table <- as.data.frame(ekegg)
  write.table(ekegg_res_table, ekegg_table_path,
              quote = FALSE, row.names = FALSE, sep = "\t")
}