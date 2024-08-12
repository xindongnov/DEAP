library(clusterProfiler, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(optparse)
library(stringr)

option_list <- list(
  make_option(c("-t", "--title"), type = "character", default = FALSE,
              help = "Figure title"),
  make_option(c("-l", "--list"), type = "character", default = FALSE,
              help = "gene list"),
  make_option(c("-s", "--species"), type = "character", default = 'hg',
              help = "species, hg or mm"),
  make_option(c("-o", "--ont"), type = "character", default = 'BP',
              help = "species, CC or BP or MF"),
  make_option(c("-f", "--figure_path"), type = "character", default = 'GO.png',
              help = "figure path"),
  make_option(c("-p", "--table_path"), type = "character", default = 'GO.csv',
              help = "table path")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

title = opt$title
gene_list = opt$list
species = opt$species
ont = opt$ont
figure_path = opt$figure_path
table_path = opt$table_path

gene_list <- strsplit(gene_list, ',')[[1]]

if (species == 'hs') {
  gene_entriz <- bitr(gene_list, fromType = "SYMBOL", toType = c("ENTREZID"),
                      OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene_entriz$ENTREZID, OrgDb = org.Hs.eg.db, ont = ont,
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 1,
                 qvalueCutoff = 1,
                 keyType = 'ENTREZID')
}

if (species == 'mm'){
  gene_entriz <- bitr(gene_list, fromType = "SYMBOL", toType = c("ENTREZID"),
                      OrgDb = "org.Mm.eg.db")
  go <- enrichGO(gene_entriz$ENTREZID, OrgDb = org.Mm.eg.db, ont = ont,
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 1,
                 qvalueCutoff = 1,
                 keyType = 'ENTREZID')
}

p <- dotplot(go, showCategory = 10, title = title)

ggsave(figure_path, p, width = 5, height = 5, scale = 1.5) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))

res_table <- as.data.frame(go)
write.csv(res_table, table_path, quote = FALSE, row.names = FALSE)
