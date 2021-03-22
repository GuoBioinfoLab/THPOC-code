
# Library -----------------------------------------------------------------

library(magrittr)
library(org.Hs.eg.db)

# Load data ---------------------------------------------------------------

panel <- readr::read_rds(file = "data/rda/panel.rds.gz")

# Functions ---------------------------------------------------------------


# Enrichment --------------------------------------------------------------

panel_entrez <- clusterProfiler::bitr(
  geneID = panel,
  fromType = "ENSEMBL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
) %>%
  dplyr::filter(! ENTREZID %in% c("109504726", "100529261"))

ggo <- clusterProfiler::groupGO(
  gene = panel_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  level = 3,
  readable = TRUE
)


kegg <- clusterProfiler::enrichKEGG(
  gene = panel_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

clusterProfiler::enrichMKEGG(gene = panel_entrez$ENTREZID, organism = "hsa")

# Save image --------------------------------------------------------------

save.image(file = "data/rda/04-enrich.rda")
