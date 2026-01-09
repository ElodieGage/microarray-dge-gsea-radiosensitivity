# ==========================================
# 03_GSEA_analysis.R
# Gene Set Enrichment Analysis using fGSEA
# ==========================================

# Libraries
library(dplyr)
library(fgsea)
library(biomaRt)

# Load DEG results
DEG_results <- read.csv(
  "results/tables/DEG_results.csv",
  row.names = 1
)


# ---- BioMart annotation ----

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)

# Fetch HGNC symbols for Affymetrix probe IDs
anno <- getBM(
  attributes = c("affy_hg_u133a_2", "hgnc_symbol"),
  filters = "affy_hg_u133a_2",
  values = rownames(DEG_results),
  mart = ensembl
)

# Remove empty and duplicated gene symbols
anno_clean <- anno %>%
  filter(hgnc_symbol != "") %>%
  distinct(hgnc_symbol, .keep_all = TRUE)

# Merge annotation with DEG table
DEG_anno <- DEG_results %>%
  tibble::rownames_to_column("probe_id") %>%
  inner_join(
    anno_clean,
    by = c("probe_id" = "affy_hg_u133a_2")
  )


# ---- Gene ranking for GSEA ----

gene_ranks <- DEG_anno$t
names(gene_ranks) <- DEG_anno$hgnc_symbol

# Sort ranks (required by fgsea)
gene_ranks <- sort(gene_ranks, decreasing = TRUE)


# ---- Load MSigDB pathways ----

kegg_pathways <- gmtPathways("data/c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt")
reactome_pathways <- gmtPathways("data/c2.cp.reactome.v2025.1.Hs.symbols.gmt")


# ---- Run fgsea ----

set.seed(2)

fgsea_kegg <- fgseaMultilevel(
  pathways = kegg_pathways,
  stats = gene_ranks,
  minSize = 15,
  maxSize = 500,
  eps = 0.0
)

fgsea_reactome <- fgseaMultilevel(
  pathways = reactome_pathways,
  stats = gene_ranks,
  minSize = 15,
  maxSize = 500,
  eps = 0.0
)

# Filter significant pathways 
fgsea_kegg_sig <- fgsea_kegg %>%
  filter(padj < 0.001 & abs(NES) >= 1.5)

fgsea_reactome_sig <- fgsea_reactome %>%
  filter(padj < 0.001 & abs(NES) >= 1.5)


# Save results
write.csv(
  fgsea_kegg,
  "results/tables/GSEA_KEGG_all.csv",
  row.names = FALSE
)

write.csv(
  fgsea_reactome,
  "results/tables/GSEA_REACTOME_all.csv",
  row.names = FALSE
)

write.csv(
  fgsea_kegg_sig,
  "results/tables/GSEA_KEGG_significant.csv",
  row.names = FALSE
)

write.csv(
  fgsea_reactome_sig,
  "results/tables/GSEA_REACTOME_significant.csv",
  row.names = FALSE
)
