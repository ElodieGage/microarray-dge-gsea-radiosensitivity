# ==============================================
# 02_volcano_and_PCA.R
# Annotation and visualisation of DEG results
# ==============================================

# Libraries
library(ggplot2)
library(ggrepel)

# Load DEG results (from limma)
DEG_results <- read.csv(
  "results/tables/DEG_results.csv",
  row.names = 1
)

# Load: ExprMat_log, FeatureMat, pDataMat
load("data/processed/expression_objects.RData")


# ---- Annotation ----

# Merge limma results with feature annotation
res_annot <- merge(
  DEG_results,
  FeatureMat[, c("ID", "Gene Symbol")],
  by.x = "row.names",
  by.y = "ID",
  all.x = TRUE
)

# Clean up row names
rownames(res_annot) <- res_annot$Row.names
res_annot$Row.names <- NULL

# Remove probes without gene symbols
res_annot <- res_annot[!is.na(res_annot$`Gene Symbol`), ]


# ---- Volcano Plot ----

output <- res_annot
output$DiffExpressed <- "Not Sig"
output$DiffExpressed[output$logFC > 1 & output$adj.P.Val < 0.05] <- "Up"
output$DiffExpressed[output$logFC < -1 & output$adj.P.Val < 0.05] <- "Down"

# Select top 15 DEGs
top15table <- head(output[order(output$adj.P.Val), ], 15)

# Label top 15 DEGs
output$DElabel <- ifelse(
  rownames(output) %in% rownames(top15table),
  ifelse(!is.na(output$`Gene Symbol`),
         output$`Gene Symbol`,
         rownames(output)),
  ""
)

# Volcano plot
ggplot(
  data = subset(output, !is.na(adj.P.Val)),
  aes(x = logFC,
      y = -log10(adj.P.Val),
      col = DiffExpressed,
      label = DElabel)
) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  theme_minimal() +
  scale_color_manual(values = c("darkblue", "grey", "red")) +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    x = "Log2 Fold Change",
    y = "-Log10(adjusted p-value)",
    color = "Gene Regulation"
  )

ggsave(
  "results/figures/VolcanoPlot.png",
  width = 10,
  height = 7,
  dpi = 300
)


# ---- PCA ----

pca <- prcomp(t(ExprMat_log), scale. = FALSE)

pca_df <- data.frame(
  Sample = colnames(ExprMat_log),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = pDataMat$group
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 3) +
  theme_bw() +
  labs(
    title = "PCA Plot of Samples (All Genes)",
    x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)"),
    color = "Condition"
  )

ggsave(
  "results/figures/PCAplot.png",
  width = 10,
  height = 7,
  dpi = 300
)
