# ===============================
# 01_limma_DEG_analysis.R
# Differential expression analysis using limma
# ===============================

# Libraries
library(limma)
library(dplyr)

# Log2 transform expression matrix
ExprMat_log <- log2(ExprMat + 1)

# Define experimental groups
pDataMat$group <- factor(pDataMat$title,
                         levels = c("Radioresistant", "Radiosensitive"))

# Design matrix
design <- model.matrix(~0 + group, data = pDataMat)
colnames(design) <- levels(pDataMat$group)

# Fit linear model
fit <- lmFit(ExprMat_log, design)

# Define contrast
contrast <- makeContrasts(Radiosensitive_vs_Radioresistant = Radiosensitive - Radioresistant,
                          levels = design)

# Apply contrast and empirical Bayes
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Extract all DEGs
DEG_results <- topTable(fit2,
                        number = Inf,
                        adjust.method = "BH")

# Save results
write.csv(DEG_results,
          file = "results/tables/DEG_results.csv",
          row.names = TRUE)

# Top 15 DEGs (by adjusted p-value)
top15_DEGs <- DEG_results %>%
  arrange(adj.P.Val) %>%
  head(15)

write.csv(top15_DEGs,
          file = "results/tables/DEG_top15.csv",
          row.names = TRUE)
