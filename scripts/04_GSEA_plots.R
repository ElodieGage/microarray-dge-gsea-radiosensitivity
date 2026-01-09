# ===================================
# 04_GSEA_plots.R
# Visualisation of GSEA results
# ===================================

# Libraries
library(ggplot2)
library(dplyr)

# ---- Load GSEA results ----
kegg_sig <- read.csv(
  "results/tables/GSEA_KEGG_significant.csv"
)

reactome_sig <- read.csv(
  "results/tables/GSEA_REACTOME_significant.csv"
)

# Add pathway database labels
kegg_sig$pathwayDatabase <- "KEGG"
reactome_sig$pathwayDatabase <- "REACTOME"

gsea_combined <- bind_rows(kegg_sig, reactome_sig)

# Define direction of enrichment
gsea_combined <- gsea_combined %>%
  mutate(
    Enriched_in = ifelse(NES > 0, "Radiosensitive", "Radioresistant")
  )


# ---- Assign biological process categories ----

gsea_combined <- gsea_combined %>%
  mutate(
    biologicalProcess = case_when(
      
      # Immune system
      grepl("INTERFERON|CYTOKINE|ANTIGEN|ANTIVIRAL|CHEMOKINE", pathway, ignore.case = TRUE) ~
        "Immune System",

      # DNA-related processes
      grepl("REPLICATION|REPAIR|TRANSLESION", pathway, ignore.case = TRUE) ~
        "DNA Maintenance",

      # Cell cycle
      grepl("CYCLIN|MITOTIC|CELL_CYCLE", pathway, ignore.case = TRUE) ~
        "Cell Cycle",

      # Chromatin and epigenetics
      grepl("CHROMATIN|METHYLATION|HISTONE|PRC2|HDAC|SIRT", pathway, ignore.case = TRUE) ~
        "Chromatin Organisation",

      # Protein & RNA metabolism
      grepl("TRANSLATION|RIBOSOME|MRNA|RRNA", pathway, ignore.case = TRUE) ~
        "Gene Expression",

      # Developmental / keratinisation
      grepl("KERATIN|CORNIFIED|DIFFERENTIATION", pathway, ignore.case = TRUE) ~
        "Developmental Biology",

      TRUE ~ "Other"
      ))


# ---- Order biological processes ----

gsea_combined$biologicalProcess <- factor(
  gsea_combined$biologicalProcess,
  levels = c(
    "Immune System",
    "DNA Maintenance",
    "Cell Cycle",
    "Chromatin Organisation",
    "Gene Expression",
    "Developmental Biology",
    "Other"
  ))


# ---- Plot enrichment by condition ----

p_condition <- ggplot(
  gsea_combined,
  aes(x = NES, y = biologicalProcess, colour = Enriched_in)
) +
  geom_jitter(height = 0.2, width = 0, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  theme_bw() +
  labs(
    title = "Differential pathway enrichment in radiosensitivity",
    subtitle = "Pathways grouped by biological process",
    x = "Normalized Enrichment Score (NES)",
    y = "Biological Process",
    colour = "Enriched in"
  ) +
  scale_colour_manual(
    values = c(
      "Radiosensitive" = "red",
      "Radioresistant" = "darkblue"
    ))

ggsave(
  "results/figures/GSEA_enrichment_by_condition.png",
  p_condition,
  width = 8,
  height = 5,
  dpi = 300
)

# ---- Plot enrichment by database ----
p_database <- ggplot(
  gsea_combined,
  aes(x = NES, y = biologicalProcess, colour = pathwayDatabase)
) +
  geom_jitter(height = 0.2, width = 0, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  theme_bw() +
  labs(
    title = "Pathway enrichment by database",
    subtitle = "KEGG and Reactome pathways",
    x = "Normalized Enrichment Score (NES)",
    y = "Biological Process",
    colour = "Database"
  ) +
  scale_colour_manual(
    values = c(
      "KEGG" = "darkorange",
      "REACTOME" = "purple4"
    ))

ggsave(
  "results/figures/GSEA_enrichment_by_database.png",
  p_database,
  width = 8,
  height = 5,
  dpi = 300
)
