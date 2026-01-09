# **Differential Expression and GSEA Analysis of Radiosensitivity**

## **Project Overview**

This project investigates genomic differences associated with cellular radiosensitivity using microarray expression data. Differential gene expression analysis was performed using the limma package, followed by pathway-level interpretation using Gene Set Enrichment Analysis (GSEA).

The aim of the analysis was to identify biological processes and signalling pathways that separate radiosensitive from radioresistant samples.


## **Analysis Pipeline**

The analysis was carried out in R. The main steps are described below:


## **Repository Structure**
```
microarray-dge-gsea-radiosensitivity/
  data/ 
  scripts/
    01_preprocessing.R
    02_limma_DE.R
    03_GSEA_analysis.R
    04_GSEA_plots.R
  results/
    figures/ 
    tables/ 
  README.md
```

## **Key Results**


## **Reproducibility**

The analysis was performed in R using the following key packages:

- limma
- fgsea
- biomaRt
- ggplot2
- dplyr

Scripts are numbered to indicate the recommended order of execution.


## Attribution

This project was completed as part of a bioinformatics course. 
