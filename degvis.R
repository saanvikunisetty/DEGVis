#install.packages("tidyverse")
#install.packages("pheatmap")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("DESeq2")
#BiocManager::install("airway", force = TRUE)

#library(tidyverse)
#library(DESeq2)
#library(airway)
#library(pheatmap)

data("airway")
library(SummarizedExperiment)
#airway dataset
#how human airway smooth muscle responds to DEX (asthma medication)

print(airway)
print(colnames(airway))
print(head(assay(airway)))
print(colData(airway))
#initial data inspection

keep_genes <- rowSums(assay(airway)) >= 10
airway_filtered <- airway[keep_genes, ]
#filter out any genes with really low total sample counts
colData(airway_filtered)$dex <- as.factor(colData(airway_filtered)$dex)
colData(airway_filtered)$cell <- as.factor(colData(airway_filtered)$cell)
#treat different sample types as categories
dim(airway)         
dim(airway_filtered)