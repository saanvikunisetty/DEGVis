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

print(airway)
print(colnames(airway))
print(head(assay(airway)))
print(colData(airway))