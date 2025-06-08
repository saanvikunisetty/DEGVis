#install.packages("tidyverse")
#install.packages("pheatmap")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("airway", force = TRUE)

#library(tidyverse)
library(DESeq2)
library(airway)
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
#overall data cleaning

dds <- DESeqDataSet(airway_filtered, design = ~ cell + dex)
#initializing the dataset
dds <- DESeq(dds)
#run the pipeline

res <- results(dds, contrast = c("dex", "trt", "untrt"))
#extract results from dds object
#compare/contrast treated vs untreated (which genes are different between both?)
resOrdered <- res[order(res$padj), ]
print(head(resOrdered))
#sort from smallest to largest p-value
#determine most certain differentially expressed genes

plotMA(res, ylim = c(-5, 5))
#minus average plot
#log ratio against average intensity
#shows difference in gene expression between different samples

library(ggplot2)
library(ggrepel)
res_df <- as.data.frame(resOrdered)
#DEseq results converted to dataframe format
res_df <- na.omit(res_df)
#filter out NA values
res_df$significant <- res_df$padj < 0.05
res_df$gene <- rownames(res_df)
#add columns for significance and gene names
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("black", "red")) +
  geom_text_repel(
    data = subset(res_df, padj < 1e-100),  # Label only most significant ones
    aes(label = gene),
    size = 3,
    max.overlaps = 10
  ) +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10 p-value"
  ) +
  theme_minimal()
#volcano plot