### This contains codes for comprehensive analysis with microarray data

library(ggplot2)
library(dplyr)
library(tidyr)

# Input the normalized data
microarray_data <- read.csv(file.choose("H26_all_genes_normalized.csv"))

# Identify the correlation coefficient and the significance level(or p-value) from the two different N-sources (alanine vs ammonium chloride) with each biological replicates

library(ggpubr)
require(gridExtra)

plot1 <- ggscatter(microarray_data, x = "NH4_1", y = "NH4_2", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "NH4_1", ylab = "NH4_2")

plot2 <- ggscatter(microarray_data, x = "NH4_2", y = "NH4_3", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "NH4_2", ylab = "NH4_3")

plot3 <- ggscatter(microarray_data, x = "NH4_1", y = "NH4_3", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "NH4_1", ylab = "NH4_3")

grid.arrange(plot1, plot2, plot3, ncol=3)


plot4 <- ggscatter(microarray_data, x = "alanine_1", y = "alanine_2", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "alanine_1", ylab = "alanine_2")

plot5 <- ggscatter(microarray_data, x = "alanine_2", y = "alanine_3", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "alanine_2", ylab = "alanine_3")

plot6 <- ggscatter(microarray_data, x = "alanine_1", y = "alanine_3", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "alanine_1", ylab = "alanine_3")

grid.arrange(plot4, plot5, plot6, ncol=3)


# Plot a volcano plot
H26_DE <- read.csv(file.choose("H26_all_genes_normalized.csv"))
head(H26_DE)
with(H26_DE, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))

# Add color points as following categories: red if both qvalue<0.05 and log2FC>1)
with(subset(H26_DE, qvalue<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

# Identify the significant genes in the changed growth condition
print(subset(H26_DE, qvalue<0.05 & abs(log2FoldChange)>1))
arrangement <- print(subset(H26_DE, qvalue<0.05 & abs(log2FoldChange)>1))
print_arrangement <- print(arrangement[,2])
write.csv(arrangement, "Sig_genes.csv")

# Label the significant genes from the calibrate plot, if needed.
library(calibrate)

with(subset(H26_DE, qvalue<0.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.4))

# Cluster the significant genes and visualized by Heatmap (hierarchical clusterings with dendextend)
library(tidyr)
library(dplyr)
library(dendextend)
library(RColorBrewer)
library(gplots)

sig.gene <- read.csv(file.choose("Significant_genes.csv"))
sig.gene.long <- gather(sig.gene, gene, expression, -condition)

sig.gene.cor <- select(sig.gene, -condition) %>% 
  cor(use="pairwise.complete.obs", method = "pearson")
sig.gene.dist <- as.dist(1 - sig.gene.cor)

sig.gene.tree <- hclust(sig.gene.dist, method="complete")
plot(sig.gene.tree)
sig.gene.dend <- as.dendrogram(sig.gene.tree)

plot(sig.gene.dend, leaflab = "none")
sig.clusters <- cutree(sig.gene.dend, k=3)

data.mtx <- as.matrix(select(sig.gene, -condition))
row.names(data.mtx) <- sig.gene$condition
transposed.data.mtx <- t(data.mtx)

color.scheme <- rev(brewer.pal(8,"RdYlBu"))

heatmap.2(transposed.data.mtx,
          Rowv = sig.gene.dend,  
          Colv = NULL,
          dendrogram = "row", 
          col = color.scheme,
          trace = "none", density.info = "none",
          scale="row",
          cexRow=0.7, cexCol=0.5, 
          margins = c(5,5),
          xlab = "Growth condition")
