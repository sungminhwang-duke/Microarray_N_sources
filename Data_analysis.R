### This contains codes for comprehensive analysis with microarray data

library(ggplot2)
library(dplyr)
library(tidyr)

# Input the normalized data
microarray_data <- read.csv("#Input .csv format")

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
H26_DE <- read.csv("#Input .csv format")
head(H26_DE)
with(H26_DE, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))

# Add color points as following categories: green if qvalue<0.05, orange if log2FC>1, red if both)
with(subset(H26_DE, qvalue<0.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(H26_DE, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(H26_DE, qvalue<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

# Identify the significant genes in the changed growth condition
print(subset(H26_DE, qvalue<0.05 & abs(log2FoldChange)>1))
arrangement <- print(subset(H26_DE, qvalue<0.05 & abs(log2FoldChange)>1))
print_arrangement <- print(arrangement[,2])
write.csv(arrangement, "#Assign where and what name of result .csv format")

# Label the significant genes from the calibrate plot
library(calibrate)

with(subset(H26_DE, qvalue<0.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.4))

# Cluster the significant genes and visualized by Heatmap
library("gplots")
library("heatmap.plus")
library("RColorBrewer")

test <- read.csv("#Input .csv format", row.names = 1)
NH4Cl <- test[grep("NH4",row.names(test)),]
Alanine <- test[grep("alanine",row.names(test)),]

# http://www.rapidtables.com/web/color/RGB_Color.htm
condition_colors <- unlist(lapply(rownames(test),function(x){
  if(grepl("alanine",x)) '#FFC0CB' #pink
  else if(grepl('NH4',x)) '#808080' #grey
}))

input <- as.matrix(t(test))
heatmap.2(input, hclust=function(x) hclust(x,method="average"), trace="none", density="none", col=bluered(20), cexRow=0.5, cexCol=0.3, margins = c(5,10),
          ColSideColors=condition_colors, scale="row", main="Less than 5% FDR")
legend(0.75, 0.9,legend=c("NH4Cl","Alanine"),fill=c('#808080','#FFC0CB'),cex=0.3)

