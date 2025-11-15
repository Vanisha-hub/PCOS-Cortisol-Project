gc()
library(GEOquery)
library(limma)
library(dplyr)
library(biomaRt)
library(org.Hs.eg.db)
library(hgug4112a.db)
library(arrayQualityMetrics)
library(annotate)
library(ggplot2)
library(reshape2)

gse_data <- getGEO("GSE137684", GSEMatrix = TRUE)
untar("GSE137684/GSE137684_RAW.tar", exdir = "GSE137684_RAW")
files <- list.files("GSE137684_RAW", pattern = "txt.gz$", full.names = TRUE)
rawData <- read.maimages(files, source = "agilent", green.only = TRUE)
bg.corrected <- backgroundCorrect(rawData, method = "normexp", offset = 50)
normData <- normalizeBetweenArrays(bg.corrected, method = "quantile")
exprsData <- log2(normData$E)
exprsData <- exprs(gse_data[[1]])
rownames(exprs(gse_data[[1]]))[1:12]

gpl <- getGEO("GPL17077")
annot <- Table(gpl)[, c("ID", "GENE_SYMBOL")]
colnames(annot) <- c("ProbeID", "GeneSymbol")

expr_df <- data.frame(ProbeID = rownames(exprsData),
                      exprsData,
                      stringsAsFactors = FALSE)
expr_annot <- merge(expr_df, annot, by = "ProbeID", all.x = TRUE)

expr_annot <- expr_annot[!is.na(expr_annot$GeneSymbol) & expr_annot$GeneSymbol != "", ]

expr_gene <- expr_annot %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  as.data.frame()

rownames(expr_gene) <- expr_gene$GeneSymbol
expr_gene <- expr_gene[, -1]

melted <- melt(expr_gene, variable.name = "Sample", value.name = "Expression")
colnames(melted) <- c("Sample", "Expression")

#Boxplot Before and After
boxplot(rawData$E,
        main = "Boxplot: Before Normalization",
        las = 2,
        col = "blue",
        ylab = "Raw intensity (log scale)")


ggplot(melted, aes(x = Sample, y = Expression)) +
  geom_boxplot(fill = "skyblue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Boxplot of normalized log2 expression per sample")
colnames(expr_gene) <- gsub("GSE137684_RAW\\.|\\.txt", "", colnames(expr_gene))

#Density plot Before and After
matplot(density(rawData$E[,1])$x, 
        apply(rawData$E, 2, function(x) density(x)$y),
        type = "l", lty = 1, col = rainbow(ncol(rawData$E)),
        xlab = "Expression intensity", ylab = "Density",
        main = "Density plot: Before Normalization")
legend("topright", legend = colnames(rawData$E), col = rainbow(ncol(rawData$E)), 
       lty = 1, cex = 0.6)

ggplot(melted, aes(x = Expression, color = Sample)) +
  geom_density(alpha = 0.6) +
  theme_bw() +
  ggtitle("Density plot of normalized log2 expression") +
  theme(legend.position = "none")

write.table(expr_gene, "GSE137684_DEG_ready.txt",
            sep = "\t", quote = FALSE, row.names = TRUE)