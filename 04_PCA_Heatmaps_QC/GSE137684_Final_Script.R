##############################################################
# GSE137684 – Agilent single-color microarray
##############################################################

# Load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma","GEOquery","sva","pheatmap","ggplot2","dplyr","ggrepel"))
library(limma)
library(GEOquery)
library(sva)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Folders Created to store files
gse_id        <- "GSE137684"
meta_file     <- "metadata/GSE137684_series_matrix.txt"
plot_dir      <- file.path("QC", gse_id)
processed_dir <- file.path("processed_data", gse_id)

dir.create(plot_dir,      recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load data from series matrix
gse_list <- getGEO(filename = meta_file)
meta     <- pData(gse_list)
feature  <- fData(gse_list)
expr     <- exprs(gse_list)

# Group from metadata only
meta$Group <- ifelse(grepl("Normal|control", meta$characteristics_ch1.2, ignore.case = TRUE),
                     "Control", "PCOS")
meta$Group <- factor(meta$Group)

# Align everything perfectly
meta <- meta[colnames(expr), ]
rownames(meta) <- colnames(expr)

cat("Matched", ncol(expr), "samples\n")

# 2. Quantile normalization
expr_norm <- normalizeBetweenArrays(expr, method = "quantile")

# Save normalized data
write.csv(expr_norm, file.path(processed_dir, "GSE137684_NormalizedExpression.csv"), row.names = TRUE)

# 3. Robust gene symbol extraction
sym_col <- grep("symbol|gene.*name", colnames(feature), ignore.case = TRUE, value = TRUE)
if (length(sym_col) == 0) sym_col <- "ID"  # fallback

feature <- feature[rownames(expr_norm), ]  # align rows
gene_symbols <- feature[[sym_col[1]]]
gene_symbols <- sub(" *//.*$", "", gene_symbols)  # clean
gene_symbols[is.na(gene_symbols) | gene_symbols == "" | gene_symbols == "---"] <- rownames(expr_norm)

rownames(expr_norm) <- make.unique(gene_symbols)

percent_mapped <- round(100 * mean(!grepl("^A_", gene_symbols)), 1)  # rough estimate

# 4. QC plots
# Boxplot (labels fit)
png(file.path(plot_dir, "Boxplot.png"), width = 5000, height = 2400, res = 300)
par(mar = c(15,5,4,2))
boxplot(expr_norm, las = 2, col = ifelse(meta$Group == "Control", "#00BFC4", "#F8766D"),
        main = "Boxplot", ylab = "Expression", cex.axis = 0.9)
dev.off()

# Density plot (fixed name)
png(file.path(plot_dir, "DensityPlot.png"), width = 3200, height = 2400, res = 300)
group_colors <- ifelse(meta$Group == "Control", "#00BFC4", "#F8766D")
limma::plotDensities(expr_norm, col = group_colors, legend = FALSE)
legend("topright", legend = c("Control","PCOS"), col = c("#00BFC4","#F8766D"), lwd = 2)
dev.off()

# 5. Outlier detection
pca_temp <- prcomp(t(expr_norm), scale. = TRUE)
distances <- sqrt(rowSums(pca_temp$x[,1:2]^2))
outliers  <- names(which(distances > mean(distances) + 3*sd(distances)))
cat("Outliers flagged:", length(outliers), "\n")

# Keep all samples
data_clean <- expr_norm
meta_clean <- meta

# 6. Hierarchical clustering
if (ncol(data_clean) >= 2) {
  hc <- hclust(dist(t(data_clean)), method = "complete")
  png(file.path(plot_dir, "Clustering.png"), width = 3000, height = 2000, res = 300)
  plot(hc, labels = paste0(colnames(data_clean), " (", meta_clean$Group, ")"),
       main = "Hierarchical Clustering", cex = 0.8)
  dev.off()
}

# 7. Heatmap top 30
top_genes <- names(sort(apply(data_clean, 1, var, na.rm = TRUE), decreasing = TRUE))[1:30]
ann_col   <- data.frame(Group = meta_clean$Group, row.names = colnames(data_clean))

png(file.path(plot_dir, "Heatmap_Top30.png"), width = 4000, height = 2600, res = 300)
pheatmap(data_clean[top_genes, ],
         scale = "row",
         annotation_col = ann_col,
         color = colorRampPalette(c("blue","white","red"))(100),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_col = 8,
         main = "Top 30 Most Variable Genes - GSE137684")
dev.off()

# 8. PCA + Scree
pca_res <- prcomp(t(data_clean), scale. = TRUE)
pct     <- round(100*summary(pca_res)$importance[2,1:2], 1)

pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2],
                     Group = meta_clean$Group, Sample = colnames(data_clean))

png(file.path(plot_dir, "PCA.png"), width = 3400, height = 2400, res = 300)
ggplot(pca_df, aes(PC1, PC2, color = Group, label = Sample)) +
  geom_point(size = 5) + geom_text_repel(size = 3.5) +
  scale_color_manual(values = c("Control" = "#00BFC4", "PCOS" = "#F8766D")) +
  labs(title = "PCA Plot",
       x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)")) +
  theme_minimal(base_size = 14)
dev.off()

png(file.path(plot_dir, "Scree.png"), width = 3200, height = 2200, res = 300)
barplot(100*summary(pca_res)$importance[2,1:10], names.arg = paste0("PC",1:10),
        col = "steelblue", las = 2, ylab = "% Variance Explained",
        main = "Scree Plot")
dev.off()

# 9. Batch effect check + ComBat (safe)
batch_col <- grep("batch|run|plate|date", names(meta_clean), ignore.case = TRUE, value = TRUE)
batch <- if (length(batch_col)>0) factor(meta_clean[[batch_col[1]]]) else factor(rep("1", ncol(data_clean)))

# Only run ComBat if every batch has ≥2 samples
batch_ok <- all(table(batch) >= 2)
needs_combat <- nlevels(batch) > 1 && batch_ok

pca_before <- prcomp(t(data_clean), scale. = TRUE)
df_before  <- data.frame(PC1 = pca_before$x[,1], PC2 = pca_before$x[,2],
                         Group = meta_clean$Group, Batch = batch)
png(file.path(plot_dir, "Batch_Before.png"), width = 3000, height = 2200, res = 300)
ggplot(df_before, aes(PC1, PC2, color = Group, shape = Batch)) +
  geom_point(size = 5) + theme_minimal() +
  scale_color_manual(values = c("Control" = "#00BFC4", "PCOS" = "#F8766D")) +
  labs(title = "PCA Before Batch Correction")
dev.off()

if (needs_combat) {
  mod <- model.matrix(~ Group, data = meta_clean)
  data_combat <- ComBat(dat = data_clean, batch = batch, mod = mod)
  write.csv(data_combat, file.path(processed_dir, "GSE137684_ComBatCorrected.csv"))
  
  pca_after <- prcomp(t(data_combat), scale. = TRUE)
  df_after  <- data.frame(PC1 = pca_after$x[,1], PC2 = pca_after$x[,2],
                          Group = meta_clean$Group, Batch = batch)
  png(file.path(plot_dir, "Batch_After.png"), width = 3000, height = 2200, res = 300)
  ggplot(df_after, aes(PC1, PC2, color = Group, shape = Batch)) +
    geom_point(size = 5) + theme_minimal() +
    scale_color_manual(values = c("Control" = "#00BFC4", "PCOS" = "#F8766D")) +
    labs(title = "PCA After ComBat")
  dev.off()
  batch_note <- "Yes – correction applied"
} else {
  batch_note <- "No – not needed or not possible"
}

# 10. Final summary with all extra info
summary <- paste0(
  "GSE137684 finished!\n",
  "Samples: ", ncol(data_clean), "\n",
  "Probes mapped to symbols: ~", percent_mapped, "%\n",
  "Outliers flagged: ", length(outliers), "\n",
  "PC1 variance: ", pct[1], "%   PC2: ", pct[2], "%\n",
  "Batch levels: ", nlevels(batch), " → correction: ", batch_note, "\n",
  "All files saved in: ", plot_dir, " and ", processed_dir
)
writeLines(summary, file.path(plot_dir, "Conclusion.txt"))
cat(summary, "\n")
writeLines(capture.output(sessionInfo()), file.path(plot_dir, "sessionInfo.txt"))