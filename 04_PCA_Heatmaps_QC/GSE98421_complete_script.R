##############################################################
# GSE114419 – Affymetrix Human Transcriptome Array 2.0 (GPL17586)
# Simple step-by-step analysis
##############################################################

# 1. Install and load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery","oligo","arrayQualityMetrics","pd.hta.2.0","hta20transcriptcluster.db",
  "AnnotationDbi","limma","sva","pheatmap","ggplot2","dplyr","ggrepel"
))

library(GEOquery)
library(oligo)
library(arrayQualityMetrics)
library(AnnotationDbi)
library(hta20transcriptcluster.db)
library(limma)
library(sva)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)

# 2. Set folder paths
gse_id        <- "GSE114419"
raw_dir       <- file.path("Raw_data", gse_id)
meta_file     <- file.path("metadata", paste0(gse_id, "_series_matrix.txt"))
plot_dir      <- file.path("QC", gse_id)
processed_dir <- file.path("processed_data", gse_id)

dir.create(plot_dir,      recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Load sample information (metadata)
gse_list <- getGEO(filename = meta_file, getGPL = FALSE)
gse      <- gse_list
meta     <- pData(gse)

# Make Group column: PCOS or Control
meta$Group <- ifelse(grepl("PCOS", meta$title, ignore.case = TRUE) |
                       grepl("PCOS", meta$source_name_ch1, ignore.case = TRUE),
                     "PCOS", "Control")
meta$Group <- factor(meta$Group)

# 4. Read the raw .CEL files
cel_files <- list.celfiles(raw_dir, full.names = TRUE)
raw_data  <- read.celfiles(cel_files)

# 5. Quality check before normalization
arrayQualityMetrics(raw_data,
                    outdir = file.path(plot_dir, "QC_Raw"),
                    force = TRUE)

# 6. Normalize with RMA (core level for HTA 2.0)
norm_data <- rma(raw_data, target = "core")
expr      <- exprs(norm_data)

# Save normalized expression data
write.csv(expr, file.path(processed_dir, "GSE114419_NormalizedExpression.csv"), row.names = TRUE)

# 7. Quality check after normalization
arrayQualityMetrics(norm_data,
                    outdir = file.path(plot_dir, "QC_Normalized"),
                    force = TRUE)

# 8. Remove low-intensity probes
row_median <- rowMedians(expr)
filtered   <- expr[row_median > 3.5, ]

# 9. Convert probe IDs to gene symbols
probe_ids    <- rownames(filtered)
gene_symbols <- mapIds(hta20transcriptcluster.db,
                       keys = probe_ids,
                       keytype = "PROBEID",
                       column = "SYMBOL",
                       multiVals = "first")

# Keep only probes with a gene symbol
keep      <- !is.na(gene_symbols)
annotated <- filtered[keep, ]
rownames(annotated) <- gene_symbols[keep]

# Average duplicate genes
data <- limma::avereps(annotated)

# 10. Find outliers using PCA
pca <- prcomp(t(data), scale. = TRUE)
distances <- sqrt(rowSums(pca$x[,1:2]^2))
outliers  <- names(which(distances > mean(distances) + 3*sd(distances)))
message("Possible outliers: ", ifelse(length(outliers)>0, paste(outliers, collapse=", "), "None"))

# 11. Density plot
png(file.path(plot_dir, "DensityPlot.png"), width = 3200, height = 2400, res = 300)
group_colors <- ifelse(meta$Group == "Control", "#00BFC4", "#F8766D")
limma::plotDensities(data, col = group_colors, legend = FALSE,
                     main = "Density Plot")
legend("topright", legend = c("Control","PCOS"), col = c("#00BFC4","#F8766D"), lwd = 2)
dev.off()

# 12. Boxplot 
png(file.path(plot_dir, "Boxplot.png"), width = 5000, height = 2400, res = 300)
par(mar = c(15,5,4,2))
boxplot(data, las = 2, col = ifelse(meta$Group == "Control", "#00BFC4", "#F8766D"),
        main = "Boxplot", ylab = "Expression", cex.axis = 0.9)
dev.off()

# 13. Hierarchical clustering
if (ncol(data) >= 2) {
  hc <- hclust(dist(t(data)), method = "complete")
  png(file.path(plot_dir, "Clustering.png"), width = 3000, height = 2000, res = 300)
  plot(hc, labels = paste0(colnames(data), " (", meta$Group, ")"),
       main = "Hierarchical Clustering", cex = 0.8)
  dev.off()
}

# 14. Heatmap – top 30 variable genes 
top_genes <- names(sort(apply(data, 1, var, na.rm = TRUE), decreasing = TRUE))[1:30]
ann_col   <- data.frame(Group = meta$Group, row.names = colnames(data))

png(file.path(plot_dir, "Heatmap_Top30.png"), width = 4000, height = 2600, res = 300)
pheatmap(data[top_genes, ],
         scale = "row",
         annotation_col = ann_col,
         color = colorRampPalette(c("blue","white","red"))(100),
         show_rownames = TRUE,
         show_colnames = TRUE,       # shows sample names
         fontsize_col = 8,
         main = "Top 30 Most Variable Genes - GSE114419")
dev.off()

# 15. PCA plot
pca_res <- prcomp(t(data), scale. = TRUE)
pct     <- round(100*summary(pca_res)$importance[2,1:2], 1)

pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2],
                     Group = meta$Group, Sample = colnames(data))

png(file.path(plot_dir, "PCA.png"), width = 3400, height = 2400, res = 300)
ggplot(pca_df, aes(PC1, PC2, color = Group, label = Sample)) +
  geom_point(size = 5) + geom_text_repel(size = 3.5) +
  scale_color_manual(values = c("Control" = "#00BFC4", "PCOS" = "#F8766D")) +
  labs(title = "PCA Plot",
       x = paste0("PC1 (", pct[1], "%)"),
       y = paste0("PC2 (", pct[2], "%)")) +
  theme_minimal(base_size = 14)
dev.off()

# 16. Scree plot
var_all <- 100 * summary(pca_res)$importance[2, ]
png(file.path(plot_dir, "Scree.png"), width = 3200, height = 2200, res = 300)
barplot(var_all[1:10], names.arg = paste0("PC",1:10),
        col = "steelblue", las = 2, ylab = "% Variance Explained",
        main = "Scree Plot")
dev.off()

# 17. Batch effect check + ComBat
batch_col <- grep("batch|run|plate|date", names(meta), ignore.case = TRUE, value = TRUE)
batch <- if (length(batch_col)>0) factor(meta[[batch_col[1]]]) else factor(rep(1, ncol(data)))

pca_before <- prcomp(t(data), scale. = TRUE)
df_before  <- data.frame(PC1 = pca_before$x[,1], PC2 = pca_before$x[,2],
                         Group = meta$Group, Batch = batch)
png(file.path(plot_dir, "Batch_Before.png"), width = 3000, height = 2200, res = 300)
ggplot(df_before, aes(PC1, PC2, color = Group, shape = Batch)) +
  geom_point(size = 5) + theme_minimal() +
  scale_color_manual(values = c("Control" = "#00BFC4", "PCOS" = "#F8766D")) +
  labs(title = "PCA Before Batch Correction")
dev.off()

needs_combat <- nlevels(batch) > 1

if (needs_combat) {
  mod <- model.matrix(~ Group, data = meta)
  data_combat <- ComBat(dat = data, batch = batch, mod = mod)
  write.csv(data_combat, file.path(processed_dir, "GSE114419_ComBatCorrected.csv"))
  
  pca_after <- prcomp(t(data_combat), scale. = TRUE)
  df_after  <- data.frame(PC1 = pca_after$x[,1], PC2 = pca_after$x[,2],
                          Group = meta$Group, Batch = batch)
  png(file.path(plot_dir, "Batch_After.png"), width = 3000, height = 2200, res = 300)
  ggplot(df_after, aes(PC1, PC2, color = Group, shape = Batch)) +
    geom_point(size = 5) + theme_minimal() +
    scale_color_manual(values = c("Control" = "#00BFC4", "PCOS" = "#F8766D")) +
    labs(title = "PCA After ComBat")
  dev.off()
  cat("Batch correction applied\n")
} else cat("No batch effect – skipped correction\n")

# 18. Conclusion
conclusion <- paste0(
  "GSE114419 analysis finished!\n",
  "Samples: ", ncol(data), "\n",
  "Outliers flagged: ", length(outliers), "\n",
  "Batch levels found: ", nlevels(batch), " → correction ", ifelse(needs_combat, "applied", "not needed"), "\n",
  "All files saved in QC/", gse_id, " and processed_data/", gse_id
)
writeLines(conclusion, file.path(plot_dir, "Conclusion.txt"))
cat(conclusion)