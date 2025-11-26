##############################################################
# GSE137684 – Agilent single-color microarray
# Using series matrix expression data directly
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

# Folders
gse_id        <- "GSE137684"
meta_file     <- "metadata/GSE137684_series_matrix.txt"
plot_dir      <- file.path("QC", gse_id)
processed_dir <- file.path("processed_data", gse_id)

dir.create(plot_dir,      recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------
# 1. Load expression data from series matrix
# -------------------------------------------------
gse <- getGEO(filename = meta_file)
meta <- pData(gse)             # sample metadata
feature_data <- fData(gse)     # gene/probe info
expr_matrix <- exprs(gse)      # expression data

# Make sample grouping
meta$Group <- ifelse(grepl("Normal|control", meta$characteristics_ch1.2, ignore.case = TRUE),
                     "Control", "PCOS")
meta$Group <- factor(meta$Group)

# Align metadata with expression data
meta <- meta[colnames(expr_matrix), ]
rownames(meta) <- colnames(expr_matrix)
cat("Matched", ncol(expr_matrix), "samples\n")

# -------------------------------------------------
# 2. Optional normalization (quantile)
# -------------------------------------------------
expr_matrix_norm <- normalizeBetweenArrays(expr_matrix, method = "quantile")

# Save normalized expression
write.csv(expr_matrix_norm, file.path(processed_dir, "GSE137684_NormalizedExpression.csv"), row.names = TRUE)

# -------------------------------------------------
# 3. QC plots
# -------------------------------------------------
# Boxplot
png(file.path(plot_dir, "Boxplot.png"), width = 3400, height = 2400, res = 300)
par(mar = c(12, 4, 4, 2))
boxplot(expr_matrix_norm, las = 2,
        col = ifelse(meta$Group == "Control", "#00BFC4", "#F8766D"),
        main = "Boxplot", ylab = "log2 Intensity")
dev.off()

# -------------------------------------------------
# 4. Outlier detection (PCA distance)
# -------------------------------------------------
pca_temp <- prcomp(t(expr_matrix_norm), scale. = TRUE)
distances <- sqrt(rowSums(pca_temp$x[,1:2]^2))
outliers <- names(which(distances > mean(distances) + 3*sd(distances)))
cat("Outliers (flagged only):", ifelse(length(outliers) > 0, paste(outliers, collapse=", "), "None"), "\n")

# Keep all samples
data_clean <- expr_matrix_norm
meta_clean <- meta

# -------------------------------------------------
# 5. Hierarchical clustering
# -------------------------------------------------
if (ncol(data_clean) >= 2) {
  hc <- hclust(dist(t(data_clean)), method = "complete")
  png(file.path(plot_dir, "Clustering.png"), width = 3000, height = 2000, res = 300)
  plot(hc, labels = paste0(colnames(data_clean), " (", meta_clean$Group, ")"),
       main = "Hierarchical Clustering", cex = 0.8)
  dev.off()
}

# -------------------------------------------------
# 6. Heatmap of top 30 variable genes
# -------------------------------------------------
## 6a) Build robust gene symbol vector (handles multiple column names)
# Try to find a column that looks like "gene symbol"
sym_col <- grep("symbol", colnames(feature_data), ignore.case = TRUE, value = TRUE)

# If none found, try some common fallbacks
if (length(sym_col) == 0) {
  sym_col <- grep("gene.*name|description", colnames(feature_data),
                  ignore.case = TRUE, value = TRUE)
}

# Pull symbols in the exact row order of the expression matrix
gene_symbols <- rep(NA_character_, nrow(expr_matrix_norm))
if (length(sym_col) > 0) {
  # keep the first matching column
  sym_col <- sym_col[1]
  # ensure row alignment (some GEO files can differ)
  feature_data <- feature_data[rownames(expr_matrix_norm), , drop = FALSE]
  gene_symbols <- feature_data[[sym_col]]
}

# Clean symbols:
gene_symbols <- sub(" *//.*$", "", gene_symbols)              # keep text before first //
gene_symbols <- trimws(gene_symbols)
missing <- is.na(gene_symbols) | gene_symbols == "" | gene_symbols == "---"
gene_symbols[missing] <- rownames(expr_matrix_norm)

# Make unique so pheatmap can label rows
rownames(expr_matrix_norm) <- make.unique(gene_symbols)

## 6b) Choose top variable genes and build annotation
top_genes <- names(sort(apply(expr_matrix_norm, 1, var, na.rm = TRUE),
                        decreasing = TRUE))[1:30]

ann_col <- data.frame(Group = meta_clean$Group,
                      row.names = colnames(expr_matrix_norm))

labels_col <- colnames(expr_matrix_norm)
labels_col <- sub("\\.CEL(\\.gz)?$", "", labels_col)  # drop .CEL/.CEL.gz if present

## 6c) Heatmap
png(file.path(plot_dir, "Heatmap_Top30.png"), width = 3200, height = 2400, res = 300)
pheatmap(
  expr_matrix_norm[top_genes, ],
  scale = "row",
  annotation_col = ann_col,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  show_rownames = TRUE,         # <-- gene names at right
  show_colnames = TRUE,         # <-- sample IDs at bottom
  labels_col   = labels_col,    # use cleaned sample labels
  fontsize_row = 7,
  fontsize_col = 9,
  main = "Top 30 Most Variable Genes"
)
dev.off()

# -------------------------------------------------
# 7. PCA + Scree plot
# -------------------------------------------------
pca <- prcomp(t(data_clean), scale. = TRUE)
pct <- round(100*summary(pca)$importance[2,1:2], 1)
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                     Group = meta_clean$Group, Sample = colnames(data_clean))

png(file.path(plot_dir, "PCA.png"), width = 3000, height = 2200, res = 300)
ggplot(pca_df, aes(PC1, PC2, color = Group, label = Sample)) +
  geom_point(size = 5) + geom_text_repel(size = 3.5) +
  scale_color_manual(values = c("Control" = "#00BFC4", "PCOS" = "#F8766D")) +
  labs(title = "PCA Plot", x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)")) +
  theme_minimal(base_size = 14)
dev.off()

# Scree plot
png(file.path(plot_dir, "Scree.png"), width = 2800, height = 2000, res = 300)
barplot(100*summary(pca)$importance[2,1:10], names.arg = paste0("PC",1:10),
        col = "steelblue", las = 2, ylab = "% Variance")
dev.off()

# -------------------------------------------------
# 8. Batch effect check + ComBat
# -------------------------------------------------
# 8a) Detect a batch-like column from metadata (you can also set this manually)
batch_col <- grep("batch|run|plate|sentrix|array|slide|date", 
                  names(meta_clean), ignore.case = TRUE, value = TRUE)
if (length(batch_col) == 0) batch_col <- NULL

# Build batch factor (drop NAs)
if (!is.null(batch_col)) {
  batch <- factor(meta_clean[[batch_col[1]]])
} else {
  batch <- factor(rep("batch1", ncol(data_clean)))  # single level placeholder
}

# Ensure Group is a clean factor aligned to data columns
Group <- factor(meta_clean$Group)
# If Group collapsed to a single level (e.g., all Control), use ~1 design
has_two_groups  <- nlevels(Group) >= 2
has_two_batches <- nlevels(batch) >= 2

# 8b) Quick PCA-by-batch BEFORE ComBat (only if ≥2 samples)
if (ncol(data_clean) >= 2) {
  pca_before <- prcomp(t(data_clean), scale. = TRUE)
  df_before  <- data.frame(PC1 = pca_before$x[,1], PC2 = pca_before$x[,2],
                           Group = Group, Batch = batch, Sample = colnames(data_clean))
  png(file.path(plot_dir, "Batch_PCA_Before.png"), width = 3000, height = 2200, res = 300)
  print(
    ggplot(df_before, aes(PC1, PC2, color = Group, shape = Batch, label = Sample)) +
      geom_point(size = 5) +
      theme_minimal() +
      scale_color_manual(values = c("Control" = "#00BFC4", "PCOS" = "#F8766D")) +
      labs(title = "PCA Before Batch Correction")
  )
  dev.off()
}

# 8c) Run ComBat only when it is statistically valid
if (has_two_batches) {
  # Design matrix: use ~Group if two groups, else ~1
  mod <- if (has_two_groups) model.matrix(~ Group) else model.matrix(~ 1)
  # ComBat
  data_combat <- sva::ComBat(dat = as.matrix(data_clean), batch = batch, mod = mod)
  
  # Save corrected matrix
  write.csv(data_combat, file.path(processed_dir, paste0(gse_id, "_ComBatCorrected.csv")))
  
  # PCA AFTER ComBat
  if (ncol(data_combat) >= 2) {
    pca_after <- prcomp(t(data_combat), scale. = TRUE)
    df_after  <- data.frame(PC1 = pca_after$x[,1], PC2 = pca_after$x[,2],
                            Group = Group, Batch = batch, Sample = colnames(data_combat))
    png(file.path(plot_dir, "Batch_PCA_After.png"), width = 3000, height = 2200, res = 300)
    print(
      ggplot(df_after, aes(PC1, PC2, color = Group, shape = Batch, label = Sample)) +
        geom_point(size = 5) +
        theme_minimal() +
        scale_color_manual(values = c("Control" = "#00BFC4", "PCOS" = "#F8766D")) +
        labs(title = "PCA After ComBat")
    )
    dev.off()
  }
  
  message("✅ ComBat ran: batch levels = ", paste(levels(batch), collapse = ", "))
  
} else {
  # No batch correction possible
  message("ℹ️ Skipping ComBat: batch has <2 levels (found: ", nlevels(batch), ").")
  cat("Skipping ComBat: batch has <2 levels.\n",
      file = file.path(plot_dir, "QC_Summary.txt"), append = TRUE)
}

# -------------------------------------------------
# 9. Conclusion
# -------------------------------------------------
conclusion <- paste0(
  "GSE137684 analysis finished!\n",
  "- Samples: ", ncol(data_clean), "\n",
  "- Outliers flagged: ", length(outliers), "\n",
  "- Batch levels found: ", length(unique(batch)), "\n",
  "- All plots saved in: ", plot_dir, "\n",
  "- Normalized data saved in: ", processed_dir
)
writeLines(conclusion, file.path(plot_dir, "Conclusion.txt"))
cat(conclusion)
