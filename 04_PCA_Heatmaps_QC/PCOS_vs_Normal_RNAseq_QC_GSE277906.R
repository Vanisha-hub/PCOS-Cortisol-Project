# -------------------------------
# 1. Load Required Packages
# -------------------------------

cran_pkgs <- c("tidyverse", "ggrepel", "pheatmap", "factoextra")
bioc_pkgs <- c("GEOquery", "DESeq2", "Biobase", "sva")

install_if_missing <- function(pkgs, installer) {
  to_install <- setdiff(pkgs, rownames(installed.packages()))
  if (length(to_install)) installer(to_install, ask = FALSE)
}
install_if_missing(cran_pkgs, install.packages)
install_if_missing(bioc_pkgs, BiocManager::install)

invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))

# -------------------------------
# 2. Define Paths
# -------------------------------
raw_path <- "Raw_data/GSE277906/GSE277906_counts_anno.txt"
output_path <- "processed_data/GSE277906/"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# 3. Load Metadata from GEO
# -------------------------------
gse <- getGEO("GSE277906", GSEMatrix = TRUE)[[1]]
meta <- pData(gse)
meta$geo_accession <- rownames(meta)

# Define PCOS vs Normal groups
meta$Condition <- ifelse(grepl("PCOS", meta$title, ignore.case = TRUE), "PCOS", "Normal")
meta$Condition <- factor(meta$Condition, levels = c("Normal", "PCOS"))

# -------------------------------
# 4. Load and Clean Count Data
# -------------------------------
raw_df <- read.delim(raw_path, sep = "\t", header = TRUE, check.names = FALSE)
rownames(raw_df) <- raw_df[, 1]
countData <- raw_df[, -1]

# Investigate NA structure
na_summary <- colSums(is.na(countData))
if (any(na_summary > 0)) {
  message("Samples with missing values:")
  print(na_summary[na_summary > 0])
}
countData <- countData[rowSums(is.na(countData)) == 0, ]

# -------------------------------
# 5. Match Samples by geo_accession
# -------------------------------
colnames(countData) <- gsub("\\s+", "", colnames(countData))
meta$title_clean <- gsub("\\s+", "", meta$title)
match_ids <- intersect(colnames(countData), meta$title_clean)
if (length(match_ids) == 0) stop("No matching sample IDs found.")

meta <- meta[meta$title_clean %in% match_ids, ]
countData <- countData[, match_ids]
stopifnot(all(colnames(countData) == meta$title_clean))

# -------------------------------
# 6. Create DESeq2 Dataset
# -------------------------------
rownames(meta) <- meta$title_clean
dds <- DESeqDataSetFromMatrix(countData = countData, colData = meta, design = ~ Condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- estimateSizeFactors(dds)

norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, paste0(output_path, "Normalized_Counts.csv"))

vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
write.csv(vst_mat, paste0(output_path, "VST_Matrix.csv"))

# -------------------------------
# 7. PCA and Outlier Detection
# -------------------------------
pca_res <- prcomp(t(vst_mat), scale. = TRUE)
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)))[1:2]
pca_df <- data.frame(pca_res$x[, 1:2], Condition = meta$Condition, Sample = rownames(meta))

# Mahalanobis distance for outlier detection
mahal_dist <- mahalanobis(pca_res$x[, 1:3], colMeans(pca_res$x[, 1:3]), cov(pca_res$x[, 1:3]))
outliers <- names(mahal_dist[mahal_dist > quantile(mahal_dist, 0.975)])
writeLines(outliers, paste0(output_path, "Potential_Outliers.txt"))

# PCA plot
p1 <- ggplot(pca_df, aes(PC1, PC2, color = Condition, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  theme_classic(base_size = 14) +
  labs(title = "PCA of PCOS vs Normal",
       x = paste0("PC1 (", percentVar[1], "%)"),
       y = paste0("PC2 (", percentVar[2], "%)")) +
  scale_color_manual(values = c("Normal" = "#1F77B4", "PCOS" = "#D62728"))
ggsave(paste0(output_path, "PCA_Plot.png"), plot = p1, width = 8, height = 6)

# -------------------------------
# 8. Heatmap of Top Variable Genes
# -------------------------------
top_genes <- names(sort(apply(vst_mat, 1, var), decreasing = TRUE))[1:30]
heatmap_data <- vst_mat[top_genes, ]
annotation_col <- data.frame(Condition = meta$Condition)
rownames(annotation_col) <- colnames(heatmap_data)

pheatmap(heatmap_data,
         scale = "row",
         annotation_col = annotation_col,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "correlation",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top 30 Most Variable Genes",
         filename = paste0(output_path, "TopVariableGenes_Heatmap.png"),
         fontsize_row = 6,
         fontsize_col = 8)

# -------------------------------
# 9. Batch Effect Assessment
# -------------------------------
if ("batch" %in% colnames(meta)) {
  pca_df$Batch <- meta$batch
  p_batch <- ggplot(pca_df, aes(PC1, PC2, color = Batch)) +
    geom_point(size = 3) +
    theme_classic(base_size = 14) +
    labs(title = "PCA Colored by Batch")
  ggsave(paste0(output_path, "PCA_by_Batch.png"), plot = p_batch, width = 8, height = 6)
}

# -------------------------------
# 10. Batch Correction (Optional)
# -------------------------------
if ("batch" %in% colnames(meta)) {
  mod <- model.matrix(~ Condition, data = meta)
  vst_corrected <- ComBat(dat = vst_mat, batch = meta$batch, mod = mod)
  vst_mat <- vst_corrected
  write.csv(vst_mat, paste0(output_path, "VST_Matrix_BatchCorrected.csv"))
}

# -------------------------------
# 11. QC Summary
# -------------------------------
qc_summary <- list(
  dataset = "GSE277906",
  total_samples = ncol(vst_mat),
  group_distribution = table(meta$Condition),
  pca_variance_PC1 = percentVar[1],
  pca_variance_PC2 = percentVar[2],
  outliers = outliers,
  top_variable_genes = top_genes,
  batch_effect_checked = "batch" %in% colnames(meta)
)

saveRDS(qc_summary, file = paste0(output_path, "QC_Summary.rds"))
writeLines(capture.output(print(qc_summary)), paste0(output_path, "QC_Summary.txt"))

cat("GSE277906 RNA-seq preprocessing and QC completed.\n")