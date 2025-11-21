raw_path <- file.path("1_RawData")
output_path <- file.path("3_Outputs")

req_packages <- c(
  "DESeq2", "readr", "dplyr", "tibble", "ggplot2",
  "ggrepel", "pheatmap", "RColorBrewer", "rlang"
)
bioc_packages <- c("DESeq2")

if (length(bioc_packages) > 0 && !requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in req_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) {
    if (pkg %in% bioc_packages) {
      cat(paste("Installing Bioconductor package:", pkg))
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      cat(paste("Installing CRAN package:", pkg))
      install.packages(pkg, dependencies = TRUE)
    }
  }
  library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
}

# Import Metadata
colData_raw <- read_csv(file.path(raw_path, "Master_Sample_List.csv"), show_col_types = FALSE)

colData_clean <- colData_raw |>
  filter(Dataset == "GSE277906") |>
  select(Geo_Accession = Sample_ID, Sample_title, Condition = Phenotype) |>
  column_to_rownames(var = "Sample_title")

num_samples <- nrow(colData_clean)

# Batch Assignment (Currently Artificial - Replace with real data if available)
colData_clean$Batch <- factor(rep(c("Batch_A", "Batch_B"), length.out = num_samples))

colData_clean$Condition <- factor(colData_clean$Condition)
colData_clean$Condition <- relevel(colData_clean$Condition, ref = "Normal")
colData_clean$Batch <- factor(colData_clean$Batch)

cat(paste(nrow(colData_clean), "samples are found.")) 

batch_condition_table <- table(colData_clean$Batch, colData_clean$Condition)
print(batch_condition_table)
fisher_test <- fisher.test(batch_condition_table)
cat(paste("Fisher's Exact Test p-value:", round(fisher_test$p.value, 4)))

# Import count data and alignment
countData_raw <- read.delim(
  file = file.path(raw_path, "GSE277906_counts_anno.txt.gz"),
  sep = "\t",
  header = TRUE,
  row.names = 1
)
countData <- as.matrix(countData_raw)
cat(paste(nrow(countData), "genes are found.")) 

if (any(is.na(countData))) {
  stop("NA values found in count data.")
}

# Align count data columns and metadata rows
sample_names_in_count_file <- rownames(colData_clean)
countData <- countData[, intersect(colnames(countData), sample_names_in_count_file), drop = FALSE]
valid_samples <- intersect(colnames(countData), rownames(colData_clean))
countData <- countData[, valid_samples]
colData_clean <- colData_clean[valid_samples, ]
countData <- countData[, rownames(colData_clean)]

alignment_check <- all(colnames(countData) == rownames(colData_clean))
cat(paste("Sample Alignment Check:", alignment_check)) 

storage.mode(countData) <- "integer"

if (alignment_check) {
  dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = colData_clean,
    design = ~ Batch + Condition 
  )
  print(dds)
  saveRDS(dds, file.path(output_path, "dds_raw.rds"))
}

# Filtering and Normalization
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(paste(nrow(dds), "genes are filtered.")) 

dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file = file.path(output_path, "DESeq2_normalized_counts.csv"), row.names = TRUE)
write.csv(vst_mat, file = file.path(output_path, "vst_matrix.csv"), row.names = TRUE)

# OUTLIER DETECTION 
cooks_dist_mat <- assays(dds)[["cooks"]]

median_cooks_df <- data.frame(
  Sample_title = colnames(cooks_dist_mat),
  Median_Cooks = apply(cooks_dist_mat, 2, median, na.rm = TRUE)
)

threshold <- 3 * median(median_cooks_df$Median_Cooks, na.rm = TRUE)

outliers_report <- as.data.frame(colData(dds)) |>
  rownames_to_column(var = "Sample_title") |>
  left_join(median_cooks_df, by = "Sample_title") |>
  mutate(Outlier_Flagged = Median_Cooks > threshold) 

outliers <- outliers_report |>
  filter(Outlier_Flagged == TRUE) |>
  select(Condition, Batch, Sample_title, Geo_Accession, Median_Cooks, Outlier_Flagged)

if (nrow(outliers) > 0) {
  cat(paste(nrow(outliers), "Sample flagged by Median Cook's Distance (Threshold:", round(threshold, 4), ").\n"))
  write_csv(outliers, file.path(output_path, "cooks_flagged_samples.csv"))
} else {
  cat("No samples flagged by Median Cook's Distance.")
}

# QC Plots
qc_plots <- function(dds_object, vsd_object, output_path) {
  vst_mat <- assay(vsd_object)
  
  # Mean-Variance Plot
  png(file.path(output_path, "mean_dispersion_plot.png"), width = 800, height = 600)
  plotDispEsts(dds_object, 
               main = "DESeq2 Mean-Dispersion Relationship")
  dev.off()
  
  # Sample Distance Heatmap
  sample_dists <- dist(t(vst_mat))
  sample_dist_mat <- as.matrix(sample_dists)
  
  labels <- paste0(vsd_object$Condition, " | ", vsd_object$Geo_Accession)
  rownames(sample_dist_mat) <- labels
  colnames(sample_dist_mat) <- labels
  
  colors <- colorRampPalette(rev(brewer.pal(9, "Spectral")))(255)
  png(file.path(output_path, "sample_distance_heatmap.png"), width = 900, height = 700)
  pheatmap(sample_dist_mat,
           clustering_distance_rows = "euclidean", 
           clustering_distance_cols = "euclidean",
           col = colors,
           main = "Sample Distance Heatmap (VST Normalized Data)",
           fontsize = 8)
  dev.off()
}

gen_pca_plot <- function(vsd_object, output_path, intgroup_name) {
  pca_data <- plotPCA(vsd_object, intgroup = c(intgroup_name), returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = .data[[intgroup_name]])) + 
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = Geo_Accession), size = 3, max.overlaps = 50) + 
    xlab(paste0("PC1: ", percentVar[1], "% of variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% of variance")) +
    labs(
      title = paste("QC: PCA Colored by", intgroup_name),
      subtitle = "PCOS vs Normal Cumulus Cells"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.title = element_text(face = "bold"), 
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  
  ggsave(file.path(output_path, paste0("PCA_plot_", intgroup_name, ".png")), plot = pca_plot, width = 8, height = 7)
}

qc_plots(dds_object = dds, vsd_object = vsd, output_path = output_path)
gen_pca_plot(vsd_object = vsd, output_path = output_path, intgroup_name = "Condition")
gen_pca_plot(vsd_object = vsd, output_path = output_path, intgroup_name = "Batch")