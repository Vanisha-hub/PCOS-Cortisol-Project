raw_path <- file.path("1_RawData")
output_path <- file.path("3_Outputs")

req_packages <- c(
  "DESeq2", "readr", "dplyr", "tibble", "ggplot2", 
  "ggrepel", "pheatmap", "RColorBrewer"
)
bioc_packages <- c("DESeq2") 

for (pkg in req_packages) {
  library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE) 
}

# Import Metadata
colData_raw <- read_csv(file.path(raw_path, "Master_Sample_List.csv"), show_col_types = FALSE)

colData_clean <- colData_raw |>
  filter(Dataset == "GSE277906") |>
  select(Sample_ID, Sample_title, Condition = Phenotype) |>
  rename(Geo_Accession = Sample_ID) |> 
  column_to_rownames(var = "Sample_title")

colData_clean$Condition <- factor(colData_clean$Condition)
colData_clean$Condition <- relevel(colData_clean$Condition, ref = "Normal")
cat(paste(nrow(colData_clean), "samples found in metadata."))

# Import Count Data and Alignment
countData <- read.delim(
  file = file.path(raw_path, "GSE277906_counts_anno.txt.gz"),
  sep = "\t",
  header = TRUE,
  row.names = 1
)
countData <- as.matrix(countData)
cat(paste(nrow(countData), "genes are found."))

if (any(is.na(countData))) {
  stop("NA values detected in count data.")
}

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
    design = ~ Condition 
  )
  print(dds)
  saveRDS(dds, file.path(output_path, "dds_raw.rds"))
}

# FILTERING AND NORMALIZATION
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(paste(nrow(dds), "genes are filtered."))

# DESeq2 Analysis
dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file = file.path(output_path, "DESeq2_normalized_counts.csv"), row.names = TRUE)

# VST
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
write.csv(vst_mat, file = file.path(output_path, "vst_matrix.csv"), row.names = TRUE)


# QC PLOTS
generate_qc_plots <- function(dds_object, vsd_object, output_path) {
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
  
  # PCA Plot
  pca_data <- plotPCA(vsd_object, intgroup = c("Condition"), returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = Condition)) + 
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = Geo_Accession), size = 3, max.overlaps = 50) + 
    xlab(paste0("PC1: ", percentVar[1], "% of variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% of variance")) +
    labs(
      title = "QC: Principal Component Analysis (PCA)",
      subtitle = "PCOS vs. Normal Cumulus Cells"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.title = element_text(face = "bold"), 
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  
  ggsave(file.path(output_path, "PCA_plot.png"), plot = pca_plot, width = 8, height = 7)
}
generate_qc_plots(dds_object = dds, vsd_object = vsd, output_path = output_path)