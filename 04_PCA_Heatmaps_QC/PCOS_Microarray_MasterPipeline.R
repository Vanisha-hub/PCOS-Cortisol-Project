# ================================================================
# PCOS Microarray Master Pipeline (All Datasets Automated)
# ================================================================

# ----- Load Required Packages -----
packages <- c("GEOquery", "limma", "affy", "oligo", "arrayQualityMetrics",
              "AnnotationDbi", "org.Hs.eg.db", "ggplot2", "pheatmap",
              "factoextra", "ggrepel", "sva", "dplyr", "reshape2")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
  library(pkg, character.only = TRUE)
}

# ----- Journal-style theme -----
theme_journal <- theme_minimal(base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# ----- Dataset Configuration -----
datasets <- list(
  list(id = "GSE98421", platform = "Affymetrix", annot_pkg = "hgu133plus2.db"),
  list(id = "GSE114419", platform = "Affymetrix", annot_pkg = "hta20transcriptcluster.db"),
  list(id = "GSE137684", platform = "Agilent", annot_pkg = "hgug4112a.db")
)

# ================================================================
# MASTER FUNCTION
# ================================================================
process_dataset <- function(ds) {
  gse_id <- ds$id
  platform <- ds$platform
  annot_pkg <- ds$annot_pkg
  cat("\n==============================\nProcessing:", gse_id, "| Platform:", platform, "\n")
  
  # ---- Create Directories ----
  out_dir <- file.path("Processed_Data", gse_id)
  plot_dir <- file.path(out_dir, "Plots")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ------------------------------------------------
  # STEP 1: Load and Normalize Expression Data
  # ------------------------------------------------
  if (platform == "Affymetrix") {
    cel_path <- file.path("Raw_data", gse_id)
    cel_files <- list.celfiles(cel_path, full.names = TRUE)
    if (length(cel_files) == 0) stop(paste("No CEL files in", cel_path))
    
    cat("Reading CEL files for", gse_id, "\n")
    raw_data <- read.celfiles(cel_files)
    
    cat("Performing RMA normalization...\n")
    norm_data <- rma(raw_data)
    expr_data <- exprs(norm_data)
    
  } else if (platform == "Agilent") {
    tar_path <- file.path("Raw_data", gse_id)
    files <- list.files(tar_path, pattern = "txt.gz$", full.names = TRUE)
    if (length(files) == 0) stop(paste("No .txt.gz files in", tar_path))
    
    cat("Reading Agilent TXT files for", gse_id, "\n")
    raw_data <- limma::read.maimages(
      files,
      source = "agilent",
      green.only = TRUE,
      columns = list(G = "gMedianSignal", Gb = "gBGMedianSignal"),
      annotation = c("ProbeName", "GeneName")
    )
    
    cat("Performing background correction and normalization...\n")
    raw_data$E <- limma::backgroundCorrect(raw_data$E, method = "normexp")
    norm_data <- limma::normalizeBetweenArrays(raw_data, method = "quantile")
    
    expr_data <- norm_data$E
    expr_data[expr_data <= 0 | is.na(expr_data)] <- NA
    expr_data <- log2(expr_data)
    
    # Remove incomplete rows/cols, impute medians
    expr_data <- expr_data[rowSums(is.na(expr_data)) < ncol(expr_data) * 0.5, ]
    expr_data <- expr_data[, colSums(is.na(expr_data)) < nrow(expr_data) * 0.5]
    expr_data <- t(apply(expr_data, 1, function(x) {
      x[is.na(x)] <- median(x, na.rm = TRUE)
      x
    }))
    
    if (!is.null(raw_data$genes$ProbeName)) {
      rownames(expr_data) <- raw_data$genes$ProbeName
    } else if (!is.null(raw_data$genes$GeneName)) {
      rownames(expr_data) <- raw_data$genes$GeneName
    } else {
      stop("No probe or gene names found in Agilent file headers.")
    }
  }
  
  write.csv(expr_data, file.path(out_dir, paste0(gse_id, "_NormalizedData.csv")))
  
  # ------------------------------------------------
  # STEP 2: Probe → Gene Mapping
  # ------------------------------------------------
  suppressMessages(library(annot_pkg, character.only = TRUE))
  probe_ids <- rownames(expr_data)
  gene_symbols <- mapIds(get(annot_pkg),
                         keys = probe_ids, keytype = "PROBEID",
                         column = "SYMBOL", multiVals = "first")
  
  expr_df <- data.frame(PROBEID = probe_ids,
                        SYMBOL = gene_symbols, expr_data, check.names = FALSE)
  expr_df <- expr_df %>% filter(!is.na(SYMBOL))
  expr_avg <- limma::avereps(expr_df[, -(1:2)], ID = expr_df$SYMBOL)
  data <- as.data.frame(expr_avg)
  write.csv(data, file.path(out_dir, paste0(gse_id, "_Processed.csv")))
  
  # ------------------------------------------------
  # STEP 3: QC Plots (Boxplot, Density)
  # ------------------------------------------------
  png(file.path(plot_dir, "Boxplot.png"), width = 3400, height = 2400, res = 300)
  par(mar = c(18, 5, 4, 2))  # Extended margin for long x labels
  boxplot(data, main = paste("Boxplot -", gse_id),
          las = 2, col = "skyblue", ylab = "Expression Intensity",
          cex.axis = 0.8, cex.names = 0.7, xaxt = "n")
  axis(1, at = 1:ncol(data), labels = colnames(data), las = 2, cex.axis = 0.7)
  dev.off()
  
  png(file.path(plot_dir, "DensityPlot.png"), width = 3200, height = 2200, res = 300)
  matplot(density(data[, 1], na.rm = TRUE)$x,
          apply(data, 2, function(x) density(x, na.rm = TRUE)$y),
          type = "l", col = rainbow(ncol(data)),
          xlab = "Expression", ylab = "Density",
          main = paste("Density Plot -", gse_id))
  dev.off()
  
  # ------------------------------------------------
  # STEP 4: PCA & Scree Plot
  # ------------------------------------------------
  pca_res <- prcomp(t(data), scale. = TRUE)
  var_explained <- summary(pca_res)$importance[2, 1:2] * 100
  groups <- factor(rep(c("Control", "PCOS"), length.out = ncol(data)))
  pca_df <- data.frame(pca_res$x[, 1:2], Group = groups)
  
  png(file.path(plot_dir, "PCA.png"), width = 3400, height = 2400, res = 300)
  print(
    ggplot(pca_df, aes(PC1, PC2, color = Group)) +
      geom_point(size = 4) +
      geom_text_repel(aes(label = rownames(pca_df)), size = 2.8, max.overlaps = 8) +
      labs(
        title = paste("PCA of PCOS vs Control (", gse_id, ")", sep = ""),
        x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(var_explained[2], 1), "%)")
      ) +
      scale_color_manual(values = c("Control" = "#1F77B4", "PCOS" = "#D62728")) +
      theme_journal
  )
  dev.off()
  
  png(file.path(plot_dir, "ScreePlot.png"), width = 3200, height = 2200, res = 300)
  print(fviz_eig(pca_res, addlabels = TRUE,
                 barfill = "#3182bd", barcolor = "black") +
          ggtitle("Variance Explained by Principal Components"))
  dev.off()
  
  # ------------------------------------------------
  # STEP 5: Heatmap & Clustering
  # ------------------------------------------------
  top_genes <- names(sort(apply(data, 1, var), decreasing = TRUE))[1:30]
  annotation_col <- data.frame(Group = groups)
  rownames(annotation_col) <- colnames(data)
  
  png(file.path(plot_dir, "Heatmap.png"), width = 3200, height = 2200, res = 300)
  pheatmap(data[top_genes, ], scale = "row",
           annotation_col = annotation_col,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = paste("Top 30 Most Variable Genes (", gse_id, ")", sep = ""),
           fontsize_row = 6, fontsize_col = 9)
  dev.off()
  
  png(file.path(plot_dir, "Hierarchical_Clustering.png"), width = 3200, height = 2200, res = 300)
  plot(hclust(dist(t(data)), method = "complete"),
       labels = groups, main = paste("Hierarchical Clustering (", gse_id, ")", sep = ""))
  dev.off()
  
  # ------------------------------------------------
  # STEP 6: Batch Effect Check & QC Summary
  # ------------------------------------------------
  batch_correction_needed <- FALSE
  if ("batch" %in% colnames(data)) {
    mod <- model.matrix(~ groups)
    combat_data <- ComBat(dat = as.matrix(data), batch = data$batch, mod = mod)
    write.csv(combat_data, file.path(out_dir, paste0(gse_id, "_BatchCorrected.csv")))
    batch_correction_needed <- TRUE
    cat("Batch correction applied for", gse_id, "\n")
  } else {
    cat("No batch info found for", gse_id, "\n")
  }
  
  qc_summary <- list(
    dataset = gse_id,
    total_samples = ncol(data),
    potential_outliers = paste(names(which(dist(t(data)) > mean(dist(t(data))) + 2 * sd(dist(t(data))))), collapse = ", "),
    top_variable_genes_used = length(top_genes),
    batch_correction_applied = batch_correction_needed
  )
  
  capture.output(qc_summary, file = file.path(out_dir, paste0("QC_Summary_", gse_id, ".txt")))
  cat("QC summary saved for", gse_id, "\n==============================\n")
}

# ================================================================
# RUN ALL DATASETS
# ================================================================
for (ds in datasets) {
  process_dataset(ds)
}
