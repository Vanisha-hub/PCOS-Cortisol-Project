output_path <- file.path("3_Outputs")

req_packages <- c("readr", "dplyr", "tibble", "ggplot2", "ggrepel", "scales")
bioc_packages <- c("DESeq2", "apeglm")
all_packages <- c(req_packages, bioc_packages)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", dependencies = TRUE)
}

for (pkg in all_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) {
    if (pkg %in% bioc_packages) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
    } else {
      install.packages(pkg, dependencies = TRUE, quiet = TRUE)
    }
  }
  library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
}

# Define thresholds
FDR_THRESHOLD <- 0.05
LFC_THRESHOLD <- 1.0

# LOAD DATA AND RUN DESEQ()
tryCatch({
  dds <- readRDS(file.path(output_path, "dds_raw.rds"))
  cat("Successfully loaded 'dds' object.")
}, error = function(e) {
  stop(paste("ERROR: Failed to load 'dds_raw.rds'. Check path:", output_path))
})

dds <- DESeq(dds)
cat("DESeq() complete.")

# LFC Shrinkage
res_shrunk <- lfcShrink(dds, coef = "Condition_PCOS_vs_Normal", type = "apeglm")

res_shrunk_df <- as.data.frame(res_shrunk) |>
  rownames_to_column(var = "Gene_ID") |>
  mutate(
    Is_Significant = padj < FDR_THRESHOLD & !is.na(padj),
    Gene_Class = case_when(
      Is_Significant & log2FoldChange >= LFC_THRESHOLD ~ "Upregulated (FDR < 0.05 & LFC >= 1)",
      Is_Significant & log2FoldChange <= -LFC_THRESHOLD ~ "Downregulated (FDR < 0.05 & LFC <= -1)",
      TRUE ~ "Not Significant"
    )
  )

# Filtering and saving results
sig_degs <- res_shrunk_df |> filter(Gene_Class != "Not Significant")

write_csv(res_shrunk_df, file = file.path(output_path, "DESeq2_full_results_GSE277906.csv"))
write_csv(sig_degs, file = file.path(output_path, "DESeq2_sigDEGs_GSE277906.csv"))

upregulated <- sum(sig_degs$log2FoldChange > 0)
downregulated <- sum(sig_degs$log2FoldChange < 0)

cat("DEG SUMMARY (LFC|", LFC_THRESHOLD, " & FDR|", FDR_THRESHOLD, ")")
cat("Upregulated genes:  ", upregulated)
cat("Downregulated genes:", downregulated)
cat("Total significant DEGs:", nrow(sig_degs))

# Volcano Plot
genes_to_label <- sig_degs |> arrange(padj) |> head(10)

volcano_plot_ggplot <- ggplot(res_shrunk_df,
                              aes(x = log2FoldChange, y = -log10(padj), color = Gene_Class)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_point(data = genes_to_label, size = 3, alpha = 0.8) +
  geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), linetype = "dashed", color = "black", linewidth = 0.5) +
  
  geom_text_repel(data = genes_to_label,
                  aes(label = Gene_ID),
                  size = 4,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = 'grey50',
                  max.overlaps = Inf) +
  
  scale_color_manual(values = c("Downregulated (FDR < 0.05 & LFC <= -1)" = "blue",
                                "Upregulated (FDR < 0.05 & LFC >= 1)" = "red",
                                "Not Significant" = "grey50")) +
  labs(title = "Volcano Plot — PCOS vs Normal (GSE277906)",
       subtitle = paste("Filtered: |LFC| >", LFC_THRESHOLD, " & FDR <", FDR_THRESHOLD),
       x = expression(log[2]~"Fold Change"),
       y = expression(-log[10]~adjusted~italic("P")~value),
       color = "Status") +
  theme_minimal(base_size = 14)

ggsave(file.path(output_path, "Volcano_Plot.png"), plot = volcano_plot_ggplot, width=10, height=8)

# MA Plot
ma_plot_ggplot <- ggplot(res_shrunk_df,
                         aes(x = baseMean + 1, y = log2FoldChange, color = Gene_Class)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_x_log10(labels = scales::label_number()) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  
  geom_point(data = sig_degs, size = 2, alpha = 0.9) +
  
  scale_color_manual(values = c("Downregulated (FDR < 0.05 & LFC <= -1)" = "blue",
                                "Upregulated (FDR < 0.05 & LFC >= 1)" = "red",
                                "Not Significant" = "grey50")) +
  labs(title = "MA Plot — PCOS vs Normal (GSE277906)",
       x = "Mean of Normalized Counts (+1)",
       y = expression(log[2]~"Fold Change"),
       color = "Status") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(file.path(output_path, "MA_Plot.png"), plot = ma_plot_ggplot, width=8, height=8)