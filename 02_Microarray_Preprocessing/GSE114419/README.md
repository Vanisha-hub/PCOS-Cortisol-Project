# GSE114419

Subfolder for 02_Microarray_Preprocessing.


library(GEOquery)
gse <- getGEO("GSE114419")
expression_data <- exprs(gse[[1]])
print(dim(expression_data))
feature_data <- fData(gse[[1]])
print(colnames(feature_data))
print("First 10 rows of gene_assignment:")
print(head(feature_data$gene_assignment, 10))

extract_gene_symbol <- function(x) {
  if(is.na(x) || x == "") return(NA)
  parts <- strsplit(x, " // ")[[1]]
  if(length(parts) >= 2) return(parts[2])
  return(NA)
}

gene_symbols <- sapply(feature_data$gene_assignment, extract_gene_symbol)
print("First 10 gene symbols:")
print(head(gene_symbols, 10))

rownames(expression_data) <- gene_symbols

print("Expression data with gene symbols:")
print(dim(expression_data))
print("First 10 gene symbols in expression data:")
print(rownames(expression_data)[1:10])



valid_genes <- !is.na(rownames(expression_data)) & rownames(expression_data) != "---"
expression_data_clean <- expression_data[valid_genes, ]

print("After removing invalid gene symbols:")
print(dim(expression_data_clean))
print("First 10 valid gene symbols:")
print(rownames(expression_data_clean)[1:10])


groups <- factor(c("normal", "normal", "normal", "pcos", "pcos", "pcos"),
                 levels = c("normal", "pcos"))
print("Groups:")
print(groups)


library(limma)


design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
print("Design matrix:")
print(design)


fit <- lmFit(expression_data_clean, design)


contrast_matrix <- makeContrasts(pcos_vs_normal = pcos - normal, levels = design)
print("Contrast matrix:")
print(contrast_matrix)


fit_contrast <- contrasts.fit(fit, contrast_matrix)
fit_ebayes <- eBayes(fit_contrast)

print("eBayes completed. Ready to extract results.")


deg_results <- topTable(fit_ebayes, coef = "pcos_vs_normal", number = Inf, adjust.method = "BH")

print("DEG results summary:")
print(dim(deg_results))
print(head(deg_results))


deg_results$threshold <- ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated", "No")
)


deg_results$threshold <- ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated", "No")
)

print("DEG counts:")
print(table(deg_results$threshold))


deg_results$threshold <- ifelse(
  deg_results$P.Value < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$P.Value < 0.05 & deg_results$logFC < -1, "Downregulated", "No")
)

print("DEG counts (using P.Value < 0.05):")
print(table(deg_results$threshold))

# Check current data distribution (boxplots)
boxplot(expression_data_clean, main = "Before Normalization")

# Create boxplot to show normalized data
boxplot(averaged_data, 
        main = "After Normalization Verification\n(GSE114419 - PCOS vs Normal)",
        xlab = "Samples", 
        ylab = "Expression Level",
        las = 2,
        col = c("lightblue", "lightblue", "lightblue", "pink", "pink", "pink"))


library(limma)
averaged_data <- avereps(expression_data_clean, ID = rownames(expression_data_clean))

#density plot
plot(density(averaged_data[,1]), 
     col = "blue", lwd = 2,
     main = "Density Plot: Normalized Expression Data\n(GSE114419 - PCOS vs Normal)",
     xlab = "Expression Level", 
     ylab = "Density",
     xlim = range(averaged_data),
     ylim = c(0, 0.3))

colors <- rep(c("blue", "red"), each = 3)
for(i in 2:6) {
  lines(density(averaged_data[,i]), 
        col = colors[i], 
        lwd = 1.5)
}

legend("topright", 
       legend = c("Normal (n=3)", "PCOS (n=3)"),
       col = c("blue", "red"), 
       lwd = 2)

print("Summary of expression values:")
print(summary(expression_data_clean[,1]))

 
print("Mean expression per sample:")
print(colMeans(expression_data_clean))


top_degs <- deg_results[deg_results$threshold %in% c("Upregulated", "Downregulated"), ]
top_degs <- top_degs[order(top_degs$P.Value), ]
top_25 <- head(top_degs, 25)

print("Top 25 DEGs:")
print(top_25[, c("logFC", "P.Value", "threshold")])


print("Top 25 DEG row names:")
print(rownames(top_25))

print("Expression data row names (first 10):")
print(rownames(expression_data_clean)[1:10])

print("Checking row name matches:")
missing_genes <- rownames(top_25)[!rownames(top_25) %in% rownames(expression_data_clean)]
print("Missing genes:")
print(missing_genes)

print("Available genes in top_25:")
print(rownames(top_25))

print("Do they match exactly?")
print(rownames(top_25) %in% rownames(expression_data_clean))

top_25_indices <- as.numeric(rownames(top_25))
top_25_genes <- rownames(expression_data_clean)[top_25_indices]

print("Top 25 DEG gene symbols:")
print(top_25_genes)


gene_counts <- table(rownames(expression_data_clean))
duplicate_genes <- gene_counts[gene_counts > 1]

print("Genes with multiple probes:")
print(length(duplicate_genes))
print("Examples:")
print(head(duplicate_genes))


library(limma)
averaged_data <- avereps(expression_data_clean, ID = rownames(expression_data_clean))
print("After averaging multiple probes:")
print(dim(averaged_data))

# Save
write.csv(deg_results, "DEGs_Results.csv")
write.csv(top_degs[top_degs$threshold == "Upregulated", ], "Upregulated_DEGs.csv") 
write.csv(top_degs[top_degs$threshold == "Downregulated", ], "Downregulated_DEGs.csv")

print("DEG files saved!")
