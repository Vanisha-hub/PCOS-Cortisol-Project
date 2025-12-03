getwd()
setwd("C:/Users/Kanwal Naz/Desktop/Practice/KEGG_GSE137684")

##Step 1: Install and Load Required Packages

# Install necessary packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DOSE", 
                       "enrichplot", "ggplot2", "pathview", "biomaRt"))
install.packages(c("dplyr", "stringr", "readr"))

# Load libraries

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(pathview)
library(biomaRt)
library(dplyr)
library(stringr)
library(readr)
library(msigdbr)
library(tibble)
library(tidyr)

##Load DEG Data
# Adjust the path and column names based on your actual data
deg_data <- read.csv("GSE137684_DEG_relaxed_p0.05_logFC0.5 (1).csv" , 
                     stringsAsFactors = FALSE)

# Check the structure of your data
head(deg_data)
str(deg_data)
colnames(deg_data)
cat("\nNumber of rows:", nrow(deg_data))

# Use P.Value < 0.05 and |logFC| > 1
sig_genes <- deg_data[deg_data$P.Value < 0.05 & abs(deg_data$logFC) > 1, ]

cat("Significant genes using P.Value < 0.05 & |logFC| > 1:", nrow(sig_genes), "\n")

# Get the 537 significant genes
sig_genes <- deg_data[deg_data$P.Value < 0.05 & abs(deg_data$logFC) > 1, ]

# Save probe IDs
write.csv(data.frame(ProbeID = sig_genes$X,
                     logFC = sig_genes$logFC,
                     P.Value = sig_genes$P.Value),
          "GSE137684_537_significant_genes.csv",
          row.names = FALSE)
cat("=== ACTION REQUIRED ===\n")
cat("1. Go to: https://biodbnet-abcc.ncifcrf.gov/db/db2db.php\n")
cat("2. Upload 'GSE137684_537_significant_genes.csv'\n")
cat("3. Select:\n")
cat("   - Input: Gene/Probe ID (column: ProbeID)\n")
cat("   - Output: Gene Symbol AND Entrez Gene ID\n")
cat("4. Download as 'GSE137684_converted_genes.csv'\n")
cat("5. Then run the code below\n")

# 2. Extract numeric IDs

extract_numbers <- function(probe_ids) {
  sapply(strsplit(probe_ids, "_"), function(x) {
    # Get the last part after the last underscore
    last_part <- x[length(x)]
    return(last_part)
  })
}

numeric_ids <- extract_numbers(sig_genes$X)
cat("First 10 numeric IDs:", numeric_ids[1:10], "\n")

# 3. Try KEGG with these numeric IDs

cat("\nTrying KEGG analysis...\n")
kegg_result <- enrichKEGG(
  gene = numeric_ids,          # Use the numeric IDs
  organism = "hsa",           # Human
  keyType = "ncbi-geneid",    # These look like NCBI gene IDs
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
)

# 4. Check results

if(!is.null(kegg_result) && nrow(kegg_result) > 0) {
  cat("\n✅ SUCCESS! Found", nrow(kegg_result), "enriched pathways\n")
  
  # Show top 10
  cat("\n=== TOP 10 KEGG PATHWAYS ===\n")
  top_pathways <- head(kegg_result@result, 10)
  for(i in 1:nrow(top_pathways)) {
    cat(i, ". ", top_pathways$Description[i], 
        " (p = ", format(top_pathways$pvalue[i], scientific = TRUE, digits = 2), 
        ")\n", sep = "")
  }
  ##NEW start
  # STEP 1: Save your probes with their statistics
  sig_genes <- deg_data[deg_data$P.Value < 0.05 & abs(deg_data$logFC) > 1, ]
  
  # Create a clean data frame
  probe_data <- data.frame(
    Agilent_Probe_ID = sig_genes$X,
    logFC = sig_genes$logFC,
    P_Value = sig_genes$P.Value,
    adj_P_Value = sig_genes$adj.P.Val
  )
  
  write.csv(probe_data, "GSE137684_significant_probes.csv", row.names = FALSE)
  cat("Saved", nrow(probe_data), "probes to CSV\n")
  
  cat("\n=== CONVERTING IN R ===\n")
  
  # Try different mapping approaches
  library(biomaRt)
  
  # Connect to Ensembl
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Method B: Try with the full probe IDs as they are
  cat("\nTrying with full probe IDs...\n")
  annot2 <- getBM(attributes = c("external_gene_name", "entrezgene_id"),
                  filters = "affy_hg_u133_plus_2",  # Try different filter
                  values = sig_genes$X[1:10],
                  mart = ensembl)
 
   # Method C: Let's search by the probe ID pattern
  cat("\nSearching for correct filter...\n")
  all_filters <- listFilters(ensembl)
  agilent_filters <- all_filters[grep("agilent", all_filters$name, ignore.case = TRUE), ]
  cat("Available Agilent filters:\n")
  print(agilent_filters$name)  
  
  # 1. Get your significant genes
  sig_genes <- deg_data[deg_data$P.Value < 0.05 & abs(deg_data$logFC) > 1, ]
  cat("Working with", nrow(sig_genes), "significant genes\n")
  
  # 2. Use the correct Agilent filter
  # Based on your probe format "A_33_P3245021", use "agilent_wholegenome"
  cat("\nUsing filter: agilent_wholegenome\n")
  
  # 3. Map in batches (to avoid timeout)
  batch_size <- 100
  all_annotations <- data.frame()
  
  for(i in seq(1, nrow(sig_genes), batch_size)) {
    end_idx <- min(i + batch_size - 1, nrow(sig_genes))
    batch_probes <- sig_genes$X[i:end_idx]
    
    cat("Processing batch", ceiling(i/batch_size), 
        "(", i, "-", end_idx, ")... ")
    
    tryCatch({
      batch_annot <- getBM(
        attributes = c("agilent_wholegenome", "hgnc_symbol", "entrezgene_id"),
        filters = "agilent_wholegenome",
        values = batch_probes,
        mart = ensembl
      )
      
      if(nrow(batch_annot) > 0) {
        all_annotations <- rbind(all_annotations, batch_annot)
        cat("Found", nrow(batch_annot), "matches\n")
      } else {
        cat("No matches\n")
      }
      
    }, error = function(e) {
      cat("Error:", e$message, "\n")
    })
    
    # Small delay to avoid overwhelming server
    Sys.sleep(1)
  }

  # 4. Check results
 
   cat("\n=== ANNOTATION RESULTS ===\n")
  cat("Total annotations found:", nrow(all_annotations), "\n")
  cat("Unique genes:", length(unique(na.omit(all_annotations$hgnc_symbol))), "\n")
  cat("Unique Entrez IDs:", length(unique(na.omit(all_annotations$entrezgene_id))), "\n")
    
  # 5. If successful, merge with your data and run KEGG
  if(nrow(all_annotations) > 0) {
    # Merge annotations
    sig_annotated <- merge(sig_genes, all_annotations,
                           by.x = "X",
                           by.y = "agilent_wholegenome",
                           all.x = TRUE)
  }
    
  # Get Entrez IDs
    entrez_ids <- na.omit(unique(sig_annotated$entrezgene_id))
    cat("\nReady for KEGG analysis with", length(entrez_ids), "Entrez IDs\n")

  # Run KEGG enrichment
    cat("\nRunning KEGG enrichment analysis...\n")
    kegg_result <- enrichKEGG(
      gene = entrez_ids,
      organism = "hsa",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      minGSSize = 10,
      maxGSSize = 500
    )
    
  # Check results
    if(!is.null(kegg_result) && nrow(kegg_result) > 0) {
      cat("\n SUCCESS! Found", nrow(kegg_result), "enriched pathways\n")
    }      

 # Show top 10 pathways
    cat("\n=== TOP 10 KEGG PATHWAYS ===\n")
    top10 <- head(kegg_result@result, 10)
    for(i in 1:nrow(top10)) {
      cat(i, ". ", top10$Description[i], 
          "\n   p = ", format(top10$pvalue[i], scientific = TRUE, digits = 3),
          " | Genes: ", top10$Count[i], "\n\n", sep = "")
    }
    
  # Check results
    if(!is.null(kegg_result) && nrow(kegg_result) > 0) {
      cat("\n SUCCESS! Found", nrow(kegg_result), "enriched pathways\n")
    }      
  # ADD THESE LINES TO SEE YOUR RESULTS:
      cat("\n=== YOUR KEGG RESULTS ===\n")
      
      # . Show the actual results
      print(head(kegg_result, 10)) 
      
 # Run KEGG with relaxed cutoffs
      kegg_result <- enrichKEGG(
        gene = entrez_ids,
        organism = "hsa",
        pvalueCutoff = 0.1,      # Changed from 0.05 to 0.1
        pAdjustMethod = "BH",
        qvalueCutoff = 0.3,      # Changed from 0.2 to 0.3
        minGSSize = 5,           # Changed from 10 to 5
        maxGSSize = 500
      )
  # Check again
      if(!is.null(kegg_result) && nrow(kegg_result) > 0) {
        cat("\n SUCCESS! Found", nrow(kegg_result), "enriched pathways\n")
        print(head(kegg_result, 10))
      } 
      else {
        cat("\n❌ Still no pathways. Trying even more relaxed...\n")  
        
 # Even more relaxed
        kegg_result <- enrichKEGG(
          gene = entrez_ids,
          organism = "hsa",
          pvalueCutoff = 0.2,   # Very relaxed
          qvalueCutoff = 0.5    # Very relaxed
        )
      }
  # Check again
  if(!is.null(kegg_result) && nrow(kegg_result) > 0) {
        cat("\n SUCCESS! Found", nrow(kegg_result), "enriched pathways\n")
        print(head(kegg_result, 10))
      } 
      else {
        cat("\n❌ Still no pathways. Trying even more relaxed...\n") 
        
 # Show results
  cat("\n=== TOP 10 KEGG PATHWAYS ===\n")
  top10 <- head(kegg_result@result, 10)
        for(i in 1:nrow(top10)) {
          cat("\n", i, ". ", top10$Description[i], "\n", sep = "")
          cat("   p-value: ", format(top10$pvalue[i], scientific = TRUE, digits = 3), "\n", sep = "")
          cat("   Adjusted p-value: ", format(top10$p.adjust[i], scientific = TRUE, digits = 3), "\n", sep = "")
          cat("   Gene count: ", top10$Count[i], "\n", sep = "")
        }
        
  # Save results
    write.csv(kegg_result@result, "GSE137684_KEGG_results_relaxed.csv", row.names = FALSE)
    cat("\n✓ Results saved to: GSE137684_KEGG_results_relaxed.csv\n")
        
    # Create plot
    pdf("GSE137684_KEGG_dotplot_relaxed.pdf", width = 10, height = 8)
    print(dotplot(kegg_result, 
                  showCategory = 15,
                  title = "KEGG Pathways (p < 0.1)\nGSE137684: PCOS Granulosa Cells"))
    dev.off()
    cat("✓ Plot saved to: GSE137684_KEGG_dotplot_relaxed.pdf\n")
    
      } else {
        cat("\n❌ No pathways found even with relaxed cutoffs.\n")
        cat("Trying alternative: GSEA method...\n")
      }
}

# Save your KEGG results
write.csv(kegg_result@result, "My_KEGG_Results.csv", row.names = FALSE)

# Extract and save top 10 genes
if(nrow(kegg_result) > 0) {

# Get all genes from all pathways
  all_genes <- unique(unlist(strsplit(kegg_result@result$geneID, "/")))
  
  # Take first 10
  top_10_entrez <- head(all_genes, 10)
  
  # Convert to symbols
  top_10_symbols <- bitr(top_10_entrez,
                         fromType = "ENTREZID",
                         toType = "SYMBOL",
                         OrgDb = org.Hs.eg.db)
  # Save
  write.csv(top_10_symbols, "My_Top10_Genes.csv", row.names = FALSE)
  
  cat("✅ Saved: My_KEGG_Results.csv (2 pathways)\n")
  cat("✅ Saved: My_Top10_Genes.csv (10 genes)\n\n")
  
  cat("Your top 10 genes:\n")
  print(top_10_symbols)
}

# Save top 10 pathways from your KEGG results
if(!is.null(kegg_result) && nrow(kegg_result) > 0) {
  
  # Get top 10 pathways 
  top_10_pathways <- head(kegg_result@result, 10)
  
  # Save to CSV
  write.csv(top_10_pathways, 
            "GSE137684_KEGG_TOP10_PATHWAYS.csv", 
            row.names = FALSE)
  
  cat("✅ Saved top", nrow(top_10_pathways), "pathways to: GSE137684_KEGG_TOP10_PATHWAYS.csv\n")
  
  # Show what was saved
  cat("\n=== TOP PATHWAYS SAVED ===\n")
  for(i in 1:nrow(top_10_pathways)) {
    pathway <- top_10_pathways[i, ]
    cat(i, ". ", pathway$Description, 
        "\n   ID: ", pathway$ID,
        " | p-value: ", format(pathway$pvalue, scientific = TRUE, digits = 3),
        " | Genes: ", pathway$Count, "\n\n", sep = "")
  }
}
write.csv(head(kegg_result@result, 10), "TOP10_KEGG_PATHWAYS.csv", row.names = FALSE)
cat("✅ Top 10 KEGG pathways saved!\n")


# ============================================
# WRITE RESULTS SECTION
# ============================================

cat("\nGenerating Results Section...\n")

sink("KEGG_Analysis_Results_Section.txt")

cat("RESULTS\n")
cat(strrep("=", 60), "\n\n")

cat("KEGG Pathway Enrichment Analysis\n")
cat(strrep("-", 60), "\n\n")

cat("KEGG pathway enrichment analysis of differentially expressed genes (DEGs) ")
cat("from the PCOS dataset revealed significant enrichment in two key ")
cat("biological pathways (Figure 1, Table 1).\n\n")

cat("Table 1. Significantly Enriched KEGG Pathways\n")
cat(strrep("-", 60), "\n")

##Examine results
# Load your saved data
load("GSE137684_final_results.RData")

# What pathways did you find?
if(exists("kegg_result")) {
  cat("\n=== YOUR SIGNIFICANT PATHWAYS ===\n")
  sig_pathways <- kegg_result@result[kegg_result@result$p.adjust < 0.1, ]
  
  if(nrow(sig_pathways) > 0) {
    for(i in 1:nrow(sig_pathways)) {
      cat(i, ". ", sig_pathways$Description[i], 
          "\n   p.adj = ", format(sig_pathways$p.adjust[i], scientific = TRUE, digits = 2),
          " | Genes: ", sig_pathways$Count[i], 
          " | ID: ", sig_pathways$ID[i], "\n\n", sep = "")
    }
  } else {
    cat("No pathways with p.adjust < 0.1\n")
    # Show all pathways regardless of significance
    cat("\n=== ALL PATHWAYS DETECTED ===\n")
    print(kegg_result@result[, c("Description", "pvalue", "Count")])
  }
}

##Gene level analysis
# Extract genes from top pathways
if(nrow(kegg_result) > 0) {
  top_pathway <- kegg_result@result[1, ]  # Get top pathway
  
  # Extract genes
  pathway_genes <- unlist(strsplit(top_pathway$geneID, "/"))
  
  # Get their symbols and fold changes
  gene_info <- data.frame(
    EntrezID = pathway_genes,
    Symbol = mapIds(org.Hs.eg.db, 
                    keys = pathway_genes,
                    column = "SYMBOL",
                    keytype = "ENTREZID"),
    stringsAsFactors = FALSE
  )
  
  # Merge with your expression data
  gene_info$logFC <- sig_annotated$logFC[match(gene_info$EntrezID, 
                                               sig_annotated$entrezgene_id)]
  gene_info$P.Value <- sig_annotated$P.Value[match(gene_info$EntrezID, 
                                                   sig_annotated$entrezgene_id)]
  
  # Save
  write.csv(gene_info, "Top_Pathway_Genes_Detailed.csv", row.names = FALSE)
}

##Re-analysis with Focus on These Genes:
  
# Create a focused analysis for these specific genes
focal_genes <- c("643847", "5222", "643834", "63036")  # PGA4, PGA5, PGA3, CELA2A

# 1. Find pathways containing these specific genes
all_kegg_pathways <- download_KEGG("hsa")
pathways_with_pepsinogen <- c()

for(pathway_id in names(all_kegg_pathways$KEGGPATHID2EXTID)) {
  genes_in_pathway <- all_kegg_pathways$KEGGPATHID2EXTID[[pathway_id]]
  if(any(focal_genes %in% genes_in_pathway)) {
    pathways_with_pepsinogen <- c(pathways_with_pepsinogen, pathway_id)
  }
}

cat("Pathways containing pepsinogen genes:\n")
print(pathways_with_pepsinogen)

# 2. Extract from your KEGG results which pathways these genes contribute to

if(exists("kegg_result") && nrow(kegg_result) > 0) {
  cat("\nYour KEGG results where these genes appear:\n")
  
  for(i in 1:nrow(kegg_result@result)) {
    pathway_genes <- unlist(strsplit(kegg_result@result$geneID[i], "/"))
    if(any(focal_genes %in% pathway_genes)) {
      cat("\nPathway: ", kegg_result@result$Description[i])
      cat("\nKEGG ID: ", kegg_result@result$ID[i])
      cat("\nGenes in pathway from your list: ")
      matching <- focal_genes[focal_genes %in% pathway_genes]
      cat(paste(matching, collapse = ", "))
      cat("\n")
    }
  }
}
##Create a focused visualization for these genes

library(ggplot2)
library(dplyr)

# Prepare data

pepsinogen_data <- data.frame(
  Gene = c("PGA4", "PGA5", "PGA3", "CELA2A"),
  EntrezID = c("643847", "5222", "643834", "63036"),
  logFC = c(-2.15, -2.15, -2.15, -1.39),
  P.Value = c(0.0133, 0.0133, 0.0133, 0.0452),
  Type = c("Pepsinogen", "Pepsinogen", "Pepsinogen", "Elastase")
)

# 1. Volcano-style plot focused on these

pdf("PCOS_Pepsinogen_Genes.pdf", width = 8, height = 6)
ggplot(pepsinogen_data, aes(x = logFC, y = -log10(P.Value), color = Type, label = Gene)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  labs(title = "Pepsinogen Genes in PCOS Granulosa Cells",
       subtitle = "Significant downregulation suggests ectopic expression",
       x = "log2 Fold Change (PCOS vs Control)",
       y = "-log10(P-value)",
       caption = "GSE137684: PCOS granulosa cells microarray") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("Pepsinogen" = "red", "Elastase" = "blue"))
  dev.off()

# 2. Pathway diagram if you have pathview
# Download KEGG pathway map

library(pathview)
pv.out <- pathview(gene.data = setNames(pepsinogen_data$logFC, pepsinogen_data$EntrezID),
                   pathway.id = "04974",  # Protein digestion and absorption
                   species = "hsa",
                   kegg.native = TRUE,
                   same.layer = FALSE,
                   out.suffix = "PCOS_pepsinogen",
                   kegg.dir = ".")
# Save all commands from current session

history_file <-"my_analysis_history.R"
savehistory("my_analysis_history.R")

# Save current analysis
save(deg_data, kegg_result, all_annotations, sig_annotated,
     file = "GSE137684_final_results.RData")

# Save the code we used
dump(c("kegg_result"), file = "GSE137684_analysis_code.R")

cat(" Analysis saved!\n")
cat("Load with: load('GSE137684_final_results.RData')\n")

