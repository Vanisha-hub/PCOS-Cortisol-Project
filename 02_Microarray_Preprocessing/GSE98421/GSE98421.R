gc()
library(limma)
library(dplyr)
library(hgu133plus2.db)
library(AnnotationDbi)
install.packages("org.Hs.eg.db")
remove.packages("org.Hs.eg.db")
BiocManager::install("org.Hs.eg.db", version = "3.21", ask = FALSE)
library(org.Hs.eg.db)
BiocManager::install("AnnotationDbi", ask = FALSE, update = TRUE, force = TRUE)
packageVersion("AnnotationDbi")
packageVersion("org.Hs.eg.db")
gse_data_a <- getGEO("GSE98421", GSEMatrix = TRUE)
expression_data_a <- exprs(gse_data_a$GSE98421_series_matrix.txt.gz)
feature_data_a <- fData(gse_data_a$GSE98421_series_matrix.txt.gz)
phenotype_data_a <- pData(gse_data_a$GSE98421_series_matrix.txt.gz)
sum(is.na(phenotype_data_a$source_name_ch1))

untar("Raw_Data/GSE98421_RAW.tar", exdir = "Raw_Data/CEL_Files_1")
raw_data_a <- ReadAffy(celfile.path = "Raw_Data/CEL_Files_1")
raw_data_a

arrayQualityMetrics(expressionset = raw_data_a,
                    outdir = "Results/QC_Raw_Data_1",
                    force = TRUE,
                    do.logtransform = TRUE)
normalized_data_a <- rma(raw_data_a)
arrayQualityMetrics(expressionset = normalized_data_a,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)

processed_data_a <- as.data.frame(exprs(normalized_data_a))
dim(processed_data_a)

row_median_a <- rowMedians(as.matrix(processed_data_a))
hist(row_median_a,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")
threshold <- 5.7
abline(v = threshold, col = "orange", lwd = 4)
indx <- row_median_a > threshold
filtered_data <- processed_data_a [indx, ]
colnames(filtered_data) <- rownames(phenotype_data_a)
processed_data_a <- filtered_data

probe_ids <- rownames(processed_data_a)

gene_symbols <- mapIds(
 hgu133plus2.db,
 keys = probe_ids,
 keytype = "PROBEID",
 column = "SYMBOL",
 multiVals = "first"
)

symbols <- AnnotationDbi::select(hgu133plus2.db,
                                 keys = probe_ids,
                                 keytype = "PROBEID",
                                 columns = c("SYMBOL", "ENTREZID", "GENENAME"))

gene_map_df <- gene_symbols %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::rename(SYMBOL = 2)

duplicate_summary <- gene_map_df %>%
  group_by(SYMBOL)  %>%
  summarise(probes_per_gene =  n()) %>%
  arrange(desc(probes_per_gene)) 

duplicate_genes <- duplicate_summary %>%
  filter(probes_per_gene > 1)

sum(duplicate_genes$probes_per_gene)
22570 - 2631

all(gene_map_df$PROBEID == row.names(processed_data_a))

processed_data_df <- processed_data_a %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::mutate(SYMBOL = gene_symbols[PROBEID]) %>%
  dplyr::relocate(SYMBOL, .after = PROBEID)

processed_data_df <- processed_data_df %>%
  dplyr::filter(!is.na(SYMBOL))

expr_only <- processed_data_df %>%
  dplyr::select(-PROBEID, -SYMBOL)

average_data <- limma::avereps(expr_only, ID = processed_data_df$SYMBOL)
dim(average_data)

data <- as.data.frame(average_data)
data <- data.matrix(data)
str(data)
is.numeric(data)

write.csv(data, file = "data.csv")
save(feature_data_a, phenotype_data_a, processed_data_a, raw_data_a, data, 
     file = "C:/Users/ASUS/Desktop/PCOS Final Project.RData")

write.table(expr_gene, "GSE98421_DEG_ready.txt",
            sep = "\t", quote = FALSE, row.names = TRUE)