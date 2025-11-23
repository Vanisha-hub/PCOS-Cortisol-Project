# Task 4 — PCA, Heatmaps & Sample Quality Assessment

## Folder Contains
- **Microarray GSE114419/**  
- **Microarray GSE137684/**  
- **Microarray GSE98421/**  
- **RNAseq GSEGSE277906/**  

----------------------------------------------------------

## Dependencies

- BiocManager
- limma
- GEOquery
- sva
- pheatmap
- ggplot2
- dplyr
- ggrepel

----------------------------------------------------------

## Procedure
1. Loaded raw expression data and verified file paths. 
2. Assessed sample grouping (PCOS vs control).  
3. Matched metadata (meta$title_clean) with count matrix column names.
   - Reordered metadata to align with count matrix.
4. Ran PCA on each dataset.  
5. Assessed sample grouping (PCOS vs control).  
6. Generated heatmaps of top variable genes.
7. Checked hierarchical clustering.
8. Identified any outlier samples.
9. Performed batch-effect detection using sample metadata


----------------------------------------------------------

## Results
- QC plots: Boxplot, density plot.
- PCA plots for all datasets.
- Clustering: Hierarchical dendrogram.
- Heatmaps showing global expression variation for all datasets.
- Outlier notes added if any detected.
- Metadata and count matrix alignment confirmed (after cleaning)  

----------------------------------------------------------

## Remarks
- Outliers should be documented but removed ONLY after discussion.
- Always check that colnames(countData) and meta$title_clean match exactly before downstream analysis.
- Case sensitivity, underscores, and ordering can cause mismatches — normalize and reorder as needed.
- final Proofreading made by Deiby Cabuyales

