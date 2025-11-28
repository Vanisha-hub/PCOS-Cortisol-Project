# Task 3 — RNA-seq Preprocessing (GSE277906)
__________________________________
Editted and Final Files in This Folder
1) GSE277906_counts_anno.txt
Input file (raw counts).
Imported directly from the original GEO dataset.
Contains gene identifiers + count matrix.

2) GSE277906_dds.rds
DESeq2 dataset object (dds).
Includes:
Raw filtered counts
Metadata (PCOS vs Control)
Experimental design
Useful for downstream DEG analysis (Task 6).

3) GSE277906_vsd.rds
Variance Stabilized Transformation (VST) DESeq2 object (vsd).
Used for PCA, clustering, and batch-effect assessment.
Preferred for visualization and QC.

4) GSE277906_vst_matrix.csv
Numerical VST-normalized expression matrix.
Rows = genes, Columns = samples.
Used for heatmaps, PCA, distance matrices, and integrative QC.

5) mean_sd_plot.png
Mean–SD plot for detection of dependence between variance & mean.
Confirms effectiveness of variance stabilization.

6) PCA_vst.png
PCA plot generated from VST-normalized counts.
Shows sample clustering between PCOS vs Control groups.
Used to detect outliers and major variance drivers.

_________________________________
## Folder Contains
- **GSE277906_Raw_Counts/**  
- **DESeq2_Normalized/** — vst, normalized matrices  
- **QC_Plots/** — sample distance, variance plots  
- **Metadata/**  

----------------------------------------------------------

## 🧪 Procedure
1. Imported raw counts using **DESeq2**.  
2. Filtered low-expression genes.  
3. Computed **size factors**.  
4. Generated **VST-transformed** matrices.  
5. Created QC plots:
   - mean–variance trend  
   - sample-to-sample distances  
6. Verified separation of PCOS vs controls.

----------------------------------------------------------

## Results
- Clean, normalized RNA-seq expression matrix.
- VST data ready for PCA, DEG, and WGCNA.

----------------------------------------------------------

## 📝 Remarks
- Always store raw count files unmodified.
- Reproduce VST using provided scripts.

