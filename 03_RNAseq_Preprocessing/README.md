# Task 3 — RNA-seq Preprocessing (GSE277906)

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

