# Task 2 - Microarray Preprocessing

## Folder Contains
- **GSE137684/**  
- **GSE114419/**  
- **GSE98421/**  
**ALL THE NEW FILES (QC AND NORMALISED) ARE IN RESPECTIVE DATASET FOLDER**
-------------------------------------------------------

## Procedure
1. Performed **RMA background correction**.  
2. Applied **quantile normalization**.  
3. Performed **log2 transformation** when necessary.  
4. Mapped probe IDs to gene symbols using annotation databases.  
5. Removed probes with no gene mapping.  
6. Averaged multiple probes mapping to same gene.  
7. Generated QC plots before/after normalization.

----------------------------------------------------------

## Results
- Clean gene-level matrices for all microarray datasets.
- Visual QC confirms consistent normalization.
- All datasets ready for `limma` DEG analysis.

----------------------------------------------------------

## Remarks
- Store R scripts inside each dataset folder.  
- Maintain platform-specific annotation notes.

