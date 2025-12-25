# PCOS-Cortisol-Project

A collaborative repository for sharing PCOS and cortisol-related research, analysis scripts, datasets, and documentation.  
Contributors are welcome to add their folders, workflows, and results.

# PCOS–Cortisol Project: Task Allocation & Workflow

This repository contains all datasets, scripts, analyses, and documentation for the PCOS–Cortisol multi-omics integration project.  
Below is the complete folder structure and task allocation for all team members.

## 📋 Task Assignment Table

| Sno | Task Name | Description | Member Assigned |
|-----|-----------|-------------|-----------------|
| **1** | **Dataset Acquisition and Organisation** | Includes:<br><br>• Searching NCBI GEO for all required datasets (GSE137684, GSE114419, GSE98421, GSE277906)<br>• Downloading raw files (CEL files, Series Matrix, raw counts)<br>• Downloading sample metadata (age, BMI, phenotype, tissue)<br>• Organising datasets into labelled folders<br>• Preparing a master spreadsheet of samples, platforms, conditions<br>• Checking for missing samples, mislabeled groups, or platform differences<br>• Saving everything in a version-controlled folder (GitHub/Drive) | **Sri Keerthi Kotamaraju, Bhavini Babariya** |
| **2** | **Microarray Preprocessing (All 3 Microarray Datasets)** | Includes:<br><br>• Background correction (RMA)<br>• Quantile normalization<br>• Log2 transformation<br>• Probe-to-gene mapping using annotation packages<br>• Removing probes without mapping<br>• Averaging multi-probe genes<br>• QC: boxplots, density curves<br>• Pre vs post normalization QC<br>• Creating final expression matrix | **Sri Keerthi Kotamaraju, Bhavini Babariya** |
| **3** | **RNA-seq Preprocessing (GSE277906)** | Includes:<br><br>• Importing raw counts using DESeq2<br>• Filtering low-expression genes<br>• Calculating size factors<br>• VST transformation<br>• QC: mean–variance, sample distances<br>• Normalized matrix creation<br>• Checking PCOS vs control separation<br>• Metadata preparation | **Fanae Ahmed** |
| **4** | **PCA, Heatmaps & Sample Quality Assessment** | Includes:<br><br>• PCA for each dataset<br>• PCOS vs control separation check<br>• Identifying outliers<br>• Heatmaps of variable genes<br>• Hierarchical clustering<br>• Checking batch effects<br>• Recording QC summaries | **Deiby Cabuyales, Mbah Chinedu** |
| **5** | **Differential Gene Expression – Microarray (limma)** | Includes:<br><br>• Designing model matrices<br>• limma linear modeling<br>• Empirical Bayes moderation<br>• Dataset-specific log2FC thresholds<br>• FDR < 0.05 filtering<br>• Volcano plots<br>• DEG tables<br>• Ensuring consistency across datasets | **Youssra Azaf, Mbah Chinedu** |
| **6** | **Differential Gene Expression – RNA-seq (DESeq2)** | Includes:<br><br>• Running DESeq pipeline<br>• Negative binomial modeling<br>• Log2FC shrinkage (optional)<br>• FDR < 0.05<br>• Volcano + MA plots<br>• Saving results as CSV<br>• Up/downregulated gene summary | **Youssra Azaf, Fanae Ahmed** |
| **7** | **Cortisol Gene Set Integration with DEG Lists** | Includes:<br><br>• Importing curated cortisol gene list (NR3C1, FKBP5, DUSP1, etc.)<br>• Importing MSigDB Hallmark Glucocorticoid Response set<br>• Intersecting with DEGs<br>• Venn/UpSet plots<br>• Identifying consistent cortisol-linked DEGs<br>• Stress/metabolic gene tracking | **Laiba Ishtiaq, Youssra Boumait** |
| **8** | **GO Biological Process Enrichment Analysis** | Includes:<br><br>• Running clusterProfiler GO BP enrichment<br>• Filtering adj. p < 0.05<br>• Focus on inflammation, hormone biosynthesis, oxidative stress, cytokine signaling<br>• Dot plots & bar plots<br>• Writing term interpretations | **Laiba Ishtiaq, Maureen Nwosu** |
| **9** | **KEGG Pathway Enrichment Analysis** | Includes:<br><br>• Running KEGG enrichment<br>• Pathways: steroidogenesis, insulin resistance, metabolic dysregulation, inflammatory & glucocorticoid signaling<br>• Visualizations<br>• Ranking pathways linked to PCOS biology | **Kanwal Naz** |
| **10** | **STRING Protein–Protein Interaction Network** | Includes:<br><br>• Importing DEGs/cortisol-DEGs to STRING<br>• Setting confidence > 0.7<br>• Exporting TSV networks<br>• Identifying subnetworks<br>• Preparing PPI for Cytoscape | **Emily Dorado** |
| **11** | **Cytoscape Module Analysis (MCODE + CytoHubba)** | Includes:<br><br>• Importing STRING network<br>• Running MCODE clustering<br>• Hub gene identification (Degree + MCC)<br>• Annotating modules<br>• Creating publication-quality figures<br>• Selecting final hub genes | **Emily, Bihar and Vanisha** |
| **12** | **WGCNA Co-expression Network Analysis** | Includes:<br><br>• Soft-threshold selection<br>• Signed network construction<br>• Module detection<br>• Module eigengenes<br>• Correlation with PCOS & cortisol ssGSEA<br>• Identifying cortisol-linked modules<br>• Extracting hub genes<br>• Writing module summaries | **—** |
| **13** | **Integrated Biological Interpretation** | Includes:<br><br>• Combining DEGs + cortisol genes + enrichment + PPI + WGCNA<br>• Identifying consistent pathways<br>• Highlighting stress-related regulatory motifs<br>• Linking to PCOS inflammation, metabolism, steroidogenesis<br>• Writing final narrative | **Vanisha Sharma, Manar Elabd, Maureen Nwosu** |
| **14** | **Final Report, Figures, Proofreading** | Includes:<br><br>• Writing Introduction, Methods, Results, Discussion<br>• Figure captions & uniform styles<br>• Scientific accuracy check<br>• Reference formatting<br>• Spelling/grammar review<br>• Creating final PDF/presentation | **Vanisha Sharma, Manar Elabd** |

---

