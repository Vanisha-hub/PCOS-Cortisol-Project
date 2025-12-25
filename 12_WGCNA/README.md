# Task 11 PPI

# Protein–Protein Interaction (PPI) Network Analysis

This folder contains **protein–protein interaction (PPI) network analyses** performed using the STRING database for key gene sets identified in the *PCOS–Cortisol Project*. The analyses aim to uncover **functional interaction modules**, **network hubs**, and **biological organization** of cortisol-linked differentially expressed genes across multiple transcriptomic datasets.

##  Folder Contents & File Description

### Dataset-Specific PPI Clusters

#### GSE137684
- **GSE137684 Cluster 1.png**  
- **GSE137684 Cluster 2.png**  
- **GSE137684 Cluster 3.png**  

PPI network clusters identified within the GSE137684 dataset, representing functionally coherent gene modules based on interaction density.

- **GSE137684_Cluster1_schematic.png**  
  Simplified schematic representation of the primary interaction module for biological interpretation.

- **GSE137684_module_family_count.png**  
  Distribution of interacting proteins across gene families within PPI modules.

#### GSE98421
- **GSE98421 Cluster 1.png**  
- **GSE98421 Cluster 2.png**  
- **GSE98421 Cluster 3.png**  

PPI clusters derived from cortisol-associated DEGs in GSE98421.

- **GSE98421_Cluster1_schematic.png**  
  Schematic visualization summarizing dominant interaction patterns.

- **GSE98421_module_family_count.png**  
  Gene family composition of PPI modules identified in this dataset.

---

### Global PPI Overview
- **STRING_overview.png**  
  Global STRING interaction network integrating cortisol-linked genes across datasets, providing a systems-level view of protein connectivity.

---

##  Methodological Overview

- PPI networks were constructed using the **STRING database**
- High-confidence interactions were selected based on:
  - Experimental evidence
  - Curated database interactions
  - Co-expression (where applicable)
- Network clustering was performed to identify **functionally relevant modules**
- Gene family annotation was used to interpret module composition

