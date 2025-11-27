# GSE137684
Console GSE137684: 
> #Metadata table
> sample_info <- data.frame(Sample = c("GSM4084711", "GSM4084712", "GSM4084713", 
+                                      "GSM4084714", "GSM4084715", "GSM4084716",
+                                      "GSM4084717", "GSM4084718", "GSM4084719",
+                                      "GSM4084720", "GSM4084721", "GSM4084722"),
+                           Condition = c("pcos", "normal", 
+                                         "pcos", "pcos",
+                                         "pcos", "normal", 
+                                         "pcos", "normal", "normal",
+                                         "pcos", "pcos",
+                                         "pcos"),
+                           Title = c("Normoandrogenic PCOS s1", "Normal s3", 
+                                     "Hyperandrogenic PCOS s5", "Hyperandrogenic PCOS s9",
+                                     "Hyperandrogenic PCOS s10", "Normal s12",
+                                     "Hyperandrogenic PCOS s13", "	Normal s14", "Normal s15",
+                                     "Normoandrogenic PCOS s16", "Normoandrogenic PCOS s17",
+                                     "Normoandrogenic PCOS s18"))
> 
> print("1. METADATA TABLE")
[1] "1. METADATA TABLE"
> print(sample_info)
       Sample Condition                    Title
1  GSM4084711      pcos  Normoandrogenic PCOS s1
2  GSM4084712    normal                Normal s3
3  GSM4084713      pcos  Hyperandrogenic PCOS s5
4  GSM4084714      pcos  Hyperandrogenic PCOS s9
5  GSM4084715      pcos Hyperandrogenic PCOS s10
6  GSM4084716    normal               Normal s12
7  GSM4084717      pcos Hyperandrogenic PCOS s13
8  GSM4084718    normal             \tNormal s14
9  GSM4084719    normal               Normal s15
10 GSM4084720      pcos Normoandrogenic PCOS s16
11 GSM4084721      pcos Normoandrogenic PCOS s17
12 GSM4084722      pcos Normoandrogenic PCOS s18
> #Model matrix design(limma)
> group <- factor(sample_info$Condition, levels = c("normal", "pcos"))
> design <- model.matrix(~ 0 + group)
> colnames(design) <- c("Normal","PCOS")
> 
> print("2. MODEL MATRIX DESIGN")
[1] "2. MODEL MATRIX DESIGN"
> print(design)
   Normal PCOS
1       0    1
2       1    0
3       0    1
4       0    1
5       0    1
6       1    0
7       0    1
8       1    0
9       1    0
10      0    1
11      0    1
12      0    1
attr(,"assign")
[1] 1 1
attr(,"contrasts")
attr(,"contrasts")$group
[1] "contr.treatment"

> #Contrast matrix
> contrast.matrix <- makeContrasts(PCOSvsNormal = PCOS - Normal,
+                                  levels = design)
> print("3. CONTRAST MATRIX")
[1] "3. CONTRAST MATRIX"
> print(contrast.matrix)
        Contrasts
Levels   PCOSvsNormal
  Normal           -1
  PCOS              1
> #Fit limma model+topTable
> fit <- lmFit(expr_gene, design)
> fit2 <- contrasts.fit(fit, contrast.matrix)
> fit2 <- eBayes(fit2)
> 
> deg_GSE137684 <- topTable(fit2, coef = "PCOSvsNormal", number = Inf)
> 
> print("4.topTable OUTPUT (HEAD):")
[1] "4.topTable OUTPUT (HEAD):"
> head(deg_GSE137684)
              logFC   AveExpr         t      P.Value adj.P.Val         B
DEFB105B  -1.407618  1.965824 -6.716869 2.983911e-05 0.7786155 -3.432272
MCF2L-AS1  2.826067  5.288644  6.365741 4.854211e-05 0.7786155 -3.464173
FLJ30679  -1.012459  1.837790 -5.621598 1.438577e-04 0.9393607 -3.544920
JDP2       1.311107 11.005755  5.481546 1.780240e-04 0.9393607 -3.562427
PDC       -1.189534  1.987993 -5.467368 1.819344e-04 0.9393607 -3.564244
DMRTC1     1.911557  4.118429  5.325109 2.266088e-04 0.9393607 -3.582957
> 
> write.csv(deg_GSE137684,"GSE137684_DEGs.csv")

Subfolder for 02_Microarray_Preprocessing.
By Bhavini Babariya
