# GSE98421
Console GSE98421: 
> #Metadata table
> sample_info <- data.frame(Sample = c("GSM2595136", "GSM2595137", "GSM2595138", 
+                                      "GSM2595139", "GSM2595140", "GSM2595141", 
+                                      "GSM2595142", "GSM2595143"), 
+                           Condition = c("normal", "normal", "normal", "normal", 
+                                         "PCOS", "PCOS", "PCOS", "PCOS"),
+                           Title = c("normoandrogenic woman, biological rep1", 
+                           "normoandrogenic woman, biological rep2", 
+                           "normoandrogenic woman, biological rep3", 
+                           "normoandrogenic woman, biological rep4", 
+                           "PCOS woman, biological rep1", "PCOS woman, biological rep2", 
+                           "PCOS woman, biological rep3", "PCOS woman, biological rep4"))
> 
> print("1. METADATA TABLE")
[1] "1. METADATA TABLE"
> print(head(sample_info))
      Sample Condition                                  Title
1 GSM2595136    normal normoandrogenic woman, biological rep1
2 GSM2595137    normal normoandrogenic woman, biological rep2
3 GSM2595138    normal normoandrogenic woman, biological rep3
4 GSM2595139    normal normoandrogenic woman, biological rep4
5 GSM2595140      PCOS            PCOS woman, biological rep1
6 GSM2595141      PCOS            PCOS woman, biological rep2
> #Model matrix design
>  group <- factor(sample_info$Condition, levels = c("normal", "PCOS"))
>  design <- model.matrix(~ 0 + group)
>  colnames(design) <- c("Normal", "PCOS")
>  
> print("2. MODEL MATRIX DESIGN") 
[1] "2. MODEL MATRIX DESIGN"
> print(design)
  Normal PCOS
1      1    0
2      1    0
3      1    0
4      1    0
5      0    1
6      0    1
7      0    1
8      0    1
attr(,"assign")
[1] 1 1
attr(,"contrasts")
attr(,"contrasts")$group
[1] "contr.treatment"

> #Contrast matrix
>  contrast.matrix <- makeContrasts(PCOSvsNormal = PCOS - Normal,
+                                   levels = design)
> print("3. CONTRAST MATRIX")
[1] "3. CONTRAST MATRIX"
> print(contrast.matrix)
        Contrasts
Levels   PCOSvsNormal
  Normal           -1
  PCOS              1
> #Fit limma + extract DEGs
>  fit <- lmFit(data, design)
>  fit2 <- contrasts.fit(fit, contrast.matrix)
>  fit2 <- eBayes(fit2)
>  
>  deg_GSE98421 <- topTable(fit2, coef = "PCOSvsNormal", number = Inf)
>  
> print("4. topTable OUTPUT")
[1] "4. topTable OUTPUT"
> head(deg_GSE98421)
               logFC   AveExpr         t      P.Value adj.P.Val         B
PRSS12    -0.8614519  7.877965 -6.492913 0.0001191183 0.9999658 -3.121944
COL11A1   -1.2712399  8.041987 -5.830361 0.0002625899 0.9999658 -3.202086
EDIL3     -1.2948548 10.354183 -5.425899 0.0004377360 0.9999658 -3.260772
PAXIP1-DT  0.6134625  6.226603  5.379479 0.0004648472 0.9999658 -3.268053
CRISPLD2   0.7659415  7.371959  5.286319 0.0005249061 0.9999658 -3.283028
THAP7-AS1  0.4797740  6.153828  4.814706 0.0009896552 0.9999658 -3.366871
> write.csv(deg_GSE98421, "GSE98421_DEG.csv")

Subfolder for 02_Microarray_Preprocessing.
By Bhavini Babariya
