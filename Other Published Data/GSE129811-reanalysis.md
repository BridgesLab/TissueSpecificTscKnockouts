---
title: "Re-Analysis of GSE129811"
author: "Dave Bridges"
date: "March 19, 2021"
output:
  html_document:
    keep_md: yes
  pdf_document:
    keep_tex: yes
---



This is based on the Forsstrom *et al.* paper which did RNAseq on muscle biopsies from patients with mitochondrial myopathy https://doi.org/10.1016/j.cmet.2019.08.019.  This script was most recently run on Wed Jan 12 12:28:45 2022.


```r
library(readr)
counts.file <- "GSE129811_AWVM29samples.HTSeq.counts.tsv.gz"
counts.data <- read_tsv(counts.file)

patient.file <- "GSE129811_Readme_PatientControlMatching.txt"
patient.data <- read_table(patient.file)

#had to manually assemble from https://www.ncbi.nlm.nih.gov/gds/?term=GSE129811[ACCN]%20AND%20gsm[ETYP] used only time 0
mapping.table <- data.frame(Sample='AWVM23_P0',Patient='Patient3',Time=0) %>%
  add_row(Sample='AWVM20_P0',Patient='Patient4',Time=0) %>% 
  add_row(Sample='AWVM26_P0',Patient='Patient5',Time=0) %>%
  add_row(Sample='AWVM17_P0',Patient='Patient2',Time=0) %>%
  add_row(Sample='AWVM13_C0',Patient='Control9',Time=0) %>%
  add_row(Sample='AWVM15_C0',Patient='Control8',Time=0) %>%
  add_row(Sample='AWVM09_C0',Patient='Control7',Time=0) %>%
  add_row(Sample='AWVM11_C0',Patient='Control6',Time=0) %>%
  add_row(Sample='AWVM03_C0',Patient='Control5',Time=0) %>%
  add_row(Sample='AWVM07_C0',Patient='Control4',Time=0) %>%
  add_row(Sample='AWVM05_C0',Patient='Control3',Time=0) %>%
  add_row(Sample='AWVM01_C0',Patient='Control1',Time=0) %>%
  mutate(Group=case_when(grepl('^Patient',Patient) ~ 'Myopathy',
                         grepl('^Control',Patient) ~ 'Control'))
```

Reanalyzed from counts file at GSE129811_AWVM29samples.HTSeq.counts.tsv.gz.

Patient mapping is below


```r
mapping.table %>% kable
```



|Sample    |Patient  | Time|Group    |
|:---------|:--------|----:|:--------|
|AWVM23_P0 |Patient3 |    0|Myopathy |
|AWVM20_P0 |Patient4 |    0|Myopathy |
|AWVM26_P0 |Patient5 |    0|Myopathy |
|AWVM17_P0 |Patient2 |    0|Myopathy |
|AWVM13_C0 |Control9 |    0|Control  |
|AWVM15_C0 |Control8 |    0|Control  |
|AWVM09_C0 |Control7 |    0|Control  |
|AWVM11_C0 |Control6 |    0|Control  |
|AWVM03_C0 |Control5 |    0|Control  |
|AWVM07_C0 |Control4 |    0|Control  |
|AWVM05_C0 |Control3 |    0|Control  |
|AWVM01_C0 |Control1 |    0|Control  |

# DESeq Analysis

Used bulk, unpaired analysisof only the 14 time zero samples


```r
library(DESeq2)
counts.data.used <- 
  counts.data %>%
  select(GeneID,contains(mapping.table$Sample)) %>%
  tibble::column_to_rownames('GeneID')

dds <- DESeqDataSetFromMatrix(countData = counts.data.used,
                              colData = mapping.table,
                              design = ~ Group)

dds$Sample <- relevel(dds$Group, ref="Control")
dds <- DESeq(dds)
res <- results(dds)
```

# Analysis


```r
plotCounts(dds, gene='ENSG00000170290', intgroup="Group",main="SLN") #SLN
```

![](figures/gse129811-analysis-plots-1.png)<!-- -->

```r
plotCounts(dds, gene='ENSG00000128272', intgroup="Group",main="ATF4") #ATF4
```

![](figures/gse129811-analysis-plots-2.png)<!-- -->

```r
plotCounts(dds, gene='ENSG00000130513', intgroup="Group",main="GDF15") #GDF15
```

![](figures/gse129811-analysis-plots-3.png)<!-- -->

```r
plotCounts(dds, gene='ENSG00000116717', intgroup="Group",main="GADD45A") #GADD45A
```

![](figures/gse129811-analysis-plots-4.png)<!-- -->

```r
sessionInfo()
```

```
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS  10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] DESeq2_1.30.1               SummarizedExperiment_1.20.0
##  [3] Biobase_2.50.0              MatrixGenerics_1.2.1       
##  [5] matrixStats_0.61.0          GenomicRanges_1.42.0       
##  [7] GenomeInfoDb_1.26.7         IRanges_2.24.1             
##  [9] S4Vectors_0.28.1            BiocGenerics_0.36.1        
## [11] readr_2.1.1                 broom_0.7.11               
## [13] dplyr_1.0.7                 tidyr_1.1.4                
## [15] knitr_1.37                 
## 
## loaded via a namespace (and not attached):
##  [1] httr_1.4.2             sass_0.4.0             bit64_4.0.5           
##  [4] vroom_1.5.7            jsonlite_1.7.2         splines_4.0.2         
##  [7] bslib_0.3.1            assertthat_0.2.1       highr_0.9             
## [10] blob_1.2.2             GenomeInfoDbData_1.2.4 yaml_2.2.1            
## [13] pillar_1.6.4           RSQLite_2.2.9          backports_1.4.1       
## [16] lattice_0.20-45        glue_1.6.0             digest_0.6.29         
## [19] RColorBrewer_1.1-2     XVector_0.30.0         colorspace_2.0-2      
## [22] htmltools_0.5.2        Matrix_1.4-0           XML_3.99-0.8          
## [25] pkgconfig_2.0.3        genefilter_1.72.1      zlibbioc_1.36.0       
## [28] purrr_0.3.4            xtable_1.8-4           scales_1.1.1          
## [31] tzdb_0.2.0             BiocParallel_1.24.1    tibble_3.1.6          
## [34] annotate_1.68.0        ggplot2_3.3.5          generics_0.1.1        
## [37] ellipsis_0.3.2         cachem_1.0.6           cli_3.1.0             
## [40] survival_3.2-13        magrittr_2.0.1         crayon_1.4.2          
## [43] memoise_2.0.1          evaluate_0.14          fansi_1.0.0           
## [46] tools_4.0.2            hms_1.1.1              lifecycle_1.0.1       
## [49] stringr_1.4.0          locfit_1.5-9.4         munsell_0.5.0         
## [52] DelayedArray_0.16.3    AnnotationDbi_1.52.0   compiler_4.0.2        
## [55] jquerylib_0.1.4        rlang_0.4.12           grid_4.0.2            
## [58] RCurl_1.98-1.5         rstudioapi_0.13        bitops_1.0-7          
## [61] rmarkdown_2.11         gtable_0.3.0           DBI_1.1.2             
## [64] R6_2.5.1               fastmap_1.1.0          bit_4.0.4             
## [67] utf8_1.2.2             stringi_1.7.6          Rcpp_1.0.7            
## [70] vctrs_0.3.8            geneplotter_1.68.0     tidyselect_1.1.1      
## [73] xfun_0.29
```
