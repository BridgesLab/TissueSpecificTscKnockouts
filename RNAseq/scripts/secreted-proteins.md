---
title: "Analysis of Secreted Proteins from TSC Knockout Muscles"
author: "Dave Bridges"
date: "March 26, 2020"
output:
  html_document:
    highlight: tango
    keep_md: yes
    number_sections: no
    toc: yes
  pdf_document:
    highlight: tango
    keep_tex: yes
    number_sections: yes
    toc: yes
---



# Purpose

To broadly evaluate all potential myokines from mTORC1 activated muscles based on transcriptional changes from our RNAseq data 

# Experimental Details

Used RNAseq data compiled from our previous experiments  Used biomart to extract proteins with annotated signal peptides done via SignalP [@Almagro_Armenteros_2019].

# Raw Data


```r
library(readr) #loads the readr package
rnaseq.filename <- '../data/processed/Binary DESeq Results.csv' #make this a separate line, you can use any variable you want
rnaseq.analysed.filename <- '../data/processed/Binary Normalized Counts.csv'

#this loads whatever the file is into a dataframe called exp.data if it exists
rnaseq.stats <- read_csv(rnaseq.filename)
rnaseq.counts <- read_csv(rnaseq.analysed.filename)
```

## Annotated Proteins with Signal Peptides


```r
library(biomaRt)

mouse.data <- useDataset('mmusculus_gene_ensembl', mart=useMart('ensembl'))
#listFilters(mouse.data) #to locate with_signalp
#listAttributes(mouse.data)  #to locate ensembl_gene_id

signalp.data <- 
  getBM(attributes=c('ensembl_gene_id','signalp'), 
      values = rnaseq.stats$Row.names, 
      mart = mouse.data) %>%
  mutate(signap = as.factor(signalp))

signalp.genes <-
  signalp.data %>%
  filter(signalp %in% c('SignalP-noTM')) %>%
  pull(ensembl_gene_id)
```

These data can be found in **/Users/davebrid/Documents/GitHub/TissueSpecificTscKnockouts/RNAseq/scripts**.  The normalized RNAseq data can be found in a file named **../data/processed/Binary DESeq Results.csv**.  This script was most recently updated on **Thu Mar 26 18:44:25 2020**.

# Analysis

The ENSEMBL dataset with genes annotated as having a signalP annotation includes 3766 or 15.447% of all genes.


```r
rnaseq.secreted <-
  rnaseq.stats %>%
  filter(Row.names %in% signalp.genes) 

library(ggplot2)

sig.secreted <- 
  rnaseq.secreted %>%
  mutate(FC = 2^(log2FoldChange)) %>%
  filter(padj < 0.05,
         baseMean>2) %>%
  arrange(desc(abs(log2FoldChange))) 

sig.secreted %>%
  head(10) %>%
  kable(caption = "Top differentially expressed secreted proteins")
```



Table: Top differentially expressed secreted proteins

    X1  Row.names             baseMean   log2FoldChange   lfcSE    stat   pvalue   padj  external_gene_name       FC
------  -------------------  ---------  ---------------  ------  ------  -------  -----  -------------------  ------
  7784  ENSMUSG00000030483        7.32             7.24   0.679   10.67        0      0  Cyp2b10               151.1
 11042  ENSMUSG00000038508        3.23             5.53   0.750    7.38        0      0  Gdf15                  46.3
 16480  ENSMUSG00000060882        2.91             5.17   0.687    7.54        0      0  Kcnd2                  36.1
  2594  ENSMUSG00000020598       39.14             5.14   0.244   21.09        0      0  Nrcam                  35.2
 14549  ENSMUSG00000050808        3.95             4.98   0.726    6.86        0      0  Muc15                  31.5
   200  ENSMUSG00000001131       13.98             4.24   0.445    9.53        0      0  Timp1                  18.9
  5404  ENSMUSG00000026253        5.21             4.23   0.426    9.93        0      0  Chrng                  18.7
  9326  ENSMUSG00000033676        2.16             4.17   0.614    6.79        0      0  Gabrb3                 18.0
  5375  ENSMUSG00000026204        2.52             4.16   0.664    6.27        0      0  Ptprn                  17.9
 17555  ENSMUSG00000068547        2.30             4.08   0.617    6.62        0      0  Clca4a                 17.0


# Session Information


```r
sessionInfo()
```

```
## R version 3.6.2 (2019-12-12)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Catalina 10.15.3
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggplot2_3.3.0        biomaRt_2.42.0       readr_1.3.1         
## [4] knitcitations_1.0.10 dplyr_0.8.5          tidyr_1.0.2         
## [7] knitr_1.28          
## 
## loaded via a namespace (and not attached):
##  [1] progress_1.2.2       tidyselect_1.0.0     xfun_0.12           
##  [4] purrr_0.3.3          colorspace_1.4-1     vctrs_0.2.4         
##  [7] htmltools_0.4.0      stats4_3.6.2         BiocFileCache_1.10.2
## [10] yaml_2.2.1           blob_1.2.1           XML_3.99-0.3        
## [13] rlang_0.4.5          pillar_1.4.3         withr_2.1.2         
## [16] glue_1.3.2           DBI_1.1.0            rappdirs_0.3.1      
## [19] BiocGenerics_0.32.0  bit64_0.9-7          dbplyr_1.4.2        
## [22] lifecycle_0.2.0      plyr_1.8.6           stringr_1.4.0       
## [25] munsell_0.5.0        gtable_0.3.0         evaluate_0.14       
## [28] memoise_1.1.0        Biobase_2.46.0       IRanges_2.20.2      
## [31] curl_4.3             parallel_3.6.2       AnnotationDbi_1.48.0
## [34] highr_0.8            Rcpp_1.0.4           scales_1.1.0        
## [37] openssl_1.4.1        S4Vectors_0.24.3     jsonlite_1.6.1      
## [40] bit_1.1-15.2         hms_0.5.3            askpass_1.1         
## [43] digest_0.6.25        stringi_1.4.6        grid_3.6.2          
## [46] bibtex_0.4.2.2       tools_3.6.2          magrittr_1.5        
## [49] tibble_2.1.3         RSQLite_2.2.0        RefManageR_1.2.12   
## [52] crayon_1.3.4         pkgconfig_2.0.3      xml2_1.2.5          
## [55] prettyunits_1.1.1    lubridate_1.7.4      assertthat_0.2.1    
## [58] rmarkdown_2.1        httr_1.4.1           R6_2.4.1            
## [61] compiler_3.6.2
```

# Bibliography


```r
write.bibtex(file="secreted-database-references.bib")
bibliography("markdown")
```

<a
name=bib-Almagro_Armenteros_2019></a>[[1]](#cite-Almagro_Armenteros_2019)
J. J. A. Armenteros, K. D. Tsirigos, C. K. Sønderby, et al. "SignalP
5.0 improves signal peptide predictions using deep neural networks".
In: _Nature Biotechnology_ 37.4 (Feb. 2019), pp. 420-423. DOI:
[10.1038/s41587-019-0036-z](https://doi.org/10.1038%2Fs41587-019-0036-z).
URL:
[https://doi.org/10.1038/s41587-019-0036-z](https://doi.org/10.1038/s41587-019-0036-z).
