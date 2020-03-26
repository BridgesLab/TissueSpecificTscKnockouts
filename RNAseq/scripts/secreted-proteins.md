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



# Raw Data


```r
library(readr) #loads the readr package
rnaseq.filename <- '../data/processed/Binary DESeq Results.csv' #make this a separate line, you can use any variable you want
rnaseq.analysed.filename <- '../data/processed/Binary Normalized Counts.csv'

#this loads whatever the file is into a dataframe called exp.data if it exists
rnaseq.counts <- read_csv(rnaseq.filename)
rnaseq.stats <- read_csv(rnaseq.analysed.filename)
```

These data can be found in **/Users/davebrid/Documents/GitHub/TissueSpecificTscKnockouts/RNAseq/scripts**.  The normalized RNAseq data can be found in a file named **../data/processed/Binary DESeq Results.csv**.  This script was most recently updated on **Thu Mar 26 16:57:47 2020**.

# Analysis




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
## [1] readr_1.3.1 dplyr_0.8.5 tidyr_1.0.2 knitr_1.28 
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.4       magrittr_1.5     hms_0.5.3        tidyselect_1.0.0
##  [5] R6_2.4.1         rlang_0.4.5      stringr_1.4.0    tools_3.6.2     
##  [9] xfun_0.12        htmltools_0.4.0  yaml_2.2.1       digest_0.6.25   
## [13] assertthat_0.2.1 tibble_2.1.3     lifecycle_0.2.0  crayon_1.3.4    
## [17] purrr_0.3.3      vctrs_0.2.4      glue_1.3.2       evaluate_0.14   
## [21] rmarkdown_2.1    stringi_1.4.6    compiler_3.6.2   pillar_1.4.3    
## [25] pkgconfig_2.0.3
```
