---
title: "Cluster Profiler Analysis for aTSC Mammary Glands"
author: "Dave Bridges"
date: "December 21, 2020"
output:
  html_document:
    highlight: tango
    keep_md: yes
    number_sections: yes
    toc: yes
  pdf_document:
    highlight: tango
    keep_tex: yes
    number_sections: yes
    toc: yes
---



# Purpose

To use cluster profiler as part of DOSE to do gene set enrichments

# Raw Data

GSEA was run with folders put in this subfolder



# Gene Ontology




# KEGG 



# Disease Ontology 



# Disease Gene Network

From DisGeNET




# Cell Markers



# Session Information


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
##  [1] DOSE_3.16.0            clusterProfiler_3.18.1 org.Hs.eg.db_3.12.0   
##  [4] AnnotationDbi_1.52.0   IRanges_2.24.1         S4Vectors_0.28.1      
##  [7] Biobase_2.50.0         BiocGenerics_0.36.0    readr_1.4.0           
## [10] dplyr_1.0.5            tidyr_1.1.3            knitr_1.31            
## 
## loaded via a namespace (and not attached):
##  [1] enrichplot_1.10.2   bit64_4.0.5         RColorBrewer_1.1-2 
##  [4] tools_4.0.2         bslib_0.2.4         utf8_1.2.1         
##  [7] R6_2.5.0            DBI_1.1.1           colorspace_2.0-0   
## [10] tidyselect_1.1.0    gridExtra_2.3       curl_4.3           
## [13] bit_4.0.4           compiler_4.0.2      cli_2.3.1          
## [16] scatterpie_0.1.5    shadowtext_0.0.7    sass_0.3.1         
## [19] scales_1.1.1        stringr_1.4.0       digest_0.6.27      
## [22] rmarkdown_2.7       pkgconfig_2.0.3     htmltools_0.5.1.1  
## [25] fastmap_1.1.0       rlang_0.4.10        rstudioapi_0.13    
## [28] RSQLite_2.2.5       jquerylib_0.1.3     generics_0.1.0     
## [31] farver_2.1.0        jsonlite_1.7.2      vroom_1.4.0        
## [34] BiocParallel_1.24.1 GOSemSim_2.16.1     magrittr_2.0.1     
## [37] GO.db_3.12.1        Matrix_1.3-2        Rcpp_1.0.6         
## [40] munsell_0.5.0       fansi_0.4.2         viridis_0.5.1      
## [43] lifecycle_1.0.0     stringi_1.5.3       yaml_2.2.1         
## [46] ggraph_2.0.5        MASS_7.3-53.1       plyr_1.8.6         
## [49] qvalue_2.22.0       grid_4.0.2          blob_1.2.1         
## [52] ggrepel_0.9.1       DO.db_2.9           crayon_1.4.1       
## [55] lattice_0.20-41     graphlayouts_0.7.1  cowplot_1.1.1      
## [58] splines_4.0.2       hms_1.0.0           pillar_1.5.1       
## [61] fgsea_1.16.0        igraph_1.2.6        reshape2_1.4.4     
## [64] fastmatch_1.1-0     glue_1.4.2          evaluate_0.14      
## [67] downloader_0.4      BiocManager_1.30.12 data.table_1.14.0  
## [70] vctrs_0.3.7         tweenr_1.0.2        gtable_0.3.0       
## [73] purrr_0.3.4         polyclip_1.10-0     assertthat_0.2.1   
## [76] cachem_1.0.4        ggplot2_3.3.3       xfun_0.22          
## [79] ggforce_0.3.3       tidygraph_1.2.0     viridisLite_0.3.0  
## [82] tibble_3.1.0        rvcheck_0.1.8       memoise_2.0.0      
## [85] ellipsis_0.3.1
```
