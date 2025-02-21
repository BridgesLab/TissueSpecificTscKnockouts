---
title: "Preparation of aTSC Mammary Gland datasets for GSEA analyses"
author: "Dave Bridges"
date: "September 25, 2020"
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

To generate files for GSEA and gene set enrichment analyses

# Raw Data

Imported DESeq analysed data

These data can be found in **/Users/davebrid/Documents/GitHub/TissueSpecificTscKnockouts/RNAseq/Mammary Gland Adipocyte Tsc1 Knockout**.  

# Analysis

# Data Entry


``` r
library(readr)
deseq.results.file <-  'DESeq2 Results.csv'
deseq.results <- read_csv(deseq.results.file) %>%
  filter(baseMean>100) # filter for base counts have to be above 25
```

# GSEA Prerank Input

Needs human gene identifiers, re-arranged in order by fold change, output into a tsv file.


``` r
human.mouse.mapping.table <- 'http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt'
mapping.table <- read_tsv(human.mouse.mapping.table) %>%
  select(`Common Organism Name`,Symbol,`DB Class Key`) 
  
wide.mapping.table <- pivot_wider(mapping.table,
                                  names_from=`Common Organism Name`,
                                  values_from = Symbol,
                                  id_cols=`DB Class Key`) %>%
  rename("Mouse"='mouse, laboratory') %>%
  mutate(Mouse=as.factor(as.character(Mouse))) %>%
    mutate(human=as.factor(as.character(human)))

mapped.data <-
  deseq.results %>%
  left_join(wide.mapping.table, by=c('symbol'='Mouse')) %>%
  filter(human !='NULL')

output.file <- 'GSEA Ranked File - Effects of aTSC Knockout.rnk'
mapped.data %>%
  arrange(-log2FoldChange) %>%
  select(human, log2FoldChange) %>%
  filter(!(is.na(log2FoldChange))) %>%
  distinct(human, .keep_all=T) %>%
  write_tsv(output.file, col_names = F)

output.file <- 'GSEA Ranked File - Effects of aTSC Knockout - Mouse.rnk'
mapped.data %>%
  arrange(-log2FoldChange) %>%
  select(symbol, log2FoldChange) %>%
  filter(!(is.na(log2FoldChange))) %>%
  distinct(symbol, .keep_all=T) %>%
  write_tsv(output.file, col_names = F)

output.file.exp <- 'GSEA Ranked File - Effects of aTSC Knockout expressed.rnk'
mapped.data %>%
  arrange(-log2FoldChange) %>%
  filter(baseMean>100) %>%
  select(human, log2FoldChange) %>%
  filter(!(is.na(log2FoldChange))) %>%
  distinct(human, .keep_all=T) %>%
  write_tsv(output.file.exp, col_names = F)
```

Used the Jax human mouse orthology tables at http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt

Ran through GSEA 4.1.0 against MSigDB 7.2

# Erichr files

For this just need a list of differentially expressed up and downregulated genes


``` r
deseq.results %>% filter(padj<0.05) %>% pull(symbol) %>% write('Differentially expressed genes.txt')
deseq.results %>% filter(padj<0.05,log2FoldChange>0) %>% pull(symbol) %>% write('Differentially expressed upregulated genes.txt')
deseq.results %>% filter(padj<0.05,log2FoldChange<0) %>% pull(symbol) %>% write('Differentially expressed downregulated genes.txt')
```


# Session Information


``` r
sessionInfo()
```

```
## R version 4.4.2 (2024-10-31)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sequoia 15.1.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/Detroit
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] readr_2.1.5 dplyr_1.1.4 tidyr_1.3.1 knitr_1.49 
## 
## loaded via a namespace (and not attached):
##  [1] bit_4.5.0         jsonlite_1.8.9    compiler_4.4.2    crayon_1.5.3     
##  [5] tidyselect_1.2.1  parallel_4.4.2    jquerylib_0.1.4   yaml_2.3.10      
##  [9] fastmap_1.2.0     R6_2.5.1          generics_0.1.3    curl_6.0.1       
## [13] tibble_3.2.1      bslib_0.8.0       pillar_1.9.0      tzdb_0.4.0       
## [17] rlang_1.1.4       utf8_1.2.4        cachem_1.1.0      xfun_0.49        
## [21] sass_0.4.9        bit64_4.5.2       cli_3.6.3         withr_3.0.2      
## [25] magrittr_2.0.3    digest_0.6.37     vroom_1.6.5       hms_1.1.3        
## [29] lifecycle_1.0.4   vctrs_0.6.5       evaluate_1.0.1    glue_1.8.0       
## [33] fansi_1.0.6       rmarkdown_2.29    purrr_1.0.2       tools_4.4.2      
## [37] pkgconfig_2.0.3   htmltools_0.5.8.1
```
