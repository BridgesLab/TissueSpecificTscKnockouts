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


```r
library(readr)
deseq.results.file <-  'DESeq2 Results.csv'
deseq.results <- read_csv(deseq.results.file)
```

# GSEA Prerank Input

Needs human gene identifiers, re-arranged in order by fold change, output into a tsv file.


```r
human.mouse.mapping.table <- 'http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt'
mapping.table <- read_tsv(human.mouse.mapping.table) %>%
  select(`Common Organism Name`,Symbol,`HomoloGene ID`) 
  
wide.mapping.table <- pivot_wider(mapping.table,
                                  names_from=`Common Organism Name`,
                                  values_from = Symbol,
                                  id_cols=`HomoloGene ID`) %>%
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


```r
deseq.results %>% filter(padj<0.05) %>% pull(symbol) %>% write('Differentially expressed genes.txt')
deseq.results %>% filter(padj<0.05,log2FoldChange>0) %>% pull(symbol) %>% write('Differentially expressed upregulated genes.txt')
deseq.results %>% filter(padj<0.05,log2FoldChange<0) %>% pull(symbol) %>% write('Differentially expressed downregulated genes.txt')
```


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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] readr_1.4.0 dplyr_1.0.5 tidyr_1.1.3 knitr_1.31 
## 
## loaded via a namespace (and not attached):
##  [1] rstudioapi_0.13   magrittr_2.0.1    hms_1.0.0         tidyselect_1.1.0 
##  [5] R6_2.5.0          rlang_0.4.10      fansi_0.4.2       stringr_1.4.0    
##  [9] tools_4.0.2       xfun_0.22         utf8_1.2.1        cli_2.3.1        
## [13] DBI_1.1.1         jquerylib_0.1.3   htmltools_0.5.1.1 ellipsis_0.3.1   
## [17] assertthat_0.2.1  yaml_2.2.1        digest_0.6.27     tibble_3.1.0     
## [21] lifecycle_1.0.0   crayon_1.4.1      purrr_0.3.4       sass_0.3.1       
## [25] vctrs_0.3.7       curl_4.3          glue_1.4.2        evaluate_0.14    
## [29] rmarkdown_2.7     stringi_1.5.3     compiler_4.0.2    bslib_0.2.4      
## [33] pillar_1.5.1      generics_0.1.0    jsonlite_1.7.2    pkgconfig_2.0.3
```
