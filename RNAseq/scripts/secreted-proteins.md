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

membrane.genes <- getBM(attributes=c('ensembl_gene_id', 'go_id'),
                        filters = 'go', 
                        values = 'GO:0016020',
                        mart = mouse.data) %>%
  pull(ensembl_gene_id)


secreted.genes <-
  rnaseq.stats %>%
  mutate(SignalP=Row.names %in% signalp.genes,
         GO_Membrane=Row.names %in% membrane.genes) %>%
  mutate(Secreted=if_else(SignalP==TRUE&GO_Membrane==FALSE, 'Secreted','Non-Secreted'))
  
secreted.genes %>%  
  count(Secreted) %>%
  kable(caption="Predicted secreted gene products")
```



Table: Predicted secreted gene products

Secreted            n
-------------  ------
Non-Secreted    19610
Secreted         1780

```r
rnaseq.secreted <-
  rnaseq.stats %>%
  filter(Row.names %in% signalp.genes)  %>% #keep proteins with signal peptide
  filter(!(Row.names %in% membrane.genes))  #remove membrane proteins

secreted.output.file <- '../data/processed/Potential Secreted Genes.csv'

rnaseq.secreted %>%
  rename('ENSEMBL Gene ID'='Row.names',
         'Gene Name'='external_gene_name') %>%
  dplyr::select(-X1) %>%
  dplyr::select('Gene Name', everything()) %>%
  arrange(-abs(log2FoldChange)) %>%
  write_csv(secreted.output.file)
```

These data can be found in **/Users/davebrid/Documents/GitHub/TissueSpecificTscKnockouts/RNAseq/scripts**.  The normalized RNAseq data can be found in a file named **../data/processed/Binary DESeq Results.csv**.  This script was most recently updated on **Thu Mar 26 20:03:20 2020**.

# Analysis

The ENSEMBL dataset with genes annotated as having a signalP annotation includes 3766 and 157001 that are filtered with the GO annotation of 0016020, cellular component - membrane.  This resulted in 1780 secreted genes.

Among 4403 genes that were differentially expressed in *Tsc1* knockout muscles, 253 were potential secreted proteins.

These differentially expressed secreted proteins can be found in ../data/processed/Potential Secreted Genes.csv. 


```r
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

    X1  Row.names             baseMean   log2FoldChange   lfcSE    stat   pvalue    padj  external_gene_name         FC
------  -------------------  ---------  ---------------  ------  ------  -------  ------  -------------------  --------
  7784  ENSMUSG00000030483        7.32             7.24   0.679   10.67    0.000   0.000  Cyp2b10               151.080
 11042  ENSMUSG00000038508        3.23             5.53   0.750    7.38    0.000   0.000  Gdf15                  46.343
   200  ENSMUSG00000001131       13.98             4.24   0.445    9.53    0.000   0.000  Timp1                  18.925
 10293  ENSMUSG00000036564       78.11             3.98   0.549    7.25    0.000   0.000  Ndrg4                  15.790
 14364  ENSMUSG00000050069       37.19            -3.90   0.892   -4.38    0.000   0.000  Grem2                   0.067
 16151  ENSMUSG00000059201       18.23            -3.82   1.349   -2.83    0.005   0.021  Lep                     0.071
  7947  ENSMUSG00000030772      115.02            -3.69   0.446   -8.27    0.000   0.000  Dkk3                    0.077
 14285  ENSMUSG00000049723        5.99             3.58   0.438    8.15    0.000   0.000  Mmp12                  11.914
  8783  ENSMUSG00000032289       21.70             3.32   0.256   12.94    0.000   0.000  Thsd4                   9.966
  4002  ENSMUSG00000023031        2.91             3.12   0.722    4.33    0.000   0.000  Cela1                   8.710

```r
gdf15.data <- filter(rnaseq.secreted, external_gene_name=="Gdf15")  
ggplot(rnaseq.secreted,
       aes(x=log2FoldChange,
           y=-log10(padj))) +
  geom_point(aes(col=padj>0.05)) +
  labs(y="P-value (-log10)",
       x="mRNA Fold Change (log2)",
       title="Volcano Plot of Potentially Secreted Genes") +
  geom_segment(
     xend = gdf15.data$log2FoldChange-0.05, yend = -log10(gdf15.data$padj)+1,
     x=5, y=25,
     arrow=arrow(length=unit(0.1,'cm'))) +
  geom_text(label="Gdf15", x=5, y=27) +
  scale_color_grey() +
  theme_classic() +
  theme(text=element_text(size=18),
        legend.position='none')
```

![](figure/potential-secreted-proteins-1.png)<!-- -->


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
##  [1] Rcpp_1.0.4           lubridate_1.7.4      prettyunits_1.1.1   
##  [4] assertthat_0.2.1     digest_0.6.25        utf8_1.1.4          
##  [7] BiocFileCache_1.10.2 R6_2.4.1             plyr_1.8.6          
## [10] stats4_3.6.2         RSQLite_2.2.0        evaluate_0.14       
## [13] httr_1.4.1           highr_0.8            pillar_1.4.3        
## [16] rlang_0.4.5          progress_1.2.2       curl_4.3            
## [19] blob_1.2.1           S4Vectors_0.24.3     rmarkdown_2.1       
## [22] labeling_0.3         RefManageR_1.2.12    stringr_1.4.0       
## [25] bit_1.1-15.2         munsell_0.5.0        compiler_3.6.2      
## [28] xfun_0.12            pkgconfig_2.0.3      askpass_1.1         
## [31] BiocGenerics_0.32.0  htmltools_0.4.0      openssl_1.4.1       
## [34] tidyselect_1.0.0     tibble_2.1.3         IRanges_2.20.2      
## [37] XML_3.99-0.3         fansi_0.4.1          crayon_1.3.4        
## [40] dbplyr_1.4.2         withr_2.1.2          rappdirs_0.3.1      
## [43] grid_3.6.2           jsonlite_1.6.1       gtable_0.3.0        
## [46] lifecycle_0.2.0      DBI_1.1.0            magrittr_1.5        
## [49] scales_1.1.0         bibtex_0.4.2.2       cli_2.0.2           
## [52] stringi_1.4.6        farver_2.0.3         xml2_1.2.5          
## [55] vctrs_0.2.4          tools_3.6.2          bit64_0.9-7         
## [58] Biobase_2.46.0       glue_1.3.2           purrr_0.3.3         
## [61] hms_0.5.3            parallel_3.6.2       yaml_2.2.1          
## [64] AnnotationDbi_1.48.0 colorspace_1.4-1     memoise_1.1.0
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
