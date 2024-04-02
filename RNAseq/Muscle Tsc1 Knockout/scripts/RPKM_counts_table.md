Making the RPKM table
========================================================


Making the RPKM Table
----------------------

The counts table was from DESeq2 using **deseq.Rmd**.  



The normalized counts are read in from the ../data/processed/Binary Normalized Counts.csv file.  The transcript database is made from biomart.


```
## Warning: Infering Exon Rankings.  If this is not what you expected, then
## please be sure that you have provided a valid attribute for
## exonRankAttributeName
```

This script writes the RPKM table out to the file ../data/processed/RPKM_count.csv


Session Information
---------------------


```r
sessionInfo()
```

```
## R version 3.1.0 (2014-04-10)
## Platform: x86_64-apple-darwin13.1.0 (64-bit)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] multicore_0.2          GenomicFeatures_1.16.2 AnnotationDbi_1.26.0  
## [4] Biobase_2.24.0         GenomicRanges_1.16.3   GenomeInfoDb_1.0.2    
## [7] IRanges_1.22.8         BiocGenerics_0.10.0    knitr_1.6             
## 
## loaded via a namespace (and not attached):
##  [1] BatchJobs_1.2           BBmisc_1.6             
##  [3] BiocParallel_0.6.1      biomaRt_2.20.0         
##  [5] Biostrings_2.32.0       bitops_1.0-6           
##  [7] brew_1.0-6              BSgenome_1.32.0        
##  [9] codetools_0.2-8         DBI_0.2-7              
## [11] digest_0.6.4            evaluate_0.5.5         
## [13] fail_1.2                foreach_1.4.2          
## [15] formatR_0.10            GenomicAlignments_1.0.1
## [17] iterators_1.0.7         plyr_1.8.1             
## [19] Rcpp_0.11.2             RCurl_1.95-4.1         
## [21] Rsamtools_1.16.0        RSQLite_0.11.4         
## [23] rtracklayer_1.24.2      sendmailR_1.1-2        
## [25] stats4_3.1.0            stringr_0.6.2          
## [27] tools_3.1.0             XML_3.98-1.1           
## [29] XVector_0.4.0           zlibbioc_1.10.0
```

