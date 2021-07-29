DESeq Analysis
========================================================


```r
input_directory <- "../data/raw/htseq-counts-exons"
sample_mapping_file <- "../data/raw/sample_mapping.csv"
annotation_file <- "../data/raw/Mus_musculus.GRCm38.75.gff"
```


This file uses the counts table in ../data/raw/htseq-counts-exons from ../data/raw/sample_mapping.csv.  This script was most recently run on Tue Jun 10 14:27:13 2014


```
## Error: error in evaluating the argument 'x' in selecting a method for function '%in%': Error: argument "groupID" is missing, with no default
```

```
## Error: subscript out of bounds
```

```
## Error: object 'dxd.binary' not found
```

```
## Error: object 'dxd.binary' not found
```

```
## Error: object 'dxd.binary' not found
```

```
## Error: object 'res.binary' not found
```



```
## Error: error in evaluating the argument 'x' in selecting a method for function 'merge': Error in as.data.frame(resOrdered.binary) : 
##   error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': Error: object 'resOrdered.binary' not found
```

```
## Error: object 'annotated_data' not found
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': Error: object 'annotated_data' not found
```

```
## Error: object 'ordered_annotated_data' not found
```

```
## Error: attempt to set 'rownames' on an object with no dimensions
```

```
## Error: attempt to set 'colnames' on an object with less than two
## dimensions
```


The calculated data was exported to ../data/processed/Binary DESeq Results.csv.  There was a total of **

```

Error in eval(expr, envir, enclos) : 
  object 'ordered_annotated_data' not found

```

** significantly different genes out of 

```

Error in eval(expr, envir, enclos) : 
  object 'ordered_annotated_data' not found

```

 that met our testing criteria.  These are broken down into up and downregulated genes below:


```
## Error: no applicable method for 'xtable' applied to an object of class
## "function"
```


Diagnostic Plots
----------------


```
## Error: error in evaluating the argument 'object' in selecting a method for function 'plotMA': Error: object 'res.binary' not found
```



```
## KernSmooth 2.23 loaded
## Copyright M. P. Wand 1997-2009
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:IRanges':
## 
##     space
## 
## The following object is masked from 'package:stats':
## 
##     lowess
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'rownames': Error in subset(as.data.frame(res.binary), padj < 0.05) : 
##   error in evaluating the argument 'x' in selecting a method for function 'subset': Error in as.data.frame(res.binary) : 
##   error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': Error: object 'res.binary' not found
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'counts': Error: object 'dds.binary' not found
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'sizeFactors': Error: object 'dds.binary' not found
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'sizeFactors': Error: object 'dds.binary' not found
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'assay': Error: object 'rld.binary' not found
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'assay': Error: object 'vsd.binary' not found
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'assay': Error: object 'rld.binary' not found
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'assay': Error: object 'vsd.binary' not found
```


Plotting Estimated Dispersions
-------------------------------

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'plotDispEsts': Error: object 'dds.binary' not found
```



Principal Component Analysis
------------------------------


```
## Error: error in evaluating the argument 'x' in selecting a method for function 'assay': Error: object 'rld.binary' not found
```



Session Information
-------------------


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
##  [1] gplots_2.13.0           RColorBrewer_1.0-5     
##  [3] xtable_1.7-3            biomaRt_2.20.0         
##  [5] DEXSeq_1.10.3           BiocParallel_0.6.0     
##  [7] DESeq2_1.4.5            RcppArmadillo_0.4.300.0
##  [9] Rcpp_0.11.1             GenomicRanges_1.16.3   
## [11] GenomeInfoDb_1.0.2      IRanges_1.22.6         
## [13] Biobase_2.24.0          BiocGenerics_0.10.0    
## [15] knitr_1.5              
## 
## loaded via a namespace (and not attached):
##  [1] annotate_1.42.0      AnnotationDbi_1.26.0 BatchJobs_1.2       
##  [4] BBmisc_1.6           Biostrings_2.32.0    bitops_1.0-6        
##  [7] brew_1.0-6           caTools_1.17         codetools_0.2-8     
## [10] DBI_0.2-7            digest_0.6.4         evaluate_0.5.5      
## [13] fail_1.2             foreach_1.4.2        formatR_0.10        
## [16] gdata_2.13.3         genefilter_1.46.1    geneplotter_1.42.0  
## [19] grid_3.1.0           gtools_3.4.0         hwriter_1.3         
## [22] iterators_1.0.7      KernSmooth_2.23-12   lattice_0.20-29     
## [25] locfit_1.5-9.1       plyr_1.8.1           RCurl_1.95-4.1      
## [28] Rsamtools_1.16.0     RSQLite_0.11.4       sendmailR_1.1-2     
## [31] splines_3.1.0        statmod_1.4.19       stats4_3.1.0        
## [34] stringr_0.6.2        survival_2.37-7      tools_3.1.0         
## [37] XML_3.98-1.1         XVector_0.4.0        zlibbioc_1.10.0
```

