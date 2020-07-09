---
title: "Analysis of AMPK Knockout Data at Sacrifice"
author: "Katherine Kistler, Cody Cousineau, JeAnna Redd and Dave Bridges"
date: "July 5, 2019"
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

# Experimental Details

Animals were sacrificed at approximately 2PM in the fed state

# Raw Data



These data can be found in **/Users/davebrid/Documents/GitHub/TissueSpecificTscKnockouts/Mouse Data/Liver AMPK Ketogenic Diet/All Figures/Sacrifice Data** in a file named **AMPK KD Sacrifice Data.csv**.  This script was most recently updated on **Thu Jul  9 13:12:41 2020**.

# Analysis

# Sacrificed Animals

This is for animals where we have any sacrifice data.


Table: Animals in each group of this cohort

Sex   Diet      Injection     n
----  --------  ----------  ---
F     Control   Cre           7
F     Control   GFP           8
F     Keto      Cre           8
F     Keto      GFP           7
M     Control   Cre           7
M     Control   GFP          10
M     Keto      Cre           4
M     Keto      GFP          16

## Serum Levels

Serum was collected via retro-orbital bleed in the fed state.  Both glucose and ketone body levels were assessed by stick assays.

### Fed Blood Glucose

![Fed glucose levels](figures/glucose-boxplot-1.png)

![Fed glucose levels](figures/glucose-barplot-1.png)


Table: Fed glucose levels

term            estimate   std.error   statistic   p.value
-------------  ---------  ----------  ----------  --------
(Intercept)       177.74        8.88      20.007     0.000
SexM               -1.47        8.76      -0.167     0.868
DietKeto           -2.79        8.62      -0.323     0.748
InjectionGFP       -5.05        8.98      -0.562     0.576

### Fed Ketone Bodies

![Ketone body levels.](figures/ketone-boxplot-1.png)

![Ketone body levels](figures/ketone-barplot-1.png)


Table: Inguinal adipose tissue weights

term            estimate   std.error   statistic   p.value
-------------  ---------  ----------  ----------  --------
(Intercept)        2.748       0.598       4.593     0.000
SexM              -1.707       0.531      -3.214     0.005
DietKeto          -0.564       0.547      -1.031     0.317
InjectionGFP      -0.179       0.467      -0.383     0.707

## Fat Pad Weights

### Inguinal Adipose Tissue

![Inguinal adipose tissue mass](figures/iwat-boxplot-1.png)

![Inguinal adipose tissue mass](figures/iwat-barplot-1.png)


Table: Inguinal adipose tissue weights

term            estimate   std.error   statistic   p.value
-------------  ---------  ----------  ----------  --------
(Intercept)       173.76        19.3       9.006     0.000
SexM               10.03        19.7       0.509     0.613
DietKeto           39.82        19.1       2.079     0.042
InjectionGFP       -4.06        20.1      -0.202     0.840

### Gonadal Adipose Tissue

![Gonadal adipose tissue mass](figures/gwat-boxplot-1.png)

![Gonadal adipose tissue mass](figures/gwat-barplot-1.png)


Table: Gonadal adipose tissue weights

term            estimate   std.error   statistic   p.value
-------------  ---------  ----------  ----------  --------
(Intercept)       316.61        33.8       9.361     0.000
SexM               -9.38        34.1      -0.275     0.784
DietKeto          115.48        33.4       3.462     0.001
InjectionGFP      -10.96        34.9      -0.314     0.755


## Muscle Weights

### Gastrocnemius

![Inguinal adipose tissue mass](figures/gastroc-boxplot-1.png)

![Gastroc tissue mass](figures/gastroc-barplot-1.png)


Table: Inguinal adipose tissue weights

term            estimate   std.error   statistic   p.value
-------------  ---------  ----------  ----------  --------
(Intercept)       101.68        5.83      17.431     0.000
SexM               28.35        5.88       4.822     0.000
DietKeto           -3.61        5.75      -0.627     0.533
InjectionGFP       10.89        6.02       1.809     0.075

### Quadriceps

![Quadriceps mass](figures/quad-boxplot-1.png)

![Quadriceps mass](figures/quad-barplot-1.png)


Table: Quadricep weights

term            estimate   std.error   statistic   p.value
-------------  ---------  ----------  ----------  --------
(Intercept)      134.078        6.46      20.765     0.000
SexM              29.419        6.51       4.521     0.000
DietKeto           0.521        6.37       0.082     0.935
InjectionGFP       7.876        6.67       1.182     0.242



# Session Information


```r
sessionInfo()
```

```
## R version 4.0.0 (2020-04-24)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Catalina 10.15.5
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
## [1] broom_0.5.6     ggplot2_3.3.0   forcats_0.5.0   lubridate_1.7.8
## [5] readr_1.3.1     dplyr_0.8.5     tidyr_1.0.3     knitr_1.28     
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.4.6     pillar_1.4.4     compiler_4.0.0   highr_0.8       
##  [5] tools_4.0.0      digest_0.6.25    lattice_0.20-41  nlme_3.1-147    
##  [9] evaluate_0.14    lifecycle_0.2.0  tibble_3.0.1     gtable_0.3.0    
## [13] pkgconfig_2.0.3  rlang_0.4.6      yaml_2.2.1       xfun_0.13       
## [17] withr_2.2.0      stringr_1.4.0    generics_0.0.2   vctrs_0.2.4     
## [21] hms_0.5.3        grid_4.0.0       tidyselect_1.0.0 glue_1.4.0      
## [25] R6_2.4.1         rmarkdown_2.1    purrr_0.3.4      farver_2.0.3    
## [29] magrittr_1.5     backports_1.1.6  scales_1.1.0     ellipsis_0.3.0  
## [33] htmltools_0.4.0  assertthat_0.2.1 colorspace_1.4-1 labeling_0.3    
## [37] stringi_1.4.6    munsell_0.5.0    crayon_1.3.4
```
