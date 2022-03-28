---
title: "Analysis of mTSC1 Knockout CLAMS Experiments"
author: "Erin Stephenson and Dave Bridges"
date: "January 31, 2019"
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

To evaluate energy expenditure and other parameters in muscle _Tsc1_ knockout mice.  This script was most recently updated on **Mon Mar 28 11:57:56 2022**.

# Experimental Details

Mice were run in the CLAMS in several batches, and combined.

# Raw Data

## Sample Key



## Oxymax Input

There are two batches of data, baseline and after 3 months of diet.

### Baseline Data


Table: Total animals tested by genotype

|Genotype    |Sex    |  n|
|:-----------|:------|--:|
|+/+; +/+    |Female | 16|
|+/+; +/+    |Male   | 16|
|+/+; Tg/+   |Female | 11|
|+/+; Tg/+   |Male   |  8|
|fl/fl; +/+  |Female | 17|
|fl/fl; +/+  |Male   |  8|
|fl/fl; Tg/+ |Female | 18|
|fl/fl; Tg/+ |Male   | 11|
|NA          |NA     |  7|



Table: Total animals tested by knockout

|Knockout |Sex    |  n|
|:--------|:------|--:|
|Control  |Female | 44|
|Control  |Male   | 32|
|Knockout |Female | 18|
|Knockout |Male   | 11|
|NA       |NA     |  7|

The baseline raw data files can be found in Oxymax/Oxymax files by time period/Baseline.  The MRI data can be found in EchoMRI.

## VO2 Analysis

![VO2 Summary Light/Dark Boxplot](figures/vo2-analysis-light-dark-1.png)![VO2 Summary Light/Dark Boxplot](figures/vo2-analysis-light-dark-2.png)

![](figures/vo2-analysis-linegraph, fig-1.png)<!-- -->

### VO2 Summary Data

![Linegraph of VO2 Data](figures/vo2-summarized-data-1.png)![Linegraph of VO2 Data](figures/vo2-summarized-data-2.png)

## VCO2 Analysis

![VCO2 Summary Light/Dark Boxplot](figures/vco2-analysis-light-dark-1.png)![VCO2 Summary Light/Dark Boxplot](figures/vco2-analysis-light-dark-2.png)

![](figures/vco2-analysis-linegraph, fig-1.png)<!-- -->

### VCO2 Summary Data
![Linegraph of VCO2 Data](figures/vco2-summarized-data-1.png)![Linegraph of VCO2 Data](figures/vco2-summarized-data-2.png)
# Heat Production

Another way to present these data is to evaluate this by heat instead of VO2. We calculated this manually from VO2 data.  The equation for Heat production from the CLAMS is the Lusk Equation:

$$(3.815 + 1.232 * RER)*VO2$$

![Linegraph of Heat Data](figures/heat-production-1.png)![Linegraph of Heat Data](figures/heat-production-2.png)

## Heat Statistics


Table: Average changes in heat production comparing wt to knockout

|Sex    |Light/Dark | Control| Knockout| Change| Pct.Change|
|:------|:----------|-------:|--------:|------:|----------:|
|Female |Dark       |   0.470|    0.504|  0.034|      7.221|
|Female |Light      |   0.434|    0.433| -0.001|     -0.173|
|Male   |Dark       |   0.474|    0.504|  0.030|      6.395|
|Male   |Light      |   0.437|    0.439|  0.003|      0.581|

To test whether these groups are different we constructed a linear model with the following formula:

Heat ~ as.factor(Zeitgeber.Time) + Lean + Sex + `Light/Dark` + Knockout + Knockout:`Light/Dark` + (1 | Subject).  

We used this model because the base model was that Heat production changes over the day.  We asked if lean mass modified the time dependent effect, and it did (p=` anova(heat.lme.base,heat.lme.lean)$"Pr(>Chisq)"[2]`).  After adjusting for lean mass, we asked if there was any additional benefit to including the light/dark cycle in addition to the time of day, and found that there was no significant effect, so that was not included in the model (p=NA).  we added sex as a covariate which had no significant effect 0.03. We chose to keep sex in the model though as it was borderline significant.  We next added knockout to the model and found no significant effect 0.891.  Finally we asked if Sex modified the effect of the knockout and found no significant effect 0.463.

Since it appears from the figures that the elevation in energy expenditure is restricted to the awake cycle, we next asked if there was an *interaction* between genotype and the Light/Dark cycle.  Adding this interaction was highly significant 3.411&times; 10^-13^.  

The full results are shown below:


Table: Estimates and p-values from mixed linear models, excluding time of day.

|                                   | Estimate| Std..Error| t.value|   p.z|
|:----------------------------------|--------:|----------:|-------:|-----:|
|Lean                               |    0.007|      0.003|   2.057| 0.040|
|SexMale                            |   -0.045|      0.020|  -2.212| 0.027|
|KnockoutKnockout                   |    0.014|      0.014|   0.992| 0.321|
|`Light/Dark`Light:KnockoutKnockout |   -0.032|      0.004|  -7.641| 0.000|

### How would this relate to energy balance?


Table: Average changes in heat production comparing wt to knockout

|Sex    | Control| Knockout| Change| Pct.Change|
|:------|-------:|--------:|------:|----------:|
|Female |   0.452|    0.469|  0.017|       3.67|
|Male   |   0.455|    0.472|  0.016|       3.61|

Based on these calculations, we detected a 16.512mW increase in energy expenditure.  This corresponds to 0.341kcal increase in calories consumed per day.  Over the course of 30 weeks (the NCD study) this accumulates to 71.605kcal which converts to 7.956g of fat mass if there are no other adaptations.  For the HFD studies, this corresponds to a decrease over 11 weeks of 26.255kcal which converts to 2.917g of fat mass.

# RER Analysis

![RER Summary Light/Dark Boxplot](figures/rer-analysis-light-dark-1.png)

![](figures/rer-analysis-linegraph, fig-1.png)<!-- -->

### RER Summary Data

![Linegraph of RER Data](figures/rer-summarized-data-1.png)![Linegraph of RER Data](figures/rer-summarized-data-2.png)

# Carbohydrate Oxidation Analysis

Calculated as $Carbohydrate\ oxidation = (4.585 * vCO_2) - (3.226 * vO_2)$ where both units are in L/min and the output is in g/min

![Carbohydrate Oxidation Summary Light/Dark Boxplot](figures/cho-analysis-light-dark-1.png)

![](figures/cho-analysis-linegraph, fig-1.png)<!-- -->

### Carbohydrate Oxidation Summary Data

![Linegraph of carbohydrate oxidation data](figures/cho-summarized-data-1.png)![Linegraph of carbohydrate oxidation data](figures/cho-summarized-data-2.png)

### Carbohydrate Oxidation Statistics


Table: Average changes in carbohydrate oxidation comparing wt to knockout

|Sex    |Light/Dark | Control| Knockout| Change| Pct.Change|
|:------|:----------|-------:|--------:|------:|----------:|
|Female |Dark       |   0.891|    1.039|  0.148|      16.58|
|Female |Light      |   0.793|    0.746| -0.047|      -5.98|
|Male   |Dark       |   1.056|    1.011| -0.044|      -4.21|
|Male   |Light      |   0.813|    0.759| -0.053|      -6.57|

To test whether these groups are different we constructed a linear model with the following formula:

CHO Oxidation ~ as.factor(Zeitgeber.Time) + Lean + Sex + `Light/Dark` + Knockout + Knockout:`Light/Dark` + (1 | Subject).  

We used this model because the base model was that carbohydrate oxidation changes over the day.  We asked if lean mass modified the time dependent effect, and it did (p=0.058).  After adjusting for lean mass, we asked if there was any additional benefit to including the light/dark cycle in addition to the time of day, and found that there was no significant effect, so that was not included in the model (p=NA).  We added sex as a covariate which had no significant effect 0.945. We chose to keep sex in the model though.  We next added knockout to the model and found no significant effect 0.693.  Finally we asked if Sex modified the effect of the knockout and found no significant effect 0.352.

Since it appears from the figures that the elevation in energy expenditure is restricted to the awake cycle, we next asked if there was an *interaction* between genotype and the Light/Dark cycle.  Adding this interaction was highly significant 8.225&times; 10^-12^.  

The full results are shown below:


Table: Estimates and p-values from mixed linear models, excluding time of day.

|                                   | Estimate| Std..Error| t.value|   p.z|
|:----------------------------------|--------:|----------:|-------:|-----:|
|Lean                               |    0.013|      0.010|   1.303| 0.193|
|SexMale                            |   -0.002|      0.063|  -0.038| 0.970|
|KnockoutKnockout                   |    0.085|      0.044|   1.924| 0.054|
|`Light/Dark`Light:KnockoutKnockout |   -0.136|      0.019|  -7.188| 0.000|

# Lipid Oxidation Analysis

Calculated as $Lipid\ oxidation = (1.695 * vO_2) - (1.701 * vCO_2)$ where both units are in L/min and the output is in g/min

![Lipid Oxidation Summary Light/Dark Boxplot](figures/lipid-analysis-light-dark-1.png)

![](figures/lipid-analysis-linegraph, fig-1.png)<!-- -->

### Lipid Oxidation Summary Data

![Linegraph of Lipid Oxidation Data](figures/lipid-summarized-data-1.png)![Linegraph of Lipid Oxidation Data](figures/lipid-summarized-data-2.png)

### Lipid Oxidation Statistics


Table: Average changes in lipid oxidation comparing wt to knockout

|Sex    |Light/Dark | Control| Knockout| Change| Pct.Change|
|:------|:----------|-------:|--------:|------:|----------:|
|Female |Dark       |   0.389|    0.375| -0.015|      -3.80|
|Female |Light      |   0.377|    0.382|  0.006|       1.46|
|Male   |Dark       |   0.318|    0.309| -0.009|      -2.77|
|Male   |Light      |   0.370|    0.328| -0.041|     -11.15|

To test whether these groups are different we constructed a linear model with the following formula:

Lipid Oxidation ~ as.factor(Zeitgeber.Time) + Lean + Sex + `Light/Dark` + Knockout + Knockout:`Light/Dark` + (1 | Subject).  

We used this model because the base model was that lipid oxidation changes over the day.  We asked if lean mass modified the time dependent effect, but it did not (p=0.984).  We kept it in the model to be consistent with the carbohydrate oxidation.  After adjusting for lean mass, we asked if there was any additional benefit to including the light/dark cycle in addition to the time of day, and found that there was no significant effect, so that was not included in the initial model (p=NA).  We added sex as a covariate which had a highly significant effect 0. We next added knockout to the model and found no significant effect 0.378.  Finally we asked if Sex modified the effect of the knockout and found no significant effect 0.832.

The full results are shown below:


Table: Estimates and p-values from mixed linear models, excluding time of day.

|                                   | Estimate| Std..Error| t.value|   p.z|
|:----------------------------------|--------:|----------:|-------:|-----:|
|Lean                               |    0.012|      0.004|   2.781| 0.005|
|SexMale                            |   -0.102|      0.027|  -3.710| 0.000|
|KnockoutKnockout                   |   -0.018|      0.019|  -0.936| 0.349|
|`Light/Dark`Light:KnockoutKnockout |    0.002|      0.006|   0.415| 0.678|

### Lipid versus CHO Oxidation

![Comparason of lipid and carbohydrate oxidation rates in wild-type mice](figures/lipid-cho-oxidation-1.png)

## Activity Analysis

![Activity Summary Light/Dark Boxplot](figures/activity-analysis-light-dark-1.png)

![](figures/activity-analysis-linegraph, fig-1.png)<!-- -->

### Activity Summary Data

![Linegraph of Activity Data](figures/activity-summarized-data-1.png)![Linegraph of Activity Data](figures/activity-summarized-data-2.png)


# Interpretation

A brief summary of what the interpretation of these results were

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
##  [1] multcomp_1.4-18 TH.data_1.1-0   MASS_7.3-54     survival_3.2-13
##  [5] mvtnorm_1.1-3   lme4_1.1-27.1   Matrix_1.4-0    ggplot2_3.3.5  
##  [9] lubridate_1.8.0 readr_2.1.1     readxl_1.3.1    dplyr_1.0.7    
## [13] tidyr_1.1.4     knitr_1.37     
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.7       lattice_0.20-45  zoo_1.8-9        assertthat_0.2.1
##  [5] digest_0.6.29    utf8_1.2.2       R6_2.5.1         cellranger_1.1.0
##  [9] evaluate_0.14    highr_0.9        pillar_1.6.4     rlang_0.4.12    
## [13] minqa_1.2.4      jquerylib_0.1.4  nloptr_1.2.2.3   rmarkdown_2.11  
## [17] labeling_0.4.2   splines_4.0.2    stringr_1.4.0    bit_4.0.4       
## [21] munsell_0.5.0    compiler_4.0.2   xfun_0.29        pkgconfig_2.0.3 
## [25] mgcv_1.8-38      htmltools_0.5.2  tidyselect_1.1.1 tibble_3.1.6    
## [29] codetools_0.2-18 fansi_1.0.0      crayon_1.4.2     tzdb_0.2.0      
## [33] withr_2.4.3      grid_4.0.2       nlme_3.1-153     jsonlite_1.7.2  
## [37] gtable_0.3.0     lifecycle_1.0.1  DBI_1.1.2        magrittr_2.0.1  
## [41] scales_1.1.1     stringi_1.7.6    vroom_1.5.7      farver_2.1.0    
## [45] bslib_0.3.1      ellipsis_0.3.2   generics_0.1.1   vctrs_0.3.8     
## [49] sandwich_3.0-1   boot_1.3-28      tools_4.0.2      bit64_4.0.5     
## [53] glue_1.6.0       purrr_0.3.4      hms_1.1.1        parallel_4.0.2  
## [57] fastmap_1.1.0    yaml_2.2.1       colorspace_2.0-2 sass_0.4.0
```
