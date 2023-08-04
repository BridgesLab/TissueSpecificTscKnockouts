---
title: "Ketone Tolerance in CD/KD A/J Mice"
author: "Cody Cousineau, Detrick Snyder and Dave Bridges"
date: "August 8, 2019"
output:
  html_document:
    toc: yes
    keep_md: yes
  pdf_document:
    highlight: tango
    keep_tex: yes
    number_sections: no
    toc: yes
---



This script was most recently run on Fri Aug  4 09:18:09 2023 and can be found in /Users/davebrid/Documents/GitHub/TissueSpecificTscKnockouts/Mouse Data/AJ Ketogenic Diet.

# Purpose


# Experimental Details

Beta hydroxybutyrate was injected at 1g/kg into fed mice around 14:00 PM, ketone levels were assessed via tail vein injection using strips.

# Raw Data



Table: Animals who underwent an KTT

| MouseID|Diet           |  FK| age|
|-------:|:--------------|---:|---:|
|    8063|Control Diet   | 0.5|  72|
|    8064|Control Diet   | 0.4|  72|
|    8065|Control Diet   | 0.4|  72|
|    8066|Ketogenic Diet | 0.6|  72|
|    8067|Ketogenic Diet | 3.2|  72|
|    8068|Ketogenic Diet | 0.5|  72|
|    8069|Ketogenic Diet | 0.4|  72|
|    8070|Ketogenic Diet | 0.4|  72|
|    8071|Ketogenic Diet | 0.4|  72|
|    8072|Ketogenic Diet | 0.6|  72|
|    8073|Ketogenic Diet | 0.9|  72|
|    8074|Ketogenic Diet | 0.6|  72|
|    8063|Control Diet   | 0.4|  92|
|    8064|Control Diet   | 0.4|  92|
|    8065|Control Diet   | 0.5|  92|
|    8066|Ketogenic Diet | 0.6|  92|
|    8068|Ketogenic Diet | 0.8|  92|
|    8069|Ketogenic Diet | 0.7|  92|
|    8070|Ketogenic Diet | 0.9|  92|
|    8071|Ketogenic Diet | 0.8|  92|
|    8072|Ketogenic Diet | 0.7|  92|
|    8073|Ketogenic Diet | 0.8|  92|
|    8074|Ketogenic Diet | 0.7|  92|



Table: Animals evaluated by KTT in each group

|Diet           | age|  n|
|:--------------|---:|--:|
|Control Diet   |  72|  3|
|Control Diet   |  92|  3|
|Ketogenic Diet |  72|  9|
|Ketogenic Diet |  92|  8|



# Analysis



## Fed Ketone Levels

![Fed ketones, by treatment](figures/fasting-ketone-boxplot-itt-1.png)


![Fed ketones](figures/fasting-ketone-barplot-itt-ko-1.png)

### Fasting ketone Statistics


Table: ANOVA for effects of fasting ketone

|term      | df| sumsq| meansq| statistic| p.value|
|:---------|--:|-----:|------:|---------:|-------:|
|Diet      |  1| 0.596|  0.596|     1.739|   0.203|
|age       |  1| 0.028|  0.028|     0.081|   0.779|
|Diet:age  |  1| 0.010|  0.010|     0.029|   0.867|
|Residuals | 19| 6.516|  0.343|        NA|      NA|

## Ketone Tolerance Test

![Ketone tolerance tests by treatments](figures/ktt-dotplot-1.png)![Ketone tolerance tests by treatments](figures/ktt-dotplot-2.png)![Ketone tolerance tests by treatments](figures/ktt-dotplot-3.png)

![](figures/ktt-lineplot-1.png)<!-- -->

## Area Under the Curve

![](figures/ktt-auc-1.png)<!-- -->

## Normalized Area Under the Curve

![](figures/ktt-auc-norm-1.png)<!-- -->


## Normalized to fasting ketone levels (percent change)

![](figures/normalized-ktt-1.png)<!-- -->![](figures/normalized-ktt-2.png)<!-- -->

## Normlized to the Peak of Ketone levels

![](figures/ktt-peak-1.png)<!-- -->

## KTT Statistics


Table: Mixed linear model for KTT

|term            |  sumsq| meansq| NumDF| DenDF| statistic| p.value|
|:---------------|------:|------:|-----:|-----:|---------:|-------:|
|as.factor(Time) | 69.398| 11.566|     6|   141|    30.300|   0.000|
|as.factor(age)  |  1.744|  1.744|     1|   141|     4.568|   0.034|
|Diet            |  0.038|  0.038|     1|    46|     0.098|   0.755|
|Diet:age        |  0.217|  0.217|     1|   141|     0.567|   0.453|



Table: Mixed linear model for KTT at 3 weeks

|term                 |  sumsq| meansq| NumDF| DenDF| statistic| p.value|
|:--------------------|------:|------:|-----:|-----:|---------:|-------:|
|as.factor(Time)      | 808831| 134805|     6|    54|    36.164|   0.000|
|Diet                 |   2525|   2525|     1|     9|     0.677|   0.432|
|as.factor(Time):Diet |  15014|   2502|     6|    54|     0.671|   0.673|



Table: Mixed linear model for KTT at 3 weeks

|term                 | npar| AIC| BIC| logLik| deviance| statistic| df| p.value|
|:--------------------|----:|---:|---:|------:|--------:|---------:|--:|-------:|
|ktt.lme.92d.fil.null |    9| 689| 708|   -335|      671|        NA| NA|      NA|
|ktt.lme.92d.fil      |   16| 674| 708|   -321|      642|      28.6|  7|       0|



Table: Mixed linear model for KTT at 3 weeks

|names                                |      x|
|:------------------------------------|------:|
|(Intercept)                          |  100.0|
|as.factor(Time)15                    |  403.3|
|as.factor(Time)30                    |  363.3|
|as.factor(Time)45                    |  301.7|
|as.factor(Time)60                    |  215.0|
|as.factor(Time)75                    |  228.3|
|as.factor(Time)90                    |  243.3|
|DietKetogenic Diet                   |    0.0|
|as.factor(Time)15:DietKetogenic Diet |  -45.6|
|as.factor(Time)30:DietKetogenic Diet |  -50.7|
|as.factor(Time)45:DietKetogenic Diet |  -47.2|
|as.factor(Time)60:DietKetogenic Diet |  -15.5|
|as.factor(Time)75:DietKetogenic Diet |  -64.3|
|as.factor(Time)90:DietKetogenic Diet | -105.0|

# Interpretation

# References

