---
title: "aTSC PND16.5 Mammary Gland Adipocyte Sizing and Counting"
author: "Allison Meyer"
date: "2/21/2021"
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



## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:



# Purpose
write purpose here for each file

# Experimental Details

Link to the protocol used (permalink preferred) for the experiment and include any notes relevant to your analysis.  This might include specifics not in the general protocol such as cell lines, treatment doses etc.

# Raw Data




```
## # A tibble: 88 x 3
##    `Image#` MouseID `Total Area`
##       <dbl>   <dbl>        <dbl>
##  1        1    7981      3117987
##  2        2    7981      3145728
##  3        3    7981      3120420
##  4        4    7981      3145728
##  5        5    7981      3005199
##  6        6    7981      3145728
##  7        7    7981      3145728
##  8        8    7981      2264574
##  9        1    7983      2809380
## 10        2    7983      3145728
## # … with 78 more rows
```

![](figures/histogram for adipocyte area distribution per mouse and combined for all dams-1.png)<!-- -->![](figures/histogram for adipocyte area distribution per mouse and combined for all dams-2.png)<!-- -->![](figures/histogram for adipocyte area distribution per mouse and combined for all dams-3.png)<!-- -->![](figures/histogram for adipocyte area distribution per mouse and combined for all dams-4.png)<!-- -->![](figures/histogram for adipocyte area distribution per mouse and combined for all dams-5.png)<!-- -->![](figures/histogram for adipocyte area distribution per mouse and combined for all dams-6.png)<!-- -->![](figures/histogram for adipocyte area distribution per mouse and combined for all dams-7.png)<!-- -->![](figures/histogram for adipocyte area distribution per mouse and combined for all dams-8.png)<!-- -->![](figures/histogram for adipocyte area distribution per mouse and combined for all dams-9.png)<!-- -->![](figures/histogram for adipocyte area distribution per mouse and combined for all dams-10.png)<!-- -->![](figures/histogram for adipocyte area distribution per mouse and combined for all dams-11.png)<!-- -->

```
## # A tibble: 13,809 x 11
##    MouseID `Image number` `Adipocyte numb…  Area `Mean Area`   Min   Max
##      <dbl>          <dbl>            <dbl> <dbl>       <dbl> <dbl> <dbl>
##  1    7981              1                1   185        195.   185   205
##  2    7981              1                2   100        193.   182   200
##  3    7981              1                3   119        161.   153   176
##  4    7981              1                4   171        193.   171   202
##  5    7981              1                5   101        161.   150   174
##  6    7981              1                6   455        192.   169   205
##  7    7981              1                7    88        203.   196   211
##  8    7981              1                8   181        192.   178   201
##  9    7981              1                9   111        205.   197   212
## 10    7981              1               10   109        165.   156   177
## # … with 13,799 more rows, and 4 more variables: Genotype <fct>,
## #   TotalAdipocyteNumber <dbl>, ...10 <lgl>, TotalImageArea <lgl>
```

```
## # A tibble: 13,809 x 12
##    MouseID `Image number` `Adipocyte numb…  Area `Mean Area`   Min   Max
##      <dbl>          <dbl>            <dbl> <dbl>       <dbl> <dbl> <dbl>
##  1    7981              1                1   185        195.   185   205
##  2    7981              1                2   100        193.   182   200
##  3    7981              1                3   119        161.   153   176
##  4    7981              1                4   171        193.   171   202
##  5    7981              1                5   101        161.   150   174
##  6    7981              1                6   455        192.   169   205
##  7    7981              1                7    88        203.   196   211
##  8    7981              1                8   181        192.   178   201
##  9    7981              1                9   111        205.   197   212
## 10    7981              1               10   109        165.   156   177
## # … with 13,799 more rows, and 5 more variables: Genotype <fct>,
## #   TotalAdipocyteNumber <dbl>, ...10 <lgl>, TotalImageArea <lgl>, `Total
## #   Area` <dbl>
```


```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: adipocytesperarea ~ Genotype + (1 | MouseID)
##    Data: adipocyte.number
## 
## REML criterion at convergence: -1503
## 
## Scaled residuals: 
##    Min     1Q Median     3Q    Max 
## -2.926 -0.423 -0.126  0.508  3.937 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. 
##  MouseID  (Intercept) 9.95e-10 0.0000315
##  Residual             1.10e-09 0.0000332
## Number of obs: 88, groups:  MouseID, 11
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)  
## (Intercept) 0.0000267  0.0000151 8.9999993    1.77    0.110  
## GenotypeKO  0.0000446  0.0000204 8.9999993    2.19    0.057 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr)
## GenotypeKO -0.739
```

```
## # A tibble: 88 x 6
## # Groups:   Image number, MouseID [88]
##    `Image number` MouseID Genotype adipocyttesum totalarea adipocytesperarea
##             <dbl>   <dbl> <fct>            <int>     <dbl>             <dbl>
##  1              1    7981 KO                 661   3117987       0.000212   
##  2              1    7983 KO                  69   2809380       0.0000246  
##  3              1    7984 KO                   1   3034170       0.000000330
##  4              1    8161 KO                 293   2887962       0.000101   
##  5              1    8162 WT                  34   3085188       0.0000110  
##  6              1    8444 WT                  62   2882385       0.0000215  
##  7              1    8445 WT                 162   3012591       0.0000538  
##  8              1    8446 WT                 101   3053796       0.0000331  
##  9              1    8465 KO                 422   2903272       0.000145   
## 10              1    8466 KO                 176   3070488       0.0000573  
## # … with 78 more rows
```

![](figures/Adipocyte number normalized to total MG area, and stats for count and area-1.png)<!-- -->

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: Area ~ Genotype + (1 | MouseID)
##    Data: allmergeddata
## 
## REML criterion at convergence: 191228
## 
## Scaled residuals: 
##    Min     1Q Median     3Q    Max 
## -1.542 -0.698 -0.256  0.504 11.285 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  MouseID  (Intercept)  4782     69.2   
##  Residual             60367    245.7   
## Number of obs: 13809, groups:  MouseID, 11
## 
## Fixed effects:
##             Estimate Std. Error     df t value  Pr(>|t|)    
## (Intercept)   284.87      31.30   8.91    9.10 0.0000083 ***
## GenotypeKO     41.28      42.31   8.86    0.98      0.36    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr)
## GenotypeKO -0.740
```

```
## Generalized linear mixed model fit by maximum likelihood (Laplace
##   Approximation) [glmerMod]
##  Family: poisson  ( log )
## Formula: Area ~ Genotype + (1 | MouseID)
##    Data: allmergeddata
## 
##      AIC      BIC   logLik deviance df.resid 
##  2351179  2351202 -1175587  2351173    13806 
## 
## Scaled residuals: 
##    Min     1Q Median     3Q    Max 
## -19.16  -9.79  -3.48   6.83 186.94 
## 
## Random effects:
##  Groups  Name        Variance Std.Dev.
##  MouseID (Intercept) 0.0474   0.218   
## Number of obs: 13809, groups:  MouseID, 11
## 
## Fixed effects:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   5.6234     0.0966   58.20   <2e-16 ***
## GenotypeKO    0.1419     0.1310    1.08     0.28    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr)
## GenotypeKO -0.735
```



```
## # A tibble: 88 x 6
## # Groups:   Image number, MouseID [88]
##    `Image number` MouseID Genotype adipocytearea totalarea normalizedadipocytes…
##             <dbl>   <dbl> <fct>            <dbl>     <dbl>                 <dbl>
##  1              1    7981 KO              144166   3117987               4.62   
##  2              1    7983 KO               24203   2809380               0.862  
##  3              1    7984 KO                 259   3034170               0.00854
##  4              1    8161 KO              115965   2887962               4.02   
##  5              1    8162 WT               10390   3085188               0.337  
##  6              1    8444 WT               37627   2882385               1.31   
##  7              1    8445 WT               17122   3012591               0.568  
##  8              1    8446 WT               30654   3053796               1.00   
##  9              1    8465 KO              139326   2903272               4.80   
## 10              1    8466 KO               59800   3070488               1.95   
## # … with 78 more rows
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: normalizedadipocytesarea ~ Genotype + (1 | MouseID)
##    Data: adipocyte.percent
## 
## REML criterion at convergence: 278
## 
## Scaled residuals: 
##    Min     1Q Median     3Q    Max 
## -3.502 -0.353 -0.056  0.439  2.937 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  MouseID  (Intercept) 1.46     1.21    
##  Residual             1.04     1.02    
## Number of obs: 88, groups:  MouseID, 11
## 
## Fixed effects:
##             Estimate Std. Error    df t value Pr(>|t|)  
## (Intercept)    0.761      0.564 9.000    1.35    0.210  
## GenotypeKO     1.723      0.764 9.000    2.26    0.051 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr)
## GenotypeKO -0.739
```

```
## [1] 226
```

![](figures/adipocyte percent of total MG-1.png)<!-- -->

![](figures/graphs for density area of adipocytes-1.png)<!-- -->

|term  |    df| statistic| p.value|
|:-----|-----:|---------:|-------:|
|group |     1|      19.7|       0|
|      | 13807|        NA|      NA|

![](figures/graphs for density area of adipocytes-2.png)<!-- -->

```
## # A tibble: 22 x 5
## # Groups:   Genotype [2]
##    Genotype range       count totaladipocytes percentofadipocytes
##    <fct>    <fct>       <int>           <int>               <dbl>
##  1 WT       [0,100)      1019            3285               31.0 
##  2 WT       [100,200)     741            3285               22.6 
##  3 WT       [200,300)     411            3285               12.5 
##  4 WT       [300,400)     267            3285                8.13
##  5 WT       [400,500)     208            3285                6.33
##  6 WT       [500,600)     142            3285                4.32
##  7 WT       [600,700)     219            3285                6.67
##  8 WT       [700,800)      94            3285                2.86
##  9 WT       [800,900)      54            3285                1.64
## 10 WT       [900,1e+03)    53            3285                1.61
## # … with 12 more rows
```

![](figures/graphs for density area of adipocytes-3.png)<!-- -->

```
## # A tibble: 22 x 4
## # Groups:   Genotype [2]
##    Genotype range       averagepercentofadipocytes se.averagepercentofadipocytes
##    <fct>    <fct>                            <dbl>                         <dbl>
##  1 WT       [0,100)                          29.1                          6.46 
##  2 WT       [100,200)                        22.2                          2.74 
##  3 WT       [200,300)                        13.8                          2.37 
##  4 WT       [300,400)                         8.87                         1.81 
##  5 WT       [400,500)                         7.27                         1.53 
##  6 WT       [500,600)                         4.56                         0.877
##  7 WT       [600,700)                         6.49                         1.55 
##  8 WT       [700,800)                         2.96                         0.612
##  9 WT       [800,900)                         1.37                         0.521
## 10 WT       [900,1e+03)                       1.37                         0.513
## # … with 12 more rows
```

![](figures/graphs for density area of adipocytes-4.png)<!-- -->

```
## 
## Call:
## lm(formula = percentofadipocytes ~ Genotype + range, data = percentileareabymouseonly)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -15.874  -1.793  -0.205   1.714  23.960 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)           19.4305     1.8482   10.51  < 2e-16 ***
## GenotypeKO             0.0199     1.0663    0.02    0.985    
## range[100,200)         4.8847     2.4809    1.97    0.052 .  
## range[200,300)        -2.5282     2.4809   -1.02    0.310    
## range[300,400)        -7.8738     2.4809   -3.17    0.002 ** 
## range[400,500)       -11.1724     2.4809   -4.50  1.7e-05 ***
## range[500,600)       -14.1642     2.4809   -5.71  1.0e-07 ***
## range[600,700)       -12.8713     2.4809   -5.19  1.0e-06 ***
## range[700,800)       -16.3065     2.4809   -6.57  1.8e-09 ***
## range[800,900)       -17.7873     2.4809   -7.17  9.8e-11 ***
## range[900,1e+03)     -18.1396     2.5427   -7.13  1.2e-10 ***
## range[1e+03,3.5e+03) -17.7770     2.4809   -7.17  1.0e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 5.82 on 108 degrees of freedom
## Multiple R-squared:  0.654,	Adjusted R-squared:  0.619 
## F-statistic: 18.6 on 11 and 108 DF,  p-value: <2e-16
```

```
## # A tibble: 11 x 3
##    range             pval delta
##    <fct>            <dbl> <dbl>
##  1 [0,100)         0.0495 0.393
##  2 [100,200)       0.514  1.17 
##  3 [200,300)       0.0804 1.41 
##  4 [300,400)       0.0924 1.56 
##  5 [400,500)       0.435  1.25 
##  6 [500,600)       0.377  1.29 
##  7 [600,700)       0.952  1.02 
##  8 [700,800)       0.706  1.11 
##  9 [800,900)       0.407  1.39 
## 10 [900,1e+03)     0.815  0.898
## 11 [1e+03,3.5e+03) 0.604  0.719
```



![](figures/Total Adipocyte Number and Area per Mouse-1.png)<!-- -->![](figures/Total Adipocyte Number and Area per Mouse-2.png)<!-- -->![](figures/Total Adipocyte Number and Area per Mouse-3.png)<!-- -->![](figures/Total Adipocyte Number and Area per Mouse-4.png)<!-- -->![](figures/Total Adipocyte Number and Area per Mouse-5.png)<!-- -->

![](figures/adipocyte size, count grouped by genoytpe-1.png)<!-- -->![](figures/adipocyte size, count grouped by genoytpe-2.png)<!-- -->![](figures/adipocyte size, count grouped by genoytpe-3.png)<!-- -->![](figures/adipocyte size, count grouped by genoytpe-4.png)<!-- -->![](figures/adipocyte size, count grouped by genoytpe-5.png)<!-- -->


```
## 
## 	Shapiro-Wilk normality test
## 
## data:  adipocyte.mouseID$`Total Adipocyte Number`
## W = 0.9, p-value = 0.2
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  filter(adipocyte.mouseID, Genotype == "KO")$`Total Adipocyte Number`
## W = 1, p-value = 0.8
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  filter(adipocyte.mouseID, Genotype == "WT")$`Total Adipocyte Number`
## W = 0.9, p-value = 0.4
```



Table: Levene test for equality of variances

|term  | df| statistic| p.value|
|:-----|--:|---------:|-------:|
|group |  1|      2.29|   0.165|
|      |  9|        NA|      NA|



| estimate| estimate1| estimate2| statistic| p.value| parameter| conf.low| conf.high|method            |alternative |
|--------:|---------:|---------:|---------:|-------:|---------:|--------:|---------:|:-----------------|:-----------|
|    -1097|       657|      1754|      -2.2|   0.055|         9|    -2226|      31.7|Two Sample t-test |two.sided   |

```
## 
## Call:
## lm(formula = `Total Adipocyte Number` ~ Genotype, data = adipocyte.mouseID)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
##  -1578   -323    165    323   1445 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept)      657        368    1.78    0.108  
## GenotypeKO      1097        499    2.20    0.055 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 824 on 9 degrees of freedom
## Multiple R-squared:  0.349,	Adjusted R-squared:  0.277 
## F-statistic: 4.83 on 1 and 9 DF,  p-value: 0.0555
```

```
## 
## Call:
## lm(formula = TotalAdipocyteNumber ~ Genotype, data = combined.data)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
##  -1580   -323    165    327   1443 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept)      657        369    1.78    0.108  
## GenotypeKO      1099        499    2.20    0.055 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 824 on 9 degrees of freedom
##   (13798 observations deleted due to missingness)
## Multiple R-squared:  0.35,	Adjusted R-squared:  0.278 
## F-statistic: 4.84 on 1 and 9 DF,  p-value: 0.0553
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  adipocyte.mouseID$`Adipocyte Avergae Total Area`
## W = 0.8, p-value = 0.03
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  filter(adipocyte.mouseID, Genotype == "KO")$`Adipocyte Avergae Total Area`
## W = 0.9, p-value = 0.7
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  filter(adipocyte.mouseID, Genotype == "WT")$`Adipocyte Avergae Total Area`
## W = 0.9, p-value = 0.7
```



Table: Levene test for equality of variances of 

|term  | df| statistic| p.value|
|:-----|--:|---------:|-------:|
|group |  1|      8.18|   0.019|
|      |  9|        NA|      NA|



| estimate| estimate1| estimate2| statistic| p.value| parameter| conf.low| conf.high|method                  |alternative |
|--------:|---------:|---------:|---------:|-------:|---------:|--------:|---------:|:-----------------------|:-----------|
|  -376442|    187315|    563758|     -2.12|   0.081|      5.71|  -816747|     63862|Welch Two Sample t-test |two.sided   |

```
## 
## Call:
## lm(formula = `Adipocyte Avergae Total Area` ~ Genotype, data = adipocyte.mouseID)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -525817 -182736  -30665  116291  573199 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept)   187315     143504    1.31    0.224  
## GenotypeKO    376442     194305    1.94    0.085 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 321000 on 9 degrees of freedom
## Multiple R-squared:  0.294,	Adjusted R-squared:  0.216 
## F-statistic: 3.75 on 1 and 9 DF,  p-value: 0.0847
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: Area ~ Genotype + (1 | MouseID)
##    Data: combined.data
## 
## REML criterion at convergence: 191228
## 
## Scaled residuals: 
##    Min     1Q Median     3Q    Max 
## -1.542 -0.698 -0.256  0.504 11.285 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  MouseID  (Intercept)  4782     69.2   
##  Residual             60367    245.7   
## Number of obs: 13809, groups:  MouseID, 11
## 
## Fixed effects:
##             Estimate Std. Error     df t value  Pr(>|t|)    
## (Intercept)   284.87      31.30   8.91    9.10 0.0000083 ***
## GenotypeKO     41.28      42.31   8.86    0.98      0.36    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr)
## GenotypeKO -0.740
```


```
## 
## 	Shapiro-Wilk normality test
## 
## data:  adipocyte.mouseID$`Adipocyte Average Percent Area`
## W = 0.8, p-value = 0.01
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  filter(adipocyte.mouseID, Genotype == "KO")$`Adipocyte Average Percent Area`
## W = 0.9, p-value = 0.7
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  filter(adipocyte.mouseID, Genotype == "WT")$`Adipocyte Average Percent Area`
## W = 0.9, p-value = 0.6
```



Table: Levene test for equality of variances of 

|term  | df| statistic| p.value|
|:-----|--:|---------:|-------:|
|group |  1|      11.6|   0.008|
|      |  9|        NA|      NA|



| estimate| estimate1| estimate2| statistic| p.value| parameter| conf.low| conf.high|method                  |alternative |
|--------:|---------:|---------:|---------:|-------:|---------:|--------:|---------:|:-----------------------|:-----------|
|    -1.69|      0.62|      2.31|     -2.36|   0.062|      5.34|     -3.5|     0.118|Welch Two Sample t-test |two.sided   |

```
## 
## Call:
## lm(formula = `Adipocyte Average Percent Area` ~ Genotype, data = adipocyte.mouseID)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -2.158 -0.684 -0.032  0.445  2.293 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept)    0.620      0.582    1.06    0.315  
## GenotypeKO     1.691      0.789    2.14    0.061 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.3 on 9 degrees of freedom
## Multiple R-squared:  0.338,	Adjusted R-squared:  0.265 
## F-statistic:  4.6 on 1 and 9 DF,  p-value: 0.0606
```




```
## 
## 	Shapiro-Wilk normality test
## 
## data:  adipocyte.mouseID$`Sum Adipocyte Area`
## W = 0.9, p-value = 0.1
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  filter(adipocyte.mouseID, Genotype == "KO")$`Sum Adipocyte Area`
## W = 1, p-value = 1
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  filter(adipocyte.mouseID, Genotype == "WT")$`Sum Adipocyte Area`
## W = 0.9, p-value = 0.7
```



Table: Levene test for equality of variances of 

|term  | df| statistic| p.value|
|:-----|--:|---------:|-------:|
|group |  1|      4.17|   0.071|
|      |  9|        NA|      NA|



| estimate| estimate1| estimate2| statistic| p.value| parameter| conf.low| conf.high|method            |alternative |
|--------:|---------:|---------:|---------:|-------:|---------:|--------:|---------:|:-----------------|:-----------|
|  -422891|    187481|    610371|     -2.27|   0.049|         9|  -844298|     -1483|Two Sample t-test |two.sided   |

```
## 
## Call:
## lm(formula = `Sum Adipocyte Area` ~ Genotype, data = adipocyte.mouseID)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -572431  -90554  -16789  105968  526232 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept)   187481     137581    1.36    0.206  
## GenotypeKO    422891     186286    2.27    0.049 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 308000 on 9 degrees of freedom
## Multiple R-squared:  0.364,	Adjusted R-squared:  0.293 
## F-statistic: 5.15 on 1 and 9 DF,  p-value: 0.0494
```


```
## 
## 	Shapiro-Wilk normality test
## 
## data:  adipocyte.mouseID$AverageAdipocyteArea
## W = 1, p-value = 0.8
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  filter(adipocyte.mouseID, Genotype == "KO")$AverageAdipocyteArea
## W = 1, p-value = 0.8
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  filter(adipocyte.mouseID, Genotype == "WT")$AverageAdipocyteArea
## W = 1, p-value = 0.8
```



Table: Levene test for equality of variances of 

|term  | df| statistic| p.value|
|:-----|--:|---------:|-------:|
|group |  1|     0.146|   0.712|
|      |  9|        NA|      NA|



| estimate| estimate1| estimate2| statistic| p.value| parameter| conf.low| conf.high|method            |alternative |
|--------:|---------:|---------:|---------:|-------:|---------:|--------:|---------:|:-----------------|:-----------|
|    -40.4|       285|       325|     -0.95|   0.367|         9|     -137|      55.8|Two Sample t-test |two.sided   |

```
## 
## Call:
## lm(formula = AverageAdipocyteArea ~ Genotype, data = adipocyte.mouseID)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -108.43  -31.77   -4.19   35.49  109.65 
## 
## Coefficients:
##             Estimate Std. Error t value  Pr(>|t|)    
## (Intercept)    284.8       31.4    9.06 0.0000081 ***
## GenotypeKO      40.4       42.5    0.95      0.37    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 70.3 on 9 degrees of freedom
## Multiple R-squared:  0.0912,	Adjusted R-squared:  -0.00983 
## F-statistic: 0.903 on 1 and 9 DF,  p-value: 0.367
```
## Including Plots

You can also embed plots, for example:

![](figures/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
