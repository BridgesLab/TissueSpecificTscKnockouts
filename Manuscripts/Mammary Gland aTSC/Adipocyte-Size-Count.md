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
## REML criterion at convergence: -1471
## 
## Scaled residuals: 
##    Min     1Q Median     3Q    Max 
## -2.926 -0.423 -0.126  0.508  3.937 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev. 
##  MouseID  (Intercept) 1.46e-09 0.0000382
##  Residual             1.61e-09 0.0000401
## Number of obs: 88, groups:  MouseID, 11
## 
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)  
## (Intercept) 0.0000323  0.0000182 8.9999993    1.77    0.110  
## GenotypeKO  0.0000539  0.0000247 8.9999993    2.19    0.057 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr)
## GenotypeKO -0.739
```

```
## [1] 167
```

```
## # A tibble: 88 x 6
## # Groups:   Image number, MouseID [88]
##    `Image number` MouseID Genotype adipocyttesum totalarea adipocytesperarea
##             <dbl>   <dbl> <fct>            <int>     <dbl>             <dbl>
##  1              1    7981 KO                 661  2576849.       0.000257   
##  2              1    7983 KO                  69  2321802.       0.0000297  
##  3              1    7984 KO                   1  2507579.       0.000000399
##  4              1    8161 KO                 293  2386745.       0.000123   
##  5              1    8162 WT                  34  2549742.       0.0000133  
##  6              1    8444 WT                  62  2382136.       0.0000260  
##  7              1    8445 WT                 162  2489745.       0.0000651  
##  8              1    8446 WT                 101  2523798.       0.0000400  
##  9              1    8465 KO                 422  2399398.       0.000176   
## 10              1    8466 KO                 176  2537593.       0.0000694  
## # … with 78 more rows
```

```
## # A tibble: 11 x 3
## # Groups:   MouseID [11]
##    MouseID Genotype adipocyttesum
##      <dbl> <fct>            <dbl>
##  1    7981 KO          0.000100  
##  2    7983 KO          0.0000944 
##  3    7984 KO          0.00000857
##  4    8161 KO          0.000115  
##  5    8162 WT          0.0000181 
##  6    8444 WT          0.0000502 
##  7    8445 WT          0.0000391 
##  8    8446 WT          0.0000397 
##  9    8465 KO          0.000157  
## 10    8466 KO          0.0000426 
## 11    8467 WT          0.0000143
```

```
## # A tibble: 2 x 3
##   Genotype adipocytenumbergenotpye se.adipocytenumbergenotpye
##   <fct>                      <dbl>                      <dbl>
## 1 WT                          32.3                       6.89
## 2 KO                          86.2                      21.6
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
## [1] 62.5
```



```
## # A tibble: 88 x 6
## # Groups:   Image number, MouseID [88]
##    `Image number` MouseID Genotype adipocytearea totalarea normalizedadipocytes…
##             <dbl>   <dbl> <fct>            <dbl>     <dbl>                 <dbl>
##  1              1    7981 KO              144166  2576849.                5.59  
##  2              1    7983 KO               24203  2321802.                1.04  
##  3              1    7984 KO                 259  2507579.                0.0103
##  4              1    8161 KO              115965  2386745.                4.86  
##  5              1    8162 WT               10390  2549742.                0.407 
##  6              1    8444 WT               37627  2382136.                1.58  
##  7              1    8445 WT               17122  2489745.                0.688 
##  8              1    8446 WT               30654  2523798.                1.21  
##  9              1    8465 KO              139326  2399398.                5.81  
## 10              1    8466 KO               59800  2537593.                2.36  
## # … with 78 more rows
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: normalizedadipocytesarea ~ Genotype + (1 | MouseID)
##    Data: adipocyte.percent
## 
## REML criterion at convergence: 310
## 
## Scaled residuals: 
##    Min     1Q Median     3Q    Max 
## -3.502 -0.353 -0.056  0.439  2.937 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  MouseID  (Intercept) 2.14     1.46    
##  Residual             1.52     1.23    
## Number of obs: 88, groups:  MouseID, 11
## 
## Fixed effects:
##             Estimate Std. Error    df t value Pr(>|t|)  
## (Intercept)    0.920      0.682 9.000    1.35    0.210  
## GenotypeKO     2.084      0.924 9.000    2.26    0.051 .
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

```
## # A tibble: 11 x 3
## # Groups:   MouseID [11]
##    MouseID Genotype adipocytearea
##      <dbl> <fct>            <dbl>
##  1    7981 KO               3.23 
##  2    7983 KO               2.87 
##  3    7984 KO               0.184
##  4    8161 KO               4.72 
##  5    8162 WT               0.589
##  6    8444 WT               1.32 
##  7    8445 WT               0.774
##  8    8446 WT               1.56 
##  9    8465 KO               5.57 
## 10    8466 KO               1.46 
## 11    8467 WT               0.352
```

```
## # A tibble: 2 x 3
##   Genotype adipocyteareagenotpye se.adipocyteareagenotpye
##   <fct>                    <dbl>                    <dbl>
## 1 WT                       0.920                    0.227
## 2 KO                       3.00                     0.815
```

![](figures/adipocyte percent of total MG-1.png)<!-- -->

```
## [1] 2.26
```

```
## # A tibble: 13,809 x 14
## # Groups:   Image number, MouseID, Genotype [88]
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
## # … with 13,799 more rows, and 7 more variables: Genotype <fct>,
## #   TotalAdipocyteNumber <dbl>, ...10 <lgl>, TotalImageArea <lgl>, `Total
## #   Area` <dbl>, Areainum2 <dbl>, avarea <dbl>
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: Areainum2 ~ Genotype + (1 | MouseID)
##    Data: adipocyte.areaperimage
## 
## REML criterion at convergence: 185965
## 
## Scaled residuals: 
##    Min     1Q Median     3Q    Max 
## -1.542 -0.698 -0.256  0.504 11.285 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  MouseID  (Intercept)  3266     57.2   
##  Residual             41232    203.1   
## Number of obs: 13809, groups:  MouseID, 11
## 
## Fixed effects:
##             Estimate Std. Error     df t value  Pr(>|t|)    
## (Intercept)   235.43      25.86   8.91    9.10 0.0000083 ***
## GenotypeKO     34.12      34.97   8.86    0.98      0.36    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##            (Intr)
## GenotypeKO -0.740
```

```
## # A tibble: 11 x 3
## # Groups:   MouseID [11]
##    MouseID Genotype adipocytearea
##      <dbl> <fct>            <dbl>
##  1    7981 KO                265.
##  2    7983 KO                251.
##  3    7984 KO                178.
##  4    8161 KO                341.
##  5    8162 WT                269.
##  6    8444 WT                216.
##  7    8445 WT                163.
##  8    8446 WT                326.
##  9    8465 KO                294.
## 10    8466 KO                282.
## 11    8467 WT                202.
```

```
## # A tibble: 2 x 3
##   Genotype adipocyteareagenotpye se.adipocyteareagenotpye
##   <fct>                    <dbl>                    <dbl>
## 1 WT                        235.                     28.3
## 2 KO                        269.                     22.0
```

![](figures/adipocyte area graph-1.png)<!-- -->


```
## # A tibble: 13,809 x 13
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
## # … with 13,799 more rows, and 6 more variables: Genotype <fct>,
## #   TotalAdipocyteNumber <dbl>, ...10 <lgl>, TotalImageArea <lgl>, `Total
## #   Area` <dbl>, Areaum2 <dbl>
```

![](figures/graphs for density area of adipocytes-1.png)<!-- -->

|term  |    df| statistic| p.value|
|:-----|-----:|---------:|-------:|
|group |     1|      19.7|       0|
|      | 13807|        NA|      NA|



| statistic| p.value|method                             |alternative |
|---------:|-------:|:----------------------------------|:-----------|
|     0.239|       0|Two-sample Kolmogorov-Smirnov test |two-sided   |

![](figures/graphs for density area of adipocytes-2.png)<!-- -->

```
## # A tibble: 22 x 5
## # Groups:   Genotype [2]
##    Genotype range       count totaladipocytes percentofadipocytes
##    <fct>    <fct>       <int>           <int>               <dbl>
##  1 WT       [0,100)      1249            3285              38.0  
##  2 WT       [100,200)     708            3285              21.6  
##  3 WT       [200,300)     395            3285              12.0  
##  4 WT       [300,400)     260            3285               7.91 
##  5 WT       [400,500)     186            3285               5.66 
##  6 WT       [500,600)     235            3285               7.15 
##  7 WT       [600,700)     100            3285               3.04 
##  8 WT       [700,800)      60            3285               1.83 
##  9 WT       [800,900)      38            3285               1.16 
## 10 WT       [900,1e+03)    16            3285               0.487
## # … with 12 more rows
```

![](figures/graphs for density area of adipocytes-3.png)<!-- -->

```
## # A tibble: 22 x 4
## # Groups:   Genotype [2]
##    Genotype range       averagepercentofadipocytes se.averagepercentofadipocytes
##    <fct>    <fct>                            <dbl>                         <dbl>
##  1 WT       [0,100)                         35.3                           7.36 
##  2 WT       [100,200)                       22.9                           3.35 
##  3 WT       [200,300)                       12.9                           1.95 
##  4 WT       [300,400)                        9.08                          2.09 
##  5 WT       [400,500)                        6.12                          1.25 
##  6 WT       [500,600)                        6.98                          1.69 
##  7 WT       [600,700)                        2.86                          0.649
##  8 WT       [700,800)                        1.58                          0.626
##  9 WT       [800,900)                        1.20                          0.501
## 10 WT       [900,1e+03)                      0.432                         0.120
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
## -18.975  -1.748  -0.162   1.570  26.285 
## 
## Coefficients:
##                       Estimate Std. Error t value Pr(>|t|)    
## (Intercept)           24.97789    2.05107   12.18  < 2e-16 ***
## GenotypeKO             0.00147    1.19343    0.00   0.9990    
## range[100,200)         2.25776    2.75069    0.82   0.4136    
## range[200,300)        -8.38722    2.75069   -3.05   0.0029 ** 
## range[300,400)       -14.34692    2.75069   -5.22  9.1e-07 ***
## range[400,500)       -18.25888    2.75069   -6.64  1.4e-09 ***
## range[500,600)       -17.86826    2.75069   -6.50  2.7e-09 ***
## range[600,700)       -21.81366    2.75069   -7.93  2.4e-12 ***
## range[700,800)       -23.40603    2.75069   -8.51  1.2e-13 ***
## range[800,900)       -24.02490    2.81937   -8.52  1.2e-13 ***
## range[900,1e+03)     -24.54465    2.81914   -8.71  4.5e-14 ***
## range[1e+03,3.5e+03) -24.17344    2.81914   -8.57  8.9e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 6.45 on 106 degrees of freedom
## Multiple R-squared:  0.699,	Adjusted R-squared:  0.668 
## F-statistic: 22.4 on 11 and 106 DF,  p-value: <2e-16
```

```
## # A tibble: 11 x 3
##    range             pval delta
##    <fct>            <dbl> <dbl>
##  1 [0,100)         0.0601 0.464
##  2 [100,200)       0.228  1.35 
##  3 [200,300)       0.0394 1.53 
##  4 [300,400)       0.374  1.31 
##  5 [400,500)       0.575  1.18 
##  6 [500,600)       0.927  1.03 
##  7 [600,700)       0.492  1.20 
##  8 [700,800)       0.995  0.997
##  9 [800,900)       0.485  0.663
## 10 [900,1e+03)     0.983  1.01 
## 11 [1e+03,3.5e+03) 0.602  0.686
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
