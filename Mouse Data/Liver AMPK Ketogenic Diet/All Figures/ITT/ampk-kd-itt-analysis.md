---
title: "Analysis of ITT's from AMPK KD Study"
author: "Katherine Kistler, Cody Cousineau, and Dave Bridges"
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

Injected 0.75 mU/g of lean body mass to animals fasted for 6h after 2 weeks of KD/CD (or four weeks after the injection)

# Raw Data


```r
library(readxl) #loads the readr package

data.filename <- 'ITT 19-06-19.xlsx' #make this a separate line, you can use any variable you want
not.ko <- c(8830,8831,8715,8716,8717,8718,8719,8786)
exp.data <- read_excel(data.filename) %>% 
    filter(!(Mouse %in% not.ko)) #removed mice that we are not confident of their ko status

mapping.filename <- 'mapping.csv'
#this loads whatever the file is into a dataframe called exp.data if it exists
exp.data <- read_excel(data.filename, sheet='Sheet2') %>%
  select(Group:`120`) %>%
  filter(!(Group %in% c("Average","SE",NA))) %>% # removed calculated averages and errors
  separate(Group,
           into=c("Sex","Diet","Injection"), # broke group down into sex, diet and injection 
           sep=" ") %>%
  tibble::rowid_to_column("ID") %>% #index an ID column
  mutate(ID = as.factor(ID)) %>%
  mutate(Injection = relevel(as.factor(Injection), ref="GFP")) 

data.long <-
  exp.data %>%
  group_by(Sex,Diet,Injection) %>%
  gather(key=Time,value=Glucose, -Sex, -Diet, -Injection, -ID,-Mouse) %>%
  mutate(Time = as.integer(Time),
         Glucose = as.integer(Glucose))

summary.data <-
  data.long %>%
  group_by(Sex,Diet,Injection, Time) %>%
  mutate(Mouse = as.factor(Mouse)) %>%
  summarize_if(is.numeric, .funs=funs(mean=mean(., na.rm = TRUE),
                                      se=se))
```

These data can be found in **/Users/katherinekistler/Documents/GitHub/TissueSpecificTscKnockouts/Mouse Data/Liver AMPK Ketogenic Diet/All Figures/ITT** in a file named **ITT 19-06-19.xlsx** and **mapping.csv**.  This script was most recently updated on **Mon Jun 22 10:32:23 2020**.

# Number of Mice


```r
exp.data %>%
  group_by(Sex,Diet,Injection) %>%
  distinct(ID, .keep_all = T) %>%
  count %>%
  kable(caption="Animals in each group of this cohort")
```



Table: Animals in each group of this cohort

Sex   Diet      Injection     n
----  --------  ----------  ---
F     Control   GFP           8
F     Control   Cre           7
F     Keto      GFP           7
F     Keto      Cre           8
M     Control   GFP           9
M     Control   Cre           8
M     Keto      GFP           9
M     Keto      Cre          11

# Analysis


```r
library(ggplot2)

ggplot(data.long, aes(y=Glucose,
                      x=Time,
                      col=Injection)) +
  geom_point() +
  geom_smooth() +
  facet_grid(Sex~Diet) +
  labs(y="Blood Glucose (mg/dL)",
       x="Insulin (min)")
```

![Dotplot of ITT, all groups](figures/itt-points-all-1.png)



```r
data.long %>%
  group_by(Sex,Diet,Injection,Time) %>%
  summarize(Mean = mean(Glucose,na.rm=T),
            Error = se(Glucose)) %>%
  ggplot(aes(y=Mean,
             x=Time,
             col=Injection)) +
  geom_line() +
  geom_errorbar(aes(ymin=Mean-Error,
             ymax=Mean+Error)) +
  expand_limits(y=0) +
  facet_grid(Sex~Diet) +
  labs(y="Blood Glucose (mg/dL)",
       x="Insulin (min)")
```

![Lineplot of ITT, all mice](figures/itt-line-all-1.png)


```r
data.long %>%
  filter(Sex=="M") %>%
  group_by(Sex,Diet,Injection,Time) %>%
  summarize(Mean = mean(Glucose,na.rm=T),
            Error = se(Glucose)) %>%
  ggplot(aes(y=Mean,
             x=Time,
             ymin=Mean-Error,
             ymax=Mean+Error,
             col=Injection)) +
  geom_line() +
  geom_errorbar() +
  expand_limits(y=0) +
  scale_color_brewer(palette = "Set2",
                     name="", 
                     labels=c("AAV-Tbg-GFP",
                              "AAV-Tbg-Cre")) +
  facet_grid(~Diet) +
  labs(y="Blood Glucose (mg/dL)",
       x="Insulin (min)",
       title="Hepatic AMPK Knockout Mice") +
  theme_classic() +
  theme(
    text = element_text(size=24),
    legend.position = c(0.25,0.85),
    legend.key = element_blank(),
    legend.background = element_blank())
```

![Lineplot of ITT, male mice only](figures/itt-line-males-1.png)


```r
data.long %>%
  filter(Injection=="GFP") %>%
  group_by(Sex,Diet,Time) %>%
  summarize(Mean = mean(Glucose,na.rm=T),
            Error = se(Glucose)) %>%
  ggplot(aes(y=Mean,
             x=Time,
             ymin=Mean-Error,
             ymax=Mean+Error,
             col=Diet)) +
  geom_line() +
  geom_errorbar() +
  expand_limits(y=0) +
  facet_grid(~Sex) +
  labs(y="Blood Glucose (mg/dL)",
       x="Insulin (min)",
       title="Control mice only, effect of KD")
```

![Lineplot of ITT, all GFP mice](figures/itt-line-gfp-1.png)

### ITT Statistics


```r
library(lme4)
library(broom)
library(lmerTest)
itt.lme <- lmer(Glucose~as.factor(Time)+Diet+Sex+Injection+(1|ID), data=data.long)
itt.lme.null <- lmer(Glucose~as.factor(Time)+Diet+Sex+(1|ID), data=data.long)

summary(itt.lme)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: Glucose ~ as.factor(Time) + Diet + Sex + Injection + (1 | ID)
##    Data: data.long
## 
## REML criterion at convergence: 5712
## 
## Scaled residuals: 
##    Min     1Q Median     3Q    Max 
## -3.205 -0.674 -0.060  0.583  4.256 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  ID       (Intercept) 608      24.7    
##  Residual             728      27.0    
## Number of obs: 598, groups:  ID, 67
## 
## Fixed effects:
##                    Estimate Std. Error      df t value Pr(>|t|)    
## (Intercept)         136.035      7.199  94.422   18.90  < 2e-16 ***
## as.factor(Time)15   -30.284      4.662 523.097   -6.50  1.9e-10 ***
## as.factor(Time)30   -60.373      4.662 523.097  -12.95  < 2e-16 ***
## as.factor(Time)45   -66.970      4.662 523.097  -14.37  < 2e-16 ***
## as.factor(Time)60   -60.070      4.683 523.479  -12.83  < 2e-16 ***
## as.factor(Time)75   -42.282      4.683 523.479   -9.03  < 2e-16 ***
## as.factor(Time)90   -20.691      4.683 523.479   -4.42  1.2e-05 ***
## as.factor(Time)105   -6.858      4.683 523.479   -1.46     0.14    
## as.factor(Time)120    2.793      4.683 523.479    0.60     0.55    
## DietKeto             46.853      6.448  63.071    7.27  6.8e-10 ***
## SexM                  6.900      6.459  63.047    1.07     0.29    
## InjectionCre          0.515      6.437  63.045    0.08     0.94    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) a.(T)15 a.(T)3 a.(T)4 a.(T)6 a.(T)7 a.(T)9 a.(T)10
## as.fct(T)15 -0.324                                                   
## as.fct(T)30 -0.324  0.500                                            
## as.fct(T)45 -0.324  0.500   0.500                                    
## as.fct(T)60 -0.321  0.498   0.498  0.498                             
## as.fct(T)75 -0.321  0.498   0.498  0.498  0.497                      
## as.fct(T)90 -0.321  0.498   0.498  0.498  0.497  0.497               
## as.fc(T)105 -0.321  0.498   0.498  0.498  0.497  0.497  0.497        
## as.fc(T)120 -0.321  0.498   0.498  0.498  0.497  0.497  0.497  0.497 
## DietKeto    -0.415  0.000   0.000  0.000 -0.002 -0.002 -0.002 -0.002 
## SexM        -0.471  0.000   0.000  0.000  0.002  0.002  0.002  0.002 
## InjectionCr -0.414  0.000   0.000  0.000 -0.002 -0.002 -0.002 -0.002 
##             a.(T)12 DietKt SexM  
## as.fct(T)15                      
## as.fct(T)30                      
## as.fct(T)45                      
## as.fct(T)60                      
## as.fct(T)75                      
## as.fct(T)90                      
## as.fc(T)105                      
## as.fc(T)120                      
## DietKeto    -0.002               
## SexM         0.002  -0.041       
## InjectionCr -0.002  -0.072 -0.012
```

```r
anova(itt.lme) %>% 
  tidy %>% 
  kable(caption="Type III Analysis of Variance Table with Satterthwaite's method for ITT mixed linear model.")
```



Table: Type III Analysis of Variance Table with Satterthwaite's method for ITT mixed linear model.

term                   sumsq     meansq   NumDF   DenDF   statistic   p.value
----------------  ----------  ---------  ------  ------  ----------  --------
as.factor(Time)    393473.26   49184.16       8   523.3      67.555     0.000
Diet                38440.17   38440.17       1    63.1      52.798     0.000
Sex                   830.88     830.88       1    63.0       1.141     0.289
Injection               4.66       4.66       1    63.0       0.006     0.936

```r
anova(itt.lme,itt.lme.null) %>% 
  tidy %>%
  kable(caption="Chi-squared test for effects of AAV injection on ITT.")
```



Table: Chi-squared test for effects of AAV injection on ITT.

term            df    AIC    BIC   logLik   deviance   statistic   Chi.Df   p.value
-------------  ---  -----  -----  -------  ---------  ----------  -------  --------
itt.lme.null    13   5795   5852    -2884       5769          NA       NA        NA
itt.lme         14   5797   5858    -2884       5769       0.007        1     0.935

## Normalized ITT


```r
data.long.norm <-
  data.long %>%
  group_by(ID) %>%
  mutate(Glucose.norm = Glucose/Glucose[Time==0]*100)

data.long.norm %>%
  group_by(Sex,Diet,Injection,Time) %>%
  summarize(Mean = mean(Glucose.norm,na.rm=T),
            Error = se(Glucose.norm)) %>%
  ggplot(aes(y=Mean,
             x=Time,
             ymin=Mean-Error,
             ymax=Mean+Error,
             col=Injection)) +
  geom_line() +
  geom_errorbar() +
  expand_limits(y=0) +
  facet_grid(Sex~Diet) +
  labs(y="Blood Glucose (% of Basal)",
       x="Insulin (min)")
```

![Dotplot of ITT, normalized to fasting glucose](figures/itt-norm-line-all-1.png)


```r
data.long.norm %>%
  filter(Injection=="GFP") %>%
  group_by(Sex,Diet,Injection,Time) %>%
  summarize(Mean = mean(Glucose.norm,na.rm=T),
            Error = se(Glucose.norm)) %>%
  ggplot(aes(y=Mean,
             x=Time,
             ymin=Mean-Error,
             ymax=Mean+Error,
             col=Diet)) +
  geom_line() +
  geom_errorbar() +
  expand_limits(y=0) +
  facet_grid(~Sex) +
  labs(y="Blood Glucose (% of Basal)",
       x="Insulin (min)")
```

![Lineplot of ITT from GFP injected mice only, normalized to fasting glucose](figures/itt-norm-line-gfp-1.png)

## Fasting Glucose


```r
ggplot(data.long.norm %>% filter(Time==0),
       aes(y=Glucose,x=Diet,
           col=Injection)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  expand_limits(y=0) +
  facet_grid(.~Sex) +
  labs(y="Fasting Glucose (mg/dL)",
       title="Fasting Glucose")
```

![Boxplot of glucose levels after 6h fast](figures/fasting-glucose-boxplot-1.png)


```r
data.long %>% 
  filter(Time==0) %>%
  group_by(Injection,Sex,Diet) %>%
  summarize(Mean = mean(Glucose,na.rm=T),
            Error = se(Glucose)) %>%
  ggplot(aes(y=Mean,
             x=Diet,
             fill=Injection)) +
  geom_bar(stat='identity', position='dodge', width=0.75) +
  geom_errorbar(position=position_dodge(width=0.75),aes(group=Injection,
                                                        ymin=Mean-Error,
                                                        ymax=Mean+Error), width=0.5) +
  expand_limits(y=0) +
  facet_grid(.~Sex) +
  labs(y="Fasting Glucose (mg/dL)",
       title="Fasting Glucose")
```

![Barplot of glucose levels after 6h fast](figures/fasting-glucose-barplot-1.png)

### Fasting Glucose Statistics


```r
lm(Glucose ~ Diet + Sex + Injection, data = data.long %>% filter(Time==0)) %>%
  summary %>%
  tidy %>%
  kable(caption="Linear model for effects on fasting glucose levels.")
```



Table: Linear model for effects on fasting glucose levels.

term            estimate   std.error   statistic   p.value
-------------  ---------  ----------  ----------  --------
(Intercept)        138.2        6.16       22.46     0.000
DietKeto            28.3        6.11        4.63     0.000
SexM                11.0        6.12        1.80     0.077
InjectionCre        10.9        6.10        1.78     0.080

## Area Under the Curve


```r
data.long %>%
  group_by(ID, Diet, Sex, Injection) %>%
  summarize(AUC = sum(Glucose,na.rm=T)) %>%
ggplot(aes(y=AUC,x=Diet,
           col=Injection)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  expand_limits(y=0) +
  facet_grid(.~Sex) +
  labs(y="Glucose (mg/dL x minutes)",
       title="Area Under Curve")
```

![Boxplot of the area under the curves for insulin tolerance tests.](figures/auc-boxplot-1.png)


```r
data.long %>%
  group_by(ID, Diet, Sex, Injection) %>%
  summarize(AUC = sum(Glucose)) %>%
    group_by(Diet, Sex, Injection) %>%
  summarize(Mean = mean(AUC,na.rm=T),
            Error = se(AUC)) %>%
  ggplot(aes(y=Mean,
             x=Diet,
             fill=Injection)) +
  geom_bar(stat='identity', position='dodge', width=0.75) +
  geom_errorbar(position=position_dodge(width=0.75),aes(group=Injection,
                                                        ymin=Mean-Error,
                                                        ymax=Mean+Error), width=0.5) +
  expand_limits(y=0) +
  facet_grid(.~Sex) +
  labs(y="Glucose (mg/dL x minutes)",
       title="Area Under Curve")
```

![Barplot of the area under the curves for insulin tolerance tests.](figures/auc-barplot-1.png)

### AUC Statistics


```r
lm(AUC ~ Diet + Sex + Injection, 
   data = data.long %>% 
     group_by(ID, Diet, Sex, Injection) %>%
     summarize(AUC = sum(Glucose))) %>% 
  summary %>%
  tidy %>%
  kable(caption="Linear model for effects on ITT area under curve.")
```



Table: Linear model for effects on ITT area under curve.

term            estimate   std.error   statistic   p.value
-------------  ---------  ----------  ----------  --------
(Intercept)        447.2         114       3.914     0.000
DietKeto           216.8         111       1.949     0.054
SexM                26.0         112       0.233     0.816
InjectionCre        41.5         111       0.374     0.709

# Interpretation

A brief summary of what the interpretation of these results were

# Session Information


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS High Sierra 10.13.6
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] lmerTest_3.1-0 broom_0.5.2    lme4_1.1-21    Matrix_1.2-14 
## [5] ggplot2_3.1.0  readxl_1.3.1   dplyr_0.8.1    tidyr_0.8.2   
## [9] knitr_1.20    
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.1          nloptr_1.2.1        pillar_1.4.1       
##  [4] compiler_3.5.1      cellranger_1.1.0    RColorBrewer_1.1-2 
##  [7] highr_0.7           plyr_1.8.4          tools_3.5.1        
## [10] boot_1.3-20         digest_0.6.16       nlme_3.1-137       
## [13] lattice_0.20-35     evaluate_0.11       tibble_2.1.3       
## [16] gtable_0.2.0        pkgconfig_2.0.2     rlang_0.3.4        
## [19] yaml_2.2.0          withr_2.1.2         stringr_1.3.1      
## [22] generics_0.0.2      rprojroot_1.3-2     grid_3.5.1         
## [25] tidyselect_0.2.5    glue_1.3.1          R6_2.4.0           
## [28] rmarkdown_1.10      minqa_1.2.4         purrr_0.2.5        
## [31] reshape2_1.4.3      magrittr_1.5        MASS_7.3-50        
## [34] splines_3.5.1       backports_1.1.2     scales_1.0.0       
## [37] htmltools_0.3.6     assertthat_0.2.1    colorspace_1.3-2   
## [40] numDeriv_2016.8-1.1 labeling_0.3        stringi_1.2.4      
## [43] lazyeval_0.2.1      munsell_0.5.0       crayon_1.3.4
```

