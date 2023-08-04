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



This script was most recently run on Fri Aug  4 10:41:12 2023 and can be found in /Users/davebrid/Documents/GitHub/TissueSpecificTscKnockouts/Mouse Data/AJ Ketogenic Diet.

# Purpose


# Experimental Details

Beta hydroxybutyrate was injected at 1g/kg into fed mice around 14:00 PM, ketone levels were assessed via tail vein injection using strips.

# Raw Data



```r
library(readr)

raw_data_file <- "AJ Ketone Tolerance Data.csv"
kd.mouseids <- scan('KD Mouse IDs.txt')
library(readr) 
measurement.data <- read_csv(raw_data_file)

library(lubridate)
measurement.data$parsed.experiment.date <- ymd(measurement.data$experiment.date)

genotype.problems <- c() #animals where we are unsure about the genotype/treatment

ketone.data <- 
  filter(measurement.data, assay=='Plasma Glucose') %>%
  mutate(Diet = as.factor(if_else(animal.id %in% kd.mouseids, "Ketogenic Diet", "Control Diet"))) %>% 
  separate(values, sep=',', into=paste0("t",seq(0,120, by=15))) %>%
  select(age,Genotype,t0:t120,Diet,animal.id,MouseID) %>%
  mutate(ITT = if_else(is.na(t15), "No","Yes")) %>%
  mutate(FK = as.numeric(t0)/10) %>%
  mutate_at(.vars=3:11, .funs=funs(as.numeric(.)/10))

ketone.data <- 
  ketone.data %>% 
  mutate(AUC = rowSums(across(t0:t90)),
         AOC = rowSums(across(t0:t90))-(t0*7))

ketone.data.norm <-
  ketone.data %>%
  filter(ITT=="Yes") %>%
  filter(!(MouseID %in% genotype.problems)) %>% #removed animals with incorrect genotype
  mutate(t15=t15/t0*100,
         t30=t30/t0*100,
         t45=t45/t0*100,
         t60=t60/t0*100,
         t75=t75/t0*100,
         t90=t90/t0*100,
         t105=t105/t0*100,
         t120=t120/t0*100,
         t0=t0/t0*100) 

ketone.data.norm.abs <-
  ketone.data %>%
  filter(ITT=="Yes") %>%
  filter(!(MouseID %in% genotype.problems)) %>% #removed animals with incorrect genotype
  mutate(t15=t15-t0,
         t30=t30-t0,
         t45=t45-t0,
         t60=t60-t0,
         t75=t75-t0,
         t90=t90-t0,
         t105=t105-t0,
         t120=t120-t0,
         t0=t0-t0) 

ketone.data.norm.pk <-
  ketone.data %>%
  filter(ITT=="Yes") %>%
  filter(!(MouseID %in% genotype.problems)) %>% #removed animals with incorrect genotype
  mutate(t0=t0/t15*100,
         t30=t30/t15*100,
         t45=t45/t15*100,
         t60=t60/t15*100,
         t75=t75/t15*100,
         t90=t90/t15*100,
         t105=t105/t15*100,
         t120=t120/t15*100,
         t15=t15/t15*100) 



kable(ketone.data %>% 
        filter(ITT=="Yes") %>% 
        select(MouseID, Diet, FK,age) %>%
        arrange(age,MouseID, Diet),caption="Animals who underwent an KTT")
```



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

```r
ketone.data %>%
  filter(ITT=="Yes") %>% 
  group_by(Diet,age) %>%
  count %>%
  kable(caption="Animals evaluated by KTT in each group")
```



Table: Animals evaluated by KTT in each group

|Diet           | age|  n|
|:--------------|---:|--:|
|Control Diet   |  72|  3|
|Control Diet   |  92|  3|
|Ketogenic Diet |  72|  9|
|Ketogenic Diet |  92|  8|



# Analysis


```r
ketone.summary <-
  ketone.data %>%
  group_by(Diet,age) %>%
  summarize(FK.mean = mean(FK, na.rm=T),
            FK.se = se(FK),
            N = length(FK))
```

## Fed Ketone Levels


```r
library(ggplot2)
ggplot(filter(ketone.data, ITT=="Yes"),
            aes(y=FK, x=Diet)) +
  geom_boxplot() +
  facet_grid(.~age) +
  labs(title="Fed Ketones",
       y="Ketones (mg/dL)",
       x="Diet")
```

![Fed ketones, by treatment](figures/fasting-ketone-boxplot-itt-1.png)



```r
ggplot(ketone.summary,
       aes(y=FK.mean,
           ymin = FK.mean - FK.se,
           ymax = FK.mean + FK.se,
           x = Diet)) +
  facet_grid(.~age) +
  geom_bar(stat='identity', position='dodge', width=0.75) +
  geom_errorbar(position=position_dodge(width=0.75), width=0.5) +
  labs(y='Fed Blood Ketones (mg/dL)',
       x='')
```

![Fed ketones](figures/fasting-ketone-barplot-itt-ko-1.png)

### Fed ketone Statistics


```r
library(broom)

ketone.data %>% 
  filter(age=='92') %>%
  group_by(Diet) %>%
  summarize_at(.vars='t0',.funs=list(mean=~mean(.),
                             se=~se(.),
                             sd=~sd(.),
                             n=~length(.),
                             shapiro=~shapiro.test(.)$p.value)) %>%
  mutate(Pct.change=(mean-mean[Diet=='Control Diet'])/mean[Diet=='Control Diet']*100) -> fed.ketone.summary
  
kable(fed.ketone.summary, caption="Summary statistics for fed ketone levels")
```



Table: Summary statistics for fed ketone levels

|Diet           |  mean|    se|    sd|  n| shapiro| Pct.change|
|:--------------|-----:|-----:|-----:|--:|-------:|----------:|
|Control Diet   | 0.433| 0.033| 0.058|  3|   0.000|        0.0|
|Ketogenic Diet | 0.750| 0.033| 0.093|  8|   0.522|       73.1|

```r
wilcox.test(t0~Diet,data=ketone.data %>% 
  filter(age=='92')) %>% 
  tidy %>%
  kable(caption="Mann Whitney test for modifying effects of effects of diet on fed ketones at 3 weeks")
```



Table: Mann Whitney test for modifying effects of effects of diet on fed ketones at 3 weeks

| statistic| p.value|method                                            |alternative |
|---------:|-------:|:-------------------------------------------------|:-----------|
|         0|   0.017|Wilcoxon rank sum test with continuity correction |two.sided   |

```r
fk.aov <- aov(FK~Diet*age, data=ketone.data)
tidy(fk.aov) %>% kable(caption="ANOVA for modifying effects of effects of fed ketones")
```



Table: ANOVA for modifying effects of effects of fed ketones

|term      | df| sumsq| meansq| statistic| p.value|
|:---------|--:|-----:|------:|---------:|-------:|
|Diet      |  1| 0.596|  0.596|     1.739|   0.203|
|age       |  1| 0.028|  0.028|     0.081|   0.779|
|Diet:age  |  1| 0.010|  0.010|     0.029|   0.867|
|Residuals | 19| 6.516|  0.343|        NA|      NA|

## Ketone Tolerance Test


```r
time <- seq(0,120,by=15)
ktt.data.long <- 
  ketone.data %>%
  filter(ITT=="Yes") %>%
  group_by(Diet,animal.id) %>%
  gather(`t0`:`t120`,key='TimeStamp', value="ketone") %>%
  separate(TimeStamp, sep="t", into=c("Letter","Time")) %>%
  mutate(ketone = as.numeric(ketone)) %>%
  mutate(Time = as.numeric(Time))

ktt.data.long.norm <- 
  ketone.data.norm %>%
  filter(ITT=="Yes") %>%
  group_by(animal.id,Diet) %>%
  gather(`t0`:`t120`,key='TimeStamp', value="ketone") %>%
  separate(TimeStamp, sep="t", into=c("Letter","Time")) %>%
  mutate(ketone = as.numeric(ketone)) %>%
  mutate(Time = as.numeric(Time))

ktt.data.long.norm.pk <- 
  ketone.data.norm.pk %>%
  filter(ITT=="Yes") %>%
  group_by(animal.id,Diet) %>%
  gather(`t0`:`t120`,key='TimeStamp', value="ketone") %>%
  separate(TimeStamp, sep="t", into=c("Letter","Time")) %>%
  mutate(ketone = as.numeric(ketone)) %>%
  mutate(Time = as.numeric(Time))

ktt.data.long.norm.abs <- 
  ketone.data.norm.abs %>%
  filter(ITT=="Yes") %>%
  group_by(animal.id,Diet) %>%
  gather(`t0`:`t120`,key='TimeStamp', value="ketone") %>%
  separate(TimeStamp, sep="t", into=c("Letter","Time")) %>%
  mutate(ketone = as.numeric(ketone)) %>%
  mutate(Time = as.numeric(Time))

ggplot(ktt.data.long,
            aes(y=ketone, x=Time, col=Diet)) +
  geom_point() +
  facet_grid(.~age) +
  geom_smooth(method="loess") +
  labs(title="Ketone Tolerance Test",
       subtitle="4 Days of KD",
       y="Blood Ketones (mg/dL)") 
```

![Ketone tolerance tests by treatments](figures/ktt-dotplot-1.png)

```r
ggplot(ktt.data.long,
            aes(y=ketone, x=Time, col=Diet)) +
  geom_point() +
  geom_smooth(method="loess") +
  labs(title="BHB Tolerance Test",
       y="Blood Ketones (mg/dL)") 
```

![Ketone tolerance tests by treatments](figures/ktt-dotplot-2.png)

```r
ggplot(ktt.data.long.norm,
            aes(y=ketone, x=Time, col=Diet)) +
  geom_point() +
  facet_grid(.~age) +
  geom_smooth(method="loess") +
  labs(title="Ketone Tolerance Test",
       subtitle="4 Days of KD",
       y="Blood Ketones (% of Initial)") 
```

![Ketone tolerance tests by treatments](figures/ktt-dotplot-3.png)

```r
ggplot(ktt.data.long.norm.pk,
            aes(y=ketone, x=Time, col=Diet)) +
  geom_point() +
  facet_grid(.~age) +
  geom_smooth(method="loess") +
  labs(title="Ketone Tolerance Test",
       subtitle="4 Days of KD",
       y="Blood Ketones (% of Peak)") 
```

![Ketone tolerance tests by treatments](figures/ktt-dotplot-4.png)


```r
ktt.data.long %>%
  group_by(Time,Diet,age) %>%
  summarize(Average = mean(ketone),
            Error = se(ketone)) %>%
ggplot(aes(y=Average,
           ymax=Average+Error,
           ymin=Average-Error,
           x=Time, col=Diet)) +
  geom_point() +
  geom_line() +
  facet_grid(.~age) +
  geom_errorbar() +
  labs(title="BHB Tolerance Test",
       subtitle="4 Days of KD",
       y="Blood Ketones (mg/dL)") 
```

![](figures/ktt-lineplot-1.png)<!-- -->

```r
ktt.data.long %>%
  group_by(Time,Diet,age) %>%
  summarize(Average = mean(ketone),
            Error = se(ketone)) %>%
  filter(age=='92'&Time<100) %>%
ggplot(aes(y=Average,
           ymax=Average+Error,
           ymin=Average-Error,
           x=Time, lty=Diet)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  labs(title="BHB Tolerance Test",
       y="Blood Ketone Bodies (mg/dL)") +
  theme_classic() +
  theme(text=element_text(size=16),
        legend.position=c(0.9,0.8))
```

![](figures/ktt-lineplot-2.png)<!-- -->

```r
ktt.data.long.norm.abs %>%
  group_by(Time,Diet,age) %>%
  summarize(Average = mean(ketone),
            Error = se(ketone)) %>%
  filter(age=='92'&Time<100) %>%
ggplot(aes(y=Average,
           ymax=Average+Error,
           ymin=Average-Error,
           x=Time, lty=Diet)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  labs(title="BHB Tolerance Test - Baseline Subtracted",
       y="Blood Ketone Bodies (mg/dL)") +
  theme_classic() +
  theme(text=element_text(size=16),
        legend.position=c(0.9,0.8))
```

![](figures/ktt-lineplot-3.png)<!-- -->

## Area Under the Curve


```r
ketone.data %>%
  group_by(Diet,age) %>%
  summarize(mean=mean(AUC),
            se=se(AUC)) -> ktt.auc.summary

ktt.auc.summary %>%
  ggplot(aes(y=mean,
             ymin=mean-se,
             ymax=mean+se,
             x=Diet)) +
  geom_bar(stat='identity') +
  geom_errorbar(width=0.5) +
  facet_grid(~age) +
  labs(y="Area Under the Curve")
```

![](figures/ktt-auc-1.png)<!-- -->

## Normalized Area Under the Curve


```r
ketone.data %>%
  group_by(Diet,age) %>%
  summarize(mean=mean(AOC),
            se=se(AOC)) -> ktt.auc.summary

ktt.auc.summary %>%
  ggplot(aes(y=mean,
             ymin=mean-se,
             ymax=mean+se,
             x=Diet)) +
  geom_bar(stat='identity') +
  geom_errorbar(width=0.5) +
  facet_grid(~age) +
  labs(y="Area Under the Curve (Baseline Subtracted)")
```

![](figures/ktt-auc-norm-1.png)<!-- -->

```r
ktt.auc.summary %>%
  filter(age=='92') %>%
  ggplot(aes(y=mean,
             ymin=mean-se,
             ymax=mean+se,
             x=Diet)) +
  geom_bar(stat='identity') +
  geom_errorbar(width=0.5) +
  labs(title="BHB Tolerance Test", 
       subtitle="Baseline Subtracted",
       x="",
       y="Area Under the Curve") +
  theme_classic()+
  theme(text=element_text(size=16))
```

![](figures/ktt-auc-norm-2.png)<!-- -->

### Stats for Normalized Area Under the Curve


```r
ketone.data %>% 
  filter(age=='92') %>%
  group_by(Diet) %>%
  summarize_at(.vars='AOC',.funs=list(mean=~mean(.),
                             se=~se(.),
                             sd=~sd(.),
                             n=~length(.),
                             shapiro=~shapiro.test(.)$p.value)) %>%
  mutate(Pct.change=(mean-mean[Diet=='Control Diet'])/mean[Diet=='Control Diet']*100) -> fed.ketone.summary
  
kable(fed.ketone.summary, caption="Summary statistics for fed ketone levels")
```



Table: Summary statistics for fed ketone levels

|Diet           | mean|    se|   sd|  n| shapiro| Pct.change|
|:--------------|----:|-----:|----:|--:|-------:|----------:|
|Control Diet   |  7.6| 0.586| 1.01|  3|   0.672|        0.0|
|Ketogenic Diet | 10.7| 1.817| 5.14|  8|   0.083|       41.3|

```r
library(car)
leveneTest(AOC~Diet,data=ketone.data %>% 
  filter(age=='92')) %>% 
  tidy %>%
  kable(caption="Levene's test for modifying effects of effects of diet on the baseline subtracted area over the curve")
```



Table: Levene's test for modifying effects of effects of diet on the baseline subtracted area over the curve

| statistic| p.value| df| df.residual|
|---------:|-------:|--:|-----------:|
|      1.74|    0.22|  1|           9|

```r
t.test(AOC~Diet,data=ketone.data %>% 
  filter(age=='92'),
  var.equal=T) %>% 
  tidy %>%
  kable(caption="Student's t test for modifying effects of effects of diet on the baseline subtracted area over the curve")
```



Table: Student's t test for modifying effects of effects of diet on the baseline subtracted area over the curve

| estimate| estimate1| estimate2| statistic| p.value| parameter| conf.low| conf.high|method            |alternative |
|--------:|---------:|---------:|---------:|-------:|---------:|--------:|---------:|:-----------------|:-----------|
|    -3.14|       7.6|      10.7|     -1.02|   0.336|         9|    -10.1|      3.84|Two Sample t-test |two.sided   |

```r
t.test(AOC~Diet,data=ketone.data %>% 
  filter(age=='92'),
  var.equal=T,
  alternative="less") %>% 
  tidy %>%
  kable(caption="Student's t test for modifying effects of effects of diet on the baseline subtracted area over the curve")
```



Table: Student's t test for modifying effects of effects of diet on the baseline subtracted area over the curve

| estimate| estimate1| estimate2| statistic| p.value| parameter| conf.low| conf.high|method            |alternative |
|--------:|---------:|---------:|---------:|-------:|---------:|--------:|---------:|:-----------------|:-----------|
|    -3.14|       7.6|      10.7|     -1.02|   0.168|         9|     -Inf|      2.52|Two Sample t-test |less        |


## Normalized to fasting ketone levels (percent change)


```r
ktt.data.long.norm %>%
  filter(!(animal.id %in% c('23224','23225'))) %>%
  group_by(Time,Diet,age) %>%
  summarize(Average = mean(ketone),
            Error = se(ketone)) %>%
ggplot(aes(y=Average,
           ymax=Average+Error,
           ymin=Average-Error,
           x=Time, col=Diet)) +
  geom_point() +
  geom_line() +
  facet_grid(.~age) +
  geom_errorbar() +
  labs(title="Ketone Tolerance Test",
       subtitle="4 Days of KD",
       y="Blood Ketones (% of initial)") 
```

![](figures/normalized-ktt-1.png)<!-- -->

```r
ktt.data.long.norm %>%
  filter(!(animal.id %in% c('23224','23225'))) %>%
  filter(age == '92') %>%
  group_by(Time,Diet,age) %>%
  summarize(Average = mean(ketone),
            Error = se(ketone)) %>%
ggplot(aes(y=Average,
           ymax=Average+Error,
           ymin=Average-Error,
           x=Time, col=Diet)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  scale_color_manual(values=color.scheme, labels=c('Control Diet','LCHF Diet'), name="") +
  labs(title="Ketone Tolerance Test",
       subtitle="3 Weeks of LCHF",
       y="Blood Ketones (% of Basal)") +
  xlim(0,95) +
  theme_classic() +
  theme(text = element_text(size=18),
        legend.position = c(0.75,0.75))
```

![](figures/normalized-ktt-2.png)<!-- -->


## Normalized to fasting ketone levels (absolute change)


```r
ktt.data.long.norm.abs %>%
  filter(!(animal.id %in% c('23224','23225'))) %>%
  group_by(Time,Diet,age) %>%
  summarize(Average = mean(ketone),
            Error = se(ketone)) %>%
ggplot(aes(y=Average,
           ymax=Average+Error,
           ymin=Average-Error,
           x=Time, col=Diet)) +
  geom_point() +
  geom_line() +
  facet_grid(.~age) +
  geom_errorbar() +
  labs(title="Ketone Tolerance Test",
       subtitle="4 Days of KD",
       y="Blood Ketones (% of initial)") 
```

![](figures/normalized-ktt-abs-1.png)<!-- -->

```r
ktt.data.long.norm.abs %>%
  filter(!(animal.id %in% c('23224','23225'))) %>%
  filter(age == '92') %>%
  group_by(Time,Diet,age) %>%
  summarize(Average = mean(ketone),
            Error = se(ketone)) %>%
ggplot(aes(y=Average,
           ymax=Average+Error,
           ymin=Average-Error,
           x=Time, col=Diet)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  scale_color_manual(values=color.scheme, labels=c('Control Diet','LCHF Diet'), name="") +
  labs(title="Ketone Tolerance Test",
       subtitle="3 Weeks of LCHF",
       y="Blood Ketones (Above Basal)") +
  xlim(0,95) +
  theme_classic() +
  theme(text = element_text(size=18),
        legend.position = c(0.75,0.75))
```

![](figures/normalized-ktt-abs-2.png)<!-- -->


## Normlized to the Peak of Ketone levels


```r
ktt.data.long.norm.pk %>%
  group_by(Time,Diet,age) %>%
  summarize(Average = mean(ketone),
            Error = se(ketone)) %>%
ggplot(aes(y=Average,
           ymax=Average+Error,
           ymin=Average-Error,
           x=Time, col=Diet)) +
  geom_point() +
  geom_line() +
  facet_grid(.~age) +
  geom_errorbar() +
  labs(title="Ketone Tolerance Test",
       subtitle="4 Days and 3 Weeks of KD",
       y="Blood Ketones (% of Peak)") 
```

![](figures/ktt-peak-1.png)<!-- -->

## KTT Statistics



### Linear Models


```r
library(lme4)
library(lmerTest)
ktt.lme <- lmer(ketone ~ as.factor(Time) + as.factor(age) + Diet + age:Diet + (1|animal.id),
                data=ktt.data.long)

ktt.lme.92d.null <- lmer(ketone ~ as.factor(Time) + (1|animal.id),
                data=ktt.data.long %>% filter (age==92))

ktt.lme.92d <- lmer(ketone ~ as.factor(Time) + Diet + (1|animal.id),
                data=ktt.data.long %>% filter (age==92))

ktt.lme.92d.norm <- lmer(ketone ~ as.factor(Time) + Diet + (1|animal.id),
                data=ktt.data.long.norm.abs %>% filter (age==92))

ktt.lme.92d.norm.null <- lmer(ketone ~ as.factor(Time) + (1|animal.id),
                data=ktt.data.long.norm.abs %>% filter (age==92))

anova(ktt.lme.92d) %>% tidy %>% kable(caption="Mixed linear model for KTT at 3 weeks")
```



Table: Mixed linear model for KTT at 3 weeks

|term            |  sumsq| meansq| NumDF| DenDF| statistic| p.value|
|:---------------|------:|------:|-----:|-----:|---------:|-------:|
|as.factor(Time) | 42.058|  7.010|     6|    60|     32.36|   0.000|
|Diet            |  0.597|  0.597|     1|     9|      2.76|   0.131|

```r
summary(ktt.lme.92d)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: ketone ~ as.factor(Time) + Diet + (1 | animal.id)
##    Data: ktt.data.long %>% filter(age == 92)
## 
## REML criterion at convergence: 134
## 
## Scaled residuals: 
##    Min     1Q Median     3Q    Max 
## -2.422 -0.505 -0.146  0.636  2.403 
## 
## Random effects:
##  Groups    Name        Variance Std.Dev.
##  animal.id (Intercept) 0.432    0.657   
##  Residual              0.217    0.465   
## Number of obs: 77, groups:  animal.id, 11
## 
## Fixed effects:
##                    Estimate Std. Error     df t value Pr(>|t|)    
## (Intercept)           0.107      0.414 11.056    0.26     0.80    
## as.factor(Time)15     2.418      0.198 60.000   12.19  < 2e-16 ***
## as.factor(Time)30     2.127      0.198 60.000   10.72  1.4e-15 ***
## as.factor(Time)45     1.755      0.198 60.000    8.84  1.8e-12 ***
## as.factor(Time)60     1.355      0.198 60.000    6.83  5.0e-09 ***
## as.factor(Time)75     1.173      0.198 60.000    5.91  1.7e-07 ***
## as.factor(Time)90     1.055      0.198 60.000    5.31  1.7e-06 ***
## DietKetogenic Diet    0.765      0.461  9.000    1.66     0.13    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) a.(T)1 a.(T)3 a.(T)4 a.(T)6 a.(T)7 a.(T)9
## as.fct(T)15 -0.240                                          
## as.fct(T)30 -0.240  0.500                                   
## as.fct(T)45 -0.240  0.500  0.500                            
## as.fct(T)60 -0.240  0.500  0.500  0.500                     
## as.fct(T)75 -0.240  0.500  0.500  0.500  0.500              
## as.fct(T)90 -0.240  0.500  0.500  0.500  0.500  0.500       
## DietKtgncDt -0.810  0.000  0.000  0.000  0.000  0.000  0.000
```

```r
anova(ktt.lme.92d,ktt.lme.92d.null) %>% tidy %>% kable(caption="Chi squared test for non-adjusted KTT")
```



Table: Chi squared test for non-adjusted KTT

|term             | npar| AIC| BIC| logLik| deviance| statistic| df| p.value|
|:----------------|----:|---:|---:|------:|--------:|---------:|--:|-------:|
|ktt.lme.92d.null |    9| 143| 164|  -62.5|      125|        NA| NA|      NA|
|ktt.lme.92d      |   10| 142| 165|  -61.0|      122|      2.94|  1|   0.086|

```r
anova(ktt.lme.92d.norm) %>% tidy %>% kable(caption="Mixed linear model for baseline adjusted KTT at 3 weeks")
```



Table: Mixed linear model for baseline adjusted KTT at 3 weeks

|term            |  sumsq| meansq| NumDF| DenDF| statistic| p.value|
|:---------------|------:|------:|-----:|-----:|---------:|-------:|
|as.factor(Time) | 42.058|  7.010|     6|    60|     32.36|   0.000|
|Diet            |  0.224|  0.224|     1|     9|      1.03|   0.336|

```r
summary(ktt.lme.92d.norm)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: ketone ~ as.factor(Time) + Diet + (1 | animal.id)
##    Data: ktt.data.long.norm.abs %>% filter(age == 92)
## 
## REML criterion at convergence: 133
## 
## Scaled residuals: 
##    Min     1Q Median     3Q    Max 
## -2.414 -0.517 -0.141  0.635  2.410 
## 
## Random effects:
##  Groups    Name        Variance Std.Dev.
##  animal.id (Intercept) 0.393    0.627   
##  Residual              0.217    0.465   
## Number of obs: 77, groups:  animal.id, 11
## 
## Fixed effects:
##                    Estimate Std. Error     df t value Pr(>|t|)    
## (Intercept)          -0.326      0.398 11.253   -0.82     0.43    
## as.factor(Time)15     2.418      0.198 60.000   12.19  < 2e-16 ***
## as.factor(Time)30     2.127      0.198 60.000   10.72  1.4e-15 ***
## as.factor(Time)45     1.755      0.198 60.000    8.84  1.8e-12 ***
## as.factor(Time)60     1.355      0.198 60.000    6.83  5.0e-09 ***
## as.factor(Time)75     1.173      0.198 60.000    5.91  1.7e-07 ***
## as.factor(Time)90     1.055      0.198 60.000    5.31  1.7e-06 ***
## DietKetogenic Diet    0.448      0.441  9.000    1.02     0.34    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) a.(T)1 a.(T)3 a.(T)4 a.(T)6 a.(T)7 a.(T)9
## as.fct(T)15 -0.249                                          
## as.fct(T)30 -0.249  0.500                                   
## as.fct(T)45 -0.249  0.500  0.500                            
## as.fct(T)60 -0.249  0.500  0.500  0.500                     
## as.fct(T)75 -0.249  0.500  0.500  0.500  0.500              
## as.fct(T)90 -0.249  0.500  0.500  0.500  0.500  0.500       
## DietKtgncDt -0.806  0.000  0.000  0.000  0.000  0.000  0.000
```

```r
anova(ktt.lme.92d.norm,ktt.lme.92d.norm.null) %>% tidy %>% kable(caption="Chi squared test for baseline adjusted KTT")
```



Table: Chi squared test for baseline adjusted KTT

|term                  | npar| AIC| BIC| logLik| deviance| statistic| df| p.value|
|:---------------------|----:|---:|---:|------:|--------:|---------:|--:|-------:|
|ktt.lme.92d.norm.null |    9| 140| 161|  -61.1|      122|        NA| NA|      NA|
|ktt.lme.92d.norm      |   10| 141| 164|  -60.5|      121|       1.2|  1|   0.274|

# Session Information


```r
sessionInfo()
```

```
## R version 4.2.2 (2022-10-31)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur ... 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] lmerTest_3.1-3  lme4_1.1-33     Matrix_1.5-4.1  car_3.1-2      
##  [5] carData_3.0-5   broom_1.0.4     ggplot2_3.4.2   lubridate_1.9.2
##  [9] readr_2.1.4     dplyr_1.1.2     tidyr_1.3.0     knitr_1.43     
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.10         lattice_0.21-8      digest_0.6.31      
##  [4] utf8_1.2.3          R6_2.5.1            backports_1.4.1    
##  [7] evaluate_0.21       highr_0.10          pillar_1.9.0       
## [10] rlang_1.1.1         rstudioapi_0.14     minqa_1.2.5        
## [13] jquerylib_0.1.4     nloptr_2.0.3        rmarkdown_2.22     
## [16] labeling_0.4.2      splines_4.2.2       stringr_1.5.0      
## [19] bit_4.0.5           munsell_0.5.0       numDeriv_2016.8-1.1
## [22] compiler_4.2.2      xfun_0.39           pkgconfig_2.0.3    
## [25] mgcv_1.8-42         htmltools_0.5.5     tidyselect_1.2.0   
## [28] tibble_3.2.1        fansi_1.0.4         crayon_1.5.2       
## [31] tzdb_0.4.0          withr_2.5.0         MASS_7.3-60        
## [34] grid_4.2.2          nlme_3.1-162        jsonlite_1.8.5     
## [37] gtable_0.3.3        lifecycle_1.0.3     magrittr_2.0.3     
## [40] scales_1.2.1        cli_3.6.1           stringi_1.7.12     
## [43] vroom_1.6.3         cachem_1.0.8        farver_2.1.1       
## [46] bslib_0.4.2         generics_0.1.3      vctrs_0.6.2        
## [49] boot_1.3-28.1       tools_4.2.2         bit64_4.0.5        
## [52] glue_1.6.2          purrr_1.0.1         hms_1.1.3          
## [55] abind_1.4-5         parallel_4.2.2      fastmap_1.1.1      
## [58] yaml_2.3.7          timechange_0.2.0    colorspace_2.1-0   
## [61] sass_0.4.6
```

