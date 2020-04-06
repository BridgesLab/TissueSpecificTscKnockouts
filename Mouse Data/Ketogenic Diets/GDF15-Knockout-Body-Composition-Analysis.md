---
title: "Body Composition of GDF15 Knockout Mice on KD"
author: "Claire Gleason, JeAnna Redd and Dave Bridges"
date: "June 17, 2019"
output:
  html_document:
    highlight: tango
    keep_md: yes
    toc: yes
    fig_caption: yes
  pdf_document:
    highlight: tango
    keep_tex: yes
    toc: yes
    fig_caption: yes
---




This script can be found in /Users/davebrid/Documents/GitHub/TissueSpecificTscKnockouts/Mouse Data/Ketogenic Diets and was most recently run on Mon Apr  6 14:21:43 2020.


# Data Entry



This script pulled in a total of 1595 observations.  This includes the following number of animals in each treatment group.

# Enrolled Animals

This is for animals where wer have any body composition data.


Table: Animals in each group of this cohort

Sex   Diet   Genotype     n
----  -----  ---------  ---
M     NA     -/-         10
M     NA     +/+          9
F     NA     -/-         14
F     NA     +/+          8



Table: Animals in each group of this cohort that have been put on KD

Sex   Diet   Genotype     n
----  -----  ---------  ---
M     NA     -/-         10
M     NA     +/+          9
F     NA     -/-         14
F     NA     +/+          7



Table: Animals in each group of this cohort

MouseID   Sex   Diet   Genotype 
--------  ----  -----  ---------
8651      M     NA     -/-      
8652      M     NA     -/-      
8653      M     NA     -/-      
8661      M     NA     -/-      
8662      M     NA     -/-      
9153      M     NA     -/-      
9156      M     NA     -/-      
9185      M     NA     -/-      
9337      M     NA     -/-      
9400      M     NA     -/-      
8905      M     NA     +/+      
8903      M     NA     +/+      
8904      M     NA     +/+      
8901      M     NA     +/+      
9154      M     NA     +/+      
9184      M     NA     +/+      
9329      M     NA     +/+      
9401      M     NA     +/+      
9402      M     NA     +/+      
8655      F     NA     -/-      
8656      F     NA     -/-      
8659      F     NA     -/-      
8907      F     NA     -/-      
8908      F     NA     -/-      
9159      F     NA     -/-      
9186      F     NA     -/-      
9189      F     NA     -/-      
9190      F     NA     -/-      
9191      F     NA     -/-      
9336      F     NA     -/-      
9403      F     NA     -/-      
9404      F     NA     -/-      
9405      F     NA     -/-      
8654      F     NA     +/+      
8902      F     NA     +/+      
8906      F     NA     +/+      
8981      F     NA     +/+      
9158      F     NA     +/+      
9160      F     NA     +/+      
9188      F     NA     +/+      
9335      F     NA     +/+      

# Body Weight

![Scatter Plot of Body Weights](figures/gdf15-ko-body-weight-scatterplot-1.png)

![Line Plot of Individual Body Weights](figures/gdf15-ko-body-weight-individual-1.png)

![Line Plot of Body Weights](figures/gdf15-ko-body-weight-lineplot-1.png)

# Lean Mass

![Scatter Plot of Lean Mass](figures/gdf15-ko-lean-mass-scatterplot-1.png)

![Line Plot of Lean Mass](figures/gdf15-ko-lean-mass-lineplot-1.png)

# Fat Mass

![Scatter Plot of Fat Mass](figures/gdf15-ko-fat-mass-scatterplot-1.png)![Scatter Plot of Fat Mass](figures/gdf15-ko-fat-mass-scatterplot-2.png)

![Line Plot of Fat Mass](figures/gdf15-ko-fat-mass-lineplot-1.png)![Line Plot of Fat Mass](figures/gdf15-ko-fat-mass-lineplot-2.png)![Line Plot of Fat Mass](figures/gdf15-ko-fat-mass-lineplot-3.png)

## Fat Mass Linear Models


Table: Mixed linear models for fat mass changes on diet for males

                          Estimate   Std. Error     df   t value   Pr(>|t|)
-----------------------  ---------  -----------  -----  --------  ---------
(Intercept)                  1.702        0.215   30.2     7.899      0.000
Genotype+/+                  0.311        0.314   30.4     0.992      0.329
Diet.Weeks                   0.548        0.038   95.2    14.594      0.000
Genotype+/+:Diet.Weeks      -0.216        0.060   94.8    -3.583      0.001



Table: Mixed linear models for fat mass changes on diet for females

                          Estimate   Std. Error      df   t value   Pr(>|t|)
-----------------------  ---------  -----------  ------  --------  ---------
(Intercept)                  1.939        0.142    41.6    13.671      0.000
Genotype+/+                 -0.148        0.236    36.7    -0.629      0.533
Diet.Weeks                   0.211        0.035   100.3     6.072      0.000
Genotype+/+:Diet.Weeks      -0.110        0.057    98.0    -1.909      0.059



Table: Chi-squared test between models with and without the genotype term and interaction for male mice

                    Df   AIC   BIC   logLik   deviance   Chisq   Chi Df   Pr(>Chisq)
-----------------  ---  ----  ----  -------  ---------  ------  -------  -----------
male.fat.lm.null     4   221   232     -107        213      NA       NA           NA
male.fat.lm          6   213   229     -100        201    12.5        2        0.002



Table: Chi-squared test between models with and without the genotype term and interaction for female mice

                      Df   AIC   BIC   logLik   deviance   Chisq   Chi Df   Pr(>Chisq)
-------------------  ---  ----  ----  -------  ---------  ------  -------  -----------
female.fat.lm.null     4   197   208    -94.5        189      NA       NA           NA
female.fat.lm          6   194   211    -91.0        182    6.95        2        0.031

### Modifying Effect of Sex on Fat Mass Changes


Table: Mixed linear models for fat mass changes on diet including interaction variable

                               Estimate   Std. Error      df   t value   Pr(>|t|)
----------------------------  ---------  -----------  ------  --------  ---------
(Intercept)                       1.705        0.191    67.8     8.937      0.000
Genotype+/+                       0.308        0.278    68.0     1.107      0.272
Diet.Weeks                        0.547        0.036   193.2    15.134      0.000
SexF                              0.244        0.252    69.8     0.966      0.337
Genotype+/+:Diet.Weeks           -0.215        0.058   192.2    -3.706      0.000
Genotype+/+:SexF                 -0.466        0.392    66.6    -1.187      0.239
Diet.Weeks:SexF                  -0.341        0.051   195.1    -6.653      0.000
Genotype+/+:Diet.Weeks:SexF       0.109        0.083   193.0     1.306      0.193



Table: Chi-squared test between full models with and without the inclusion of sex

                   Df   AIC   BIC   logLik   deviance   Chisq   Chi Df   Pr(>Chisq)
----------------  ---  ----  ----  -------  ---------  ------  -------  -----------
sex.fat.lm.null     6   460   481     -224        448      NA       NA           NA
sex.fat.lm         10   406   440     -193        386    62.2        4        1e-12



Table: Chi-squared test betweenn models with and without the genotype term, accounting for the modifying effects of sex

                   Df   AIC   BIC   logLik   deviance   Chisq   Chi Df   Pr(>Chisq)
----------------  ---  ----  ----  -------  ---------  ------  -------  -----------
sex.fat.lm.geno     6   417   438     -203        405      NA       NA           NA
sex.fat.lm         10   406   440     -193        386    19.3        4      0.00068
