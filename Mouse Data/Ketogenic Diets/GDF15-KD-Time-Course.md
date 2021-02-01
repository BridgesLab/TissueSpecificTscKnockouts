---
title: "GDF15 Time Course ELISA"
author: "Molly C. Mulcahy"
date: "12/22/2020"
output:
  html_document:
    highlight: tango
    keep_md: yes
    number_sections: yes
    toc: yes
---




# Data Entry

```r
#Enters results from Dilution Factor = 5
my.assays.GDF15.data <- read_csv("GDF15 Time Course ELISA 1.csv",
                           col_types = cols(
                            Sample = col_character(),
                            Dilution = col_double(),
                             Well = col_factor(levels = NULL),
                            raw = col_double(),
                            Conc. = col_double()
                           ))
                   
mapping.data <-read_csv("GDF15 Time Course ELISA Mapping 1.csv",
                        col_types = cols(
                          well =	col_factor(levels = NULL),
                          ID = col_factor(levels = NULL),
                          diet	= col_factor(levels = NULL),
                          time	= col_double(),
                          study = col_factor(levels = NULL)
                        ))

full.data <- merge(mapping.data, my.assays.GDF15.data, by = "Well")%>%
  select(ID, Well, study, time, diet, study, Conc.,Raw)
non.ID<-c("blank", "pos.control")
```


# Analysis


```r
#make dataset for time course study
time.course.data<-full.data%>%
  filter(study == "time.course")%>%
  group_by(time)%>%
  summarize(average.GDF15 = mean(Conc.),
            error.GDF15 = se(Conc.))

#Plot of time course study- averaged
ggplot(time.course.data, aes(time, average.GDF15))+
  geom_point(col = "dodgerblue")+
  geom_errorbar(aes(ymin = average.GDF15 - error.GDF15, ymax = average.GDF15 + error.GDF15), col = "dodgerblue")+
  geom_line(col = "dodgerblue")+
  labs(title = "GDF15 in Ketogenic Diets", y= "GDF15 concentration (pg/mL)", x = "Time (days)")+
  expand_limits(x = 0, y = 0)+
  geom_hline(yintercept=294, linetype = "dotdash")
```

![](figures/kd-time-course-1.png)<!-- -->

```r
# time course study-Individual plots
time.course.ID.data<-full.data%>%
  filter(study=="time.course")
ggplot(time.course.ID.data, aes(time, `Conc.`, col = ID))+
  geom_point(aes(col = ID))+
  #geom_errorbar(aes(ymin = average.GDF15 - error.GDF15, ymax = average.GDF15 + error.GDF15), col = "dodgerblue")+
  geom_line(aes(col = ID))+
  labs(title = "GDF15 in Ketogenic Diets", y= "GDF15 concentration (pg/mL)", x = "Time (days)")
```

![](figures/kd-time-course-2.png)<!-- -->
