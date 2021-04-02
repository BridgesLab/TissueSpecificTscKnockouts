---
title: "Combined Barplots for aTSC Mammary Gland Studies"
author: "Dave Bridges and Noura El Habbal"
date: "August 29, 2020"
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

# Raw Data

Imported normalized counts from Salmon aligned data, analyzed via DESeq2


```r
library(readr) #loads the readr package
counts.file <- 'DESeq2 Normalized Counts.csv'
design.file <- 'Design Table.csv'

counts.table <- read_csv(counts.file) 

library(biomaRt)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
gene.mapping <- getBM(attributes=c('mgi_symbol','ensembl_gene_id'),
                      values=counts.table$X1,
                      filters='ensembl_gene_id',
                      mart=ensembl)

annotated.counts.table <- full_join(counts.table,gene.mapping, by=c('X1'='ensembl_gene_id'))

long.counts.table <- pivot_longer(annotated.counts.table,
                                  cols=starts_with('14'),
                                  names_to="Sample",
                                  values_to="TPM")

design.table <- read_csv(design.file)
long.merged.data <- 
  design.table %>%
  dplyr::select(names, Genotype,) %>%
  mutate(Genotype= relevel(as.factor(Genotype),ref="WT")) %>%
  full_join(long.counts.table, by=c('names'='Sample'))
```

These data can be found in **/Users/davebrid/Documents/GitHub/TissueSpecificTscKnockouts/RNAseq/Mammary Gland Adipocyte Tsc1 Knockout**.  The counts file is DESeq2 Normalized Counts.csv but the design file is Design Table.csv.

# Analysis

# Single Gene Graphs


```r
#use MGI symbols and add them to this list

human.mouse.mapping.table <- 'http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt'
mapping.table <- read_tsv(human.mouse.mapping.table) %>%
  dplyr::select(`Common Organism Name`,Symbol,`HomoloGene ID`) 
  
wide.mapping.table <- pivot_wider(mapping.table,
                                  names_from=`Common Organism Name`,
                                  values_from = Symbol,
                                  id_cols=`HomoloGene ID`) %>%
  rename("Mouse"='mouse, laboratory') %>%
  mutate(Mouse=as.factor(as.character(Mouse))) %>%
    mutate(human=as.factor(as.character(human)))


genes.of.interest <- c('Cxcl13','Cxcr5','Ccr6','Ccr7', 'Il12a', 'Fut7', 'Tnfsf14', 'Cr2', 'Azgp1','Aicda', 'Lef1', 'Cd19', 'Cd22', 'Siglecg','Foxp3', 'Cxcr5', 'Ntrk1', 'Sstr4', 'Pla2g2d', 'Pla2g5', 'Pla2g3','Pla2g4f','Pla2g6', 'Cyp4f18', 'Agtr2', 'Fabp3', 'Zan', 'Cntnap2', 'Foxp3', 'Pax5','Il10', 'Pcdha10','Cnr2', 'Kcna1', 'Azgp1', 'Lef1', 'Cd79a','Cd3d', 'Cd247', 'Cd4', 'Cd3g', 'Cd28', 'Cd3e', 'Ccr7', 'Lck', 'Blk', 'Itk', 'Syk', 'Myoz2', 'Mypn', 'Mybph', 'Myom2', 'Myh8', 'Myh7', 'Myh4', 'Myh1', 'Myl2','Myl1', 'MYPF', 'Acta1','Ryr1','Csrp3','Prkag3', 'Tnnt2','Tnnt3', 'Acaca', 'Srebf1', 'Acly', 'Fasn', 'Rps6kb1', 'Eif4ebp1', 'Acat1','Akt1', 'Rps14', 'Plin1', 'Plin2',  'Fabp4', 'Dgat1', 'Dgat2', 'Lpl', 'Pnpla2', 'Gpat3', 'Acsl1', 'Acsl3', 'Acsl4', 'Acsl5', 'Acsl6', 'Elovl1', 'Elovl2', 'Elovl3', 'Elovl4', 'Elovl5', 'Elovl6', 'Elovl7', 'Scd1', 'Scd2', 'Scd3', 'Scd4', 'Fads1', 'Fads2','Fads2b','Fads3','Fads6', 'Mfsd2a', 'Mfsd6', 'Mfsd8', 'Mfs10', 'Mfsd14a', 'Mfsd14b', 'Cyp2d34', 'Cyp4f18', 'Cyp1a2', 'Cyp2c50', 'Cyp4a12a', 'Cyp4a12b', 'Cyp2c50') #you can choose any gene name  but by default it show the sig genes , previously this code had c('Fasn','Wap','Csn2','Insr') UNCERTAIN about CNR3 to find homologous if it is pcdha10 or not, used Pcdha10 for mouse which has different homologene ID number that that of humans; MYPF does not exist as a gene (possible typo? when indeed it is MYLF?); 
#from the DESEQ model
deseq.results.file <-  'DESeq2 Results.csv'
deseq.results <- read_csv(deseq.results.file)
sig.genes <- filter(deseq.results, padj<0.05) %>% pull(symbol)

genes.to.plot <- c(genes.of.interest)#,sig.genes)
```


```r
library(ggplot2)
for (gene in genes.to.plot){
  barplot <- long.merged.data %>%
    filter(mgi_symbol==gene) %>%
    group_by(Genotype) %>%
    summarize(Mean = mean(TPM),
              Error =se(TPM)) %>%
    ggplot(aes(y=Mean,x=Genotype,fill=Genotype)) +
    geom_bar(stat='identity') +
    geom_errorbar(aes(ymin=Mean-Error,
                      ymax=Mean+Error),
                  width=0.5) +
    theme_bw() +
    theme(text=element_text(size=18)) +
    labs(y="Normalized Counts (TPM)",
         x="",
         title=gene)
  
  ggsave(barplot,
         device='pdf',
         filename=c(paste('figures-noura/barplots/',paste(gene,'-barplot.pdf', sep=''), sep=''))) 
  ggsave(barplot, 
         device='png',
         filename=c(paste('figures-noura/barplots/',paste(gene,'-barplot.png', sep=''), sep='')))
}
```


```r
for (gene in genes.to.plot){
  barplot <- long.merged.data %>%
    dplyr::filter(mgi_symbol==gene) %>%
    ggplot(aes(y=TPM,x=Genotype,col=Genotype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position=position_jitterdodge(jitter.width=0.25)) +
    theme_bw() +
    theme(text=element_text(size=18)) +
    labs(y="Normalized Counts (TPM)",
         x="",
         title=gene) +
    expand_limits(y=0)
  
  ggsave(barplot,
         device='pdf',
         filename=c(paste('figures-noura/boxplots/',paste(gene,'-boxplot.pdf', sep=''), sep=''))) 
  ggsave(barplot, 
         device='png',
         filename=c(paste('figures-noura/boxplots/',paste(gene,'-boxplot.png', sep=''), sep='')))
}
```

# Grouped Gene Sets


```r
gene.sets <- tibble(name="Fatty Acid Synthesis", genes = list('Fasn','Acly','Acaca'))
gene.sets <- 
  gene.sets %>% 
    add_row(name="Milk Proteins", genes=list('Wap','Csn1s1','Csn1s2a','Csn1s2b','Csn2','Csn3','Lalba','Ltf')) %>%
  add_row(name="T and B Cell Markers", genes=list('Cd3e','Cd3d','Cd3g','Ptprc','Cd19')) %>% #CD3, CD45(Ptprc), CD19
    add_row(name="B Cell Markers", genes=list('Leuf1','Ptprc','Cd19','Sdc1','Tnfrsf13b','Spn','Cd79a')) %>% #see https://www.novusbio.com/antibody-news/how-to-identify-b-cell-subsets-using-flow-cytometry CD45(Leufr), CD45R/B220(Ptprc), CD19, CD138(Sdc1), TACI(Tnfrsf13b), CD43(Spn)
  add_row(name="Activated B Cells", genes=list('Cd19','Il2ra','Tnfrsf8')) %>%
  add_row(name="Plasma Cells", genes=list('Cd27','Cd38','Cd78','Sdc1','Slamf7','Il6')) %>%  
    add_row(name="Memory B Cells", genes=list('Cd19','Cd27','Ms4a1','Cd40','Cd80','	Pdcd1lg2','Cxcr3','Cxcr5','Cxcr6')) %>%  
   add_row(name="T Cell Markers", genes=list('Cd8a','Cd4','Il2ra','Foxp3','Il7r')) %>% #CD8 (Killer), CD4 (helper), CD25(Il2ra)/FoxP3,CD127(Il7r) Tregs)
  add_row(name="Macrophage Markers", genes=list('Fcgr1','Ly6g','Itgam','Itgax','Csf1r')) %>% #CD64=Fcgr1,CD11B=Itgam,CD11C=Itgax, CD115=Csf1r 
  add_row(name="Phospholipases",
          genes=list(gene.mapping[grepl('Pla2g',gene.mapping$mgi_symbol),]$mgi_symbol)) %>%
 add_row(name="New Targets",
         genes=list('Thrsp','Agpat6','Plin2','Aox1','Aox2','Aox3','Aox4','Btn1a1','Btn1a2')) %>%
  add_row(name="Xanthine Oxidoreductases",
          genes=list('Aox1','Aox2','Aox3','Aox4','Xdh')) %>%
  add_row(name="DHA Receptors",
          genes=list('Gpr18','Gpr32','Cmklr1','Gpr37'))%>%
  add_row(name="LTB4 Receptors",
          genes=list('Ltb4r1','Ltb4r2')) %>%
  add_row(name="B-Cell Regulators",
          genes=list('Ccl28','Ccr10','Vcam1','Slc30a2','Ackr2')) #CCL28 binds to CCR10 on B-cells (see http://dx.doi.org/10.1007/s10911-010-9188-7 and its ref45). VCAM-1 blocking decreases IgA secretion (ref 52)



for (set in gene.sets$name) {
  genes <- unlist(filter(gene.sets, name==set)$genes)
    barplot <- long.merged.data %>%
    dplyr::filter(mgi_symbol %in% genes) %>%
        group_by(mgi_symbol,Genotype) %>%
    summarize(Mean = mean(TPM),
              Error =se(TPM)) %>%
    ggplot(aes(y=Mean,
               x=factor(mgi_symbol, levels=genes),
               fill=Genotype)) +
    geom_bar(stat='identity',position='dodge',width=0.75) +
    geom_errorbar(aes(ymin=Mean-Error,
                      ymax=Mean+Error),
                  width=0.5,
                  position=position_dodge(width=0.75)) +
    scale_fill_grey() +
    scale_color_grey() +
    theme_classic() +
    theme(text=element_text(size=18),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = c(0.2,0.8)) +
    labs(y="Normalized Counts (TPM)",
         x="",
         title=set) +
      scale_fill_grey()
    
  ggsave(barplot, 
         device='png',
         filename=c(paste('figures-noura/grouped-barplots/',paste(set,' Plot.png', sep=''), sep='')))
  
    ggsave(barplot, 
         device='pdf',
         filename=c(paste('figures-noura/grouped-barplots/',paste(set,' Plot.pdf', sep=''), sep='')))
}
```

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
## [1] ggplot2_3.3.3  biomaRt_2.46.3 readr_1.4.0    dplyr_1.0.5    tidyr_1.1.3   
## [6] knitr_1.31    
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.6           prettyunits_1.1.1    assertthat_0.2.1    
##  [4] digest_0.6.27        utf8_1.2.1           BiocFileCache_1.14.0
##  [7] R6_2.5.0             stats4_4.0.2         RSQLite_2.2.5       
## [10] evaluate_0.14        httr_1.4.2           pillar_1.5.1        
## [13] rlang_0.4.10         progress_1.2.2       curl_4.3            
## [16] rstudioapi_0.13      jquerylib_0.1.3      blob_1.2.1          
## [19] S4Vectors_0.28.1     rmarkdown_2.7        labeling_0.4.2      
## [22] stringr_1.4.0        bit_4.0.4            munsell_0.5.0       
## [25] compiler_4.0.2       xfun_0.22            pkgconfig_2.0.3     
## [28] askpass_1.1          BiocGenerics_0.36.0  htmltools_0.5.1.1   
## [31] openssl_1.4.3        tidyselect_1.1.0     tibble_3.1.0        
## [34] IRanges_2.24.1       XML_3.99-0.6         fansi_0.4.2         
## [37] crayon_1.4.1         dbplyr_2.1.0         withr_2.4.1         
## [40] rappdirs_0.3.3       grid_4.0.2           jsonlite_1.7.2      
## [43] gtable_0.3.0         lifecycle_1.0.0      DBI_1.1.1           
## [46] magrittr_2.0.1       scales_1.1.1         cli_2.3.1           
## [49] stringi_1.5.3        cachem_1.0.4         farver_2.1.0        
## [52] xml2_1.3.2           bslib_0.2.4          ellipsis_0.3.1      
## [55] generics_0.1.0       vctrs_0.3.7          tools_4.0.2         
## [58] bit64_4.0.5          Biobase_2.50.0       glue_1.4.2          
## [61] purrr_0.3.4          hms_1.0.0            parallel_4.0.2      
## [64] fastmap_1.1.0        yaml_2.2.1           AnnotationDbi_1.52.0
## [67] colorspace_2.0-0     memoise_2.0.0        sass_0.3.1
```
