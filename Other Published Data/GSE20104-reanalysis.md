---
title: "Re-Analysis of GSE20104 - ATF4 Overexpression Muscles"
author: "Dave Bridges"
date: "February 7, 2022"
output:
  html_document:
    keep_md: yes
  pdf_document:
    keep_tex: yes
---



This is generated via GEO2R.  Ignored fed and fasted mice


```r
# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE20104", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6096", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "01010101XXXXXXXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Control","ATF4 OE"))
levels(gs) <- groups
gset$group <- gs
gs <- relevel(gs, ref="Control")
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[2], groups[1], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_LIST","SPOT_ID","RANGE_GB","RANGE_STRAND","RANGE_START"))
#write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf) %>%
  separate(gene_assignment, into=c('Genbank','Gene.symbol','Gene Name','Other'),
           sep=" // ",
           extra="merge")
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
  ylab = "Number of genes", main = "P-adj value distribution")
```

![](figures/GSE20104-geo2r-1.png)<!-- -->

```r
# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())
```

![](figures/GSE20104-geo2r-2.png)<!-- -->

```r
# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
```

![](figures/GSE20104-geo2r-3.png)<!-- -->

```r
# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
```

[1] "ATF4.OE-Control"

```r
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
  highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
```

![](figures/GSE20104-geo2r-4.png)<!-- -->

```r
# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)
```

![](figures/GSE20104-geo2r-5.png)<!-- -->

```r
################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE20104", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
```

![](figures/GSE20104-geo2r-6.png)<!-- -->

```r
# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE20104", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")
```

![](figures/GSE20104-geo2r-7.png)<!-- -->

```r
# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 4, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=4", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
```

![](figures/GSE20104-geo2r-8.png)<!-- -->

```r
# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE20104")
```

![](figures/GSE20104-geo2r-9.png)<!-- -->

```r
#write to output file
output_file <- 'GSE20104 Analysis.csv'
write.fit(fit2, file=output_file, adjust='BH')

atf4.oe.results <- droplevels(topTable(fit2, n=Inf, adjust.method="BH")) %>%
  separate(gene_assignment, into=c('Genbank','Gene.symbol','Gene Name','Other'),
           sep=" // ",
           extra="merge") %>%
  select(Gene.symbol,logFC,P.Value,adj.P.Val) %>%
  distinct(Gene.symbol, .keep_all=T)
#annotate the probes
sig.atf4.oe.results <- atf4.oe.results %>% filter(P.Value<0.05)
sig.atf4.oe.genes <- sig.atf4.oe.results$Gene.symbol#used nominal p value, nothing sig after FDR

sig.atf4.oe.results.up <- atf4.oe.results %>% filter(P.Value<0.05&logFC>0)#used nominal p value, nothing sig after FDR
sig.atf4.oe.genes.up <- sig.atf4.oe.results.up$Gene.symbol

sig.atf4.oe.results.down <- atf4.oe.results %>% filter(P.Value<0.05&logFC<0)#used nominal p value, nothing sig after FDR
sig.atf4.oe.genes.down <- sig.atf4.oe.results.down$Gene.symbol

mTsc.datafile <- '../RNAseq/Muscle Tsc1 Knockout/data/processed/Binary DESeq Results.csv'
mtsc.data <- read.csv(mTsc.datafile)

sig.mtsc.data <- droplevels(subset(mtsc.data, padj<0.05))
sig.mtsc.data.up <- droplevels(subset(mtsc.data, padj<0.05&log2FoldChange>0))
sig.mtsc.data.down <- droplevels(subset(mtsc.data, padj<0.05&log2FoldChange<0))

sig.mtsc.genes <- sig.mtsc.data$external_gene_name
sig.mtsc.genes.up <- sig.mtsc.data.up$external_gene_name
sig.mtsc.genes.down <- sig.mtsc.data.down$external_gene_name

sig.overlap <- intersect(sig.atf4.oe.genes, sig.mtsc.genes)
sig.overlap.up <- intersect(sig.atf4.oe.genes.up, sig.mtsc.genes.up)
sig.overlap.down <- intersect(sig.atf4.oe.genes.down, sig.mtsc.genes.down)

fisher.table <-
  matrix(c(length(unique(sig.atf4.oe.genes)),
           length(unique(atf4.oe.results$Gene.symbol)), 
           length(sig.overlap),
           dim(sig.mtsc.data)[1]),
       nrow = 2,
       dimnames = list(KD.Sig = c("Sig", "Total"),
                       mTSC.Sig = c("ATF4", "mTSC")))

kable(fisher.table, caption="Contingency table for comparason of ATF4 and mTSC, all genes")
```



Table: Contingency table for comparason of ATF4 and mTSC, all genes

|      |  ATF4| mTSC|
|:-----|-----:|----:|
|Sig   |   598|  153|
|Total | 16481| 4403|

```r
fisher.table.up <-
  matrix(c(length(unique(sig.atf4.oe.genes.up)),
           length(unique(atf4.oe.results$Gene.symbol)), 
           length(sig.overlap.up),
           dim(sig.mtsc.data)[1]),
       nrow = 2,
       dimnames = list(KD.Sig = c("Sig", "Total"),
                       mTSC.Sig = c("ATF4", "mTSC")))

kable(fisher.table.up, caption="Contingency table for comparason of ATF4 and mTSC, up-regulated genes")
```



Table: Contingency table for comparason of ATF4 and mTSC, up-regulated genes

|      |  ATF4| mTSC|
|:-----|-----:|----:|
|Sig   |   352|   69|
|Total | 16481| 4403|

```r
fisher.table.down <-
  matrix(c(length(unique(sig.atf4.oe.genes.down)),
           length(unique(atf4.oe.results$Gene.symbol)), 
           length(sig.overlap.down),
           dim(sig.mtsc.data)[1]),
       nrow = 2,
       dimnames = list(KD.Sig = c("Sig", "Total"),
                       mTSC.Sig = c("ATF4", "mTSC")))

kable(fisher.table.down, caption="Contingency table for comparason of ATF4 and mTSC, down-regulated genes")
```



Table: Contingency table for comparason of ATF4 and mTSC, down-regulated genes

|      |  ATF4| mTSC|
|:-----|-----:|----:|
|Sig   |   246|   57|
|Total | 16481| 4403|

This identified 598 significantly differentially expressed genes in their analysis out of a total of 16481 genes assessed.  Of these differentially expressed genes. 153 genes overlapped with our 4403 significantly different genes from *Tsc1* knockout mice quadriceps.  This is an insignificant level of overlap (1.044 fold enrichment; p=0.682).

Focusing specifically on directionality, this is a significant level of overlap for up (1.363 times more likely; p=0.018) and dowregulated genes (1.153 times more likely; p=0.357).  

For upregulation this represents 352 ATF4 dependent genes out of 16481 genes assessed.  Of these differentially expressed genes. 69 genes overlapped with our 4403 significantly upregulated genes from *Tsc1* knockout mice quadriceps.  

For downregulation this represents 246 ATF4 dependent genes out of 16481 genes assessed.  Of these differentially expressed genes. 57 genes overlapped with our 4403 significantly downregulated genes from *Tsc1* knockout mice quadriceps.  


```r
combined.genes <-
  mtsc.data %>%
  select(external_gene_name,log2FoldChange,pvalue,padj) %>%
  left_join(atf4.oe.results, by=c('external_gene_name'='Gene.symbol')) %>%
  rename("log2FC_TSC"="log2FoldChange",
         "pval_TSC"='pvalue',
         "padj_TSC"='padj',
         "log2FC_KD"="logFC",
         'pval_KD'="P.Value",
         'padj_KD'='adj.P.Val')

with(combined.genes, lm(log2FC_TSC~log2FC_KD))  %>% tidy %>% kable(caption="Linear model for association between ATF4 fold change and TSC fold change")
```



Table: Linear model for association between ATF4 fold change and TSC fold change

|term        | estimate| std.error| statistic| p.value|
|:-----------|--------:|---------:|---------:|-------:|
|(Intercept) |    0.113|     0.007|      15.3|       0|
|log2FC_KD   |    0.688|     0.039|      17.6|       0|

```r
library(ggplot2)
combined.genes %>%
  ggplot(aes(y=log2FC_TSC,x=log2FC_KD)) +
  geom_point(size=0.1, alpha=0.1) +
  labs(y="TSC Fold Change",x="KD Fold Change") +
  geom_smooth(se=T)
```

![](figures/atf4-mtsc-gene-overlap-1.png)<!-- -->

```r
interesting.genes <- c('Sln','Gdf15') #
combined.genes %>%
  filter(external_gene_name %in% interesting.genes) %>%
  kable(caption="Selected genes")
```



Table: Selected genes

|external_gene_name | log2FC_TSC| pval_TSC| padj_TSC| log2FC_KD| pval_KD| padj_KD|
|:------------------|----------:|--------:|--------:|---------:|-------:|-------:|
|Gdf15              |       5.53|        0|        0|     0.195|   0.457|   1.000|
|Sln                |       4.25|        0|        0|     1.375|   0.001|   0.707|

```r
combined.genes %>%
  filter(external_gene_name %in% interesting.genes) %>%
  select(external_gene_name, log2FC_KD, log2FC_TSC) %>%
  group_by(external_gene_name) %>%
  pivot_longer(names_to = 'Experiment', values_to = 'Log2FC', cols=log2FC_KD:log2FC_TSC) %>%
  ggplot(aes(y=Log2FC,x=external_gene_name,
             fill=Experiment)) +
    geom_bar(stat='identity', position='dodge')
```

![](figures/atf4-mtsc-gene-overlap-2.png)<!-- -->


```r
require(venneuler)
v.diseases <- venneuler(c("ATF4"=length(sig.atf4.oe.genes), 
                 "Tsc1 Knockout Muscles"=length(sig.mtsc.genes),
                 "ATF4&Tsc1 Knockout Muscles"=length(intersect(sig.atf4.oe.genes, sig.mtsc.genes))))

v.diseases.up <- venneuler(c("ATF4"=length(sig.atf4.oe.genes.up), 
                 "Tsc1 Knockout Muscles"=length(sig.mtsc.genes.up),
                 "ATF4&Tsc1 Knockout Muscles"=length(intersect(sig.atf4.oe.genes.up, sig.mtsc.genes.up))))
v.diseases.down <- venneuler(c("ATF4"=length(sig.atf4.oe.genes.down), 
                 "Tsc1 Knockout Muscles"=length(sig.mtsc.genes.down),
                 "ATF4&Tsc1 Knockout Muscles"=length(intersect(sig.atf4.oe.genes.down, sig.mtsc.genes.down))))

plot(v.diseases, main="TSC-Dependent Transcriptional Changes")
```

![](figures/atf4-mtsc-venn-1.png)<!-- -->

```r
plot(v.diseases.up, main="TSC-Dependent and ATF4 Upregulation")
```

![](figures/atf4-mtsc-venn-2.png)<!-- -->

```r
plot(v.diseases, main="TSC-Dependent and ATF4 Downregulation")
```

![](figures/atf4-mtsc-venn-3.png)<!-- -->

```r
# library(Vennerable)
# v.list <- list(`Tsc2 Knockout MEFs` = sig.duvel.genes, `Tsc1 Knockout Muscles` = sig.mtsc.genes)
# v.data <- Venn(v.list)
# plot(v.data)
```

# Pathway Analysis


```r
library(clusterProfiler)
library(org.Mm.eg.db)
db.db.logfc <- 
  tT2 %>% 
  arrange(-logFC) %>%
  dplyr::select(logFC, Gene.symbol) 

db.db.list <- as.numeric(db.db.logfc$logFC)
names(db.db.list) <- db.db.logfc$Gene.symbol
gsea.bp <- gseGO(db.db.list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             OrgDb = org.Mm.eg.db)
library(enrichplot)
upsetplot(gsea.bp)
```

![](figures/atf4-pathway-analysis-1.png)<!-- -->

```r
dotplot(gsea.bp, orderBy="NES")
```

![](figures/atf4-pathway-analysis-2.png)<!-- -->

```r
gsea.bp <- pairwise_termsim(gsea.bp, method="JC")
emapplot(gsea.bp,
         cex_label_category=0.5,
         showCategory = 50,
         color="NES",
         min_edge=0.1)
```

![](figures/atf4-pathway-analysis-3.png)<!-- -->

```r
gsea.bp %>% as.data.frame %>% dplyr::select(NES,pvalue,p.adjust,Description, core_enrichment) %>% arrange(-NES) %>% filter(p.adjust<0.05) %>% kable(caption="Significant GO-BP pathways")
```



Table: Significant GO-BP pathways

|           |   NES| pvalue| p.adjust|Description                                                                          |core_enrichment                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
|:----------|-----:|------:|--------:|:------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|GO:0042501 |  2.01|      0|    0.043|serine phosphorylation of STAT protein                                               |Ifnz/Ifna14/Gadd45a/Ifna13/Ifna2/Ifnb1/Ifna11/Ifna9/Ifna12/Gfra2/Il24/Ifng/Ifna4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|GO:0003009 |  1.97|      0|    0.046|skeletal muscle contraction                                                          |Tnni1/Tnnc1/Chrna1/Chrnb1/Chrnd/Cav3/Chrng                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
|GO:0043043 | -1.48|      0|    0.025|peptide biosynthetic process                                                         |Ucn/Pa2g4/Ptbp2/Ptk2b/Elavl1/Eif4g3/Magoh/Sepsecs/Cdk4/Nars2/Chchd1/Rps6kb1/Lsm14b/Mtrf1/Xrn1/Nanos3/Eif3a/Tnrc6b/Barhl2/Eif2ak3/Aasdh/Mrpl35/Pld1/Mrps6/Rpl10a/Pstk/Mrps12/Ptcd3/Zcchc13/Lrpprc/Slc35a4/Eif3i/Eef1e1/Mrpl19/Mrpl37/Gclm/Eif1a/Cpeb4/Zc3h12d/Mrps16/Nck1/Yars2/Qk/Gclc/Eif4b/Ppp1r15b/Eif5a/Mrpl13/Eif5/Etf1/Rps24/Mrpl15/Tsc1/Hbs1l/Eif4a2/Qrsl1/Mrps18c/Dph3/Rpl8/Rock1/Rpl35a/Mrpl52/Rps27l/Mrps14/Rpusd4/Rpl4/Ythdf1/Mtrf1l/Mrpl21/Enc1/Eif4g1/Mtif2/Fmr1/Otud6b/Mrps9/Gfm2/Nck2/Syncrip/Yars/Rpl29/Rock2/Fastkd2/Mrpl32/Tarsl2/Rpsa/C1qbp/Cnot7/Ddx3x/Mmp7/Pum2/Eif4e3/Drg1/Cnot6/Cpeb1/Gfm1/Uhmk1/Paip2b/Dhx9/Mrps18b/Fastkd3/Gatc/Gspt1/Eif3e/Cnot1/Serp1/Larp1/Eif3d/Eif5b/Rps18/Mapk1/Usp16/Pink1/Mrps17/Cnot8/Mrpl1/Dapk3/Eif1/Eif4e/Eif2s2/Stat3/Mrpl16/Abce1/Mcts1/Hnrnpd/Prg3/Eif3f/Cdkal1/Mrpl51/Rnf139/Larp4/Lsm14a/Eif2ak2/Paip1/Rmnd1/Wars/Rpl31/Rpl21/Rpl37/Hbb-b1/S100a9                                                                                                                                                                                                                                                                                        |
|GO:0006412 | -1.54|      0|    0.004|translation                                                                          |Tnrc6a/Aco1/Ilf3/Ep300/Mrpl11/Itga2/Mettl5/Dhx36/Ireb2/Prkca/Per2/Mex3d/Ddx25/Wars2/Eif2a/Hfm1/Qk/Eif2s3x/Tpr/Ncl/Tars/Mrpl9/Fxr1/Ucn/Pa2g4/Ptbp2/Ptk2b/Elavl1/Eif4g3/Magoh/Sepsecs/Cdk4/Nars2/Chchd1/Rps6kb1/Lsm14b/Mtrf1/Xrn1/Nanos3/Eif3a/Tnrc6b/Barhl2/Eif2ak3/Mrpl35/Pld1/Mrps6/Rpl10a/Pstk/Mrps12/Ptcd3/Zcchc13/Lrpprc/Slc35a4/Eif3i/Eef1e1/Mrpl19/Mrpl37/Eif1a/Cpeb4/Zc3h12d/Mrps16/Nck1/Yars2/Qk/Eif4b/Ppp1r15b/Eif5a/Mrpl13/Eif5/Etf1/Rps24/Mrpl15/Tsc1/Hbs1l/Eif4a2/Qrsl1/Mrps18c/Dph3/Rpl8/Rock1/Rpl35a/Mrpl52/Rps27l/Mrps14/Rpusd4/Rpl4/Ythdf1/Mtrf1l/Mrpl21/Enc1/Eif4g1/Mtif2/Fmr1/Otud6b/Mrps9/Gfm2/Nck2/Syncrip/Yars/Rpl29/Rock2/Fastkd2/Mrpl32/Tarsl2/Rpsa/C1qbp/Cnot7/Ddx3x/Pum2/Eif4e3/Drg1/Cnot6/Cpeb1/Gfm1/Uhmk1/Paip2b/Dhx9/Mrps18b/Fastkd3/Gatc/Gspt1/Eif3e/Cnot1/Serp1/Larp1/Eif3d/Eif5b/Rps18/Mapk1/Usp16/Pink1/Mrps17/Cnot8/Mrpl1/Dapk3/Eif1/Eif4e/Eif2s2/Stat3/Mrpl16/Abce1/Mcts1/Hnrnpd/Prg3/Eif3f/Cdkal1/Mrpl51/Rnf139/Larp4/Lsm14a/Eif2ak2/Paip1/Rmnd1/Wars/Rpl31/Rpl21/Rpl37/Hbb-b1/S100a9                                                                                                                                                                           |
|GO:0006913 | -1.55|      0|    0.046|nucleocytoplasmic transport                                                          |Kpna4/Ranbp2/Efcab7/Xpo1/Ptgs2/Rab23/Nup160/Nup155/Pkia/Anp32a/Nsun2/Nup98/Nup50/Txnip/Ei24/Ywhae/Appl1/Eif5a/Ipo7/Apod/Prkcq/Tsc1/Ctdspl2/Ythdc1/Nup88/Ipo11/Ipo4/Rab18/Ier3/Heatr3/Kpna3/Hsp90aa1/Kpna1/Stradb/Fyttd1/Sfn/Tmco6/Dnajc27/Camk4/Thoc7/Rbm15b/Epm2a/Rpain/Uhmk1/Dhx9/Ppm1a/Nupl2/Inpp4b/Tnpo3/Ifi27/Cse1l/Mapk1/Snupn/Nup153/Pik3r1/Stat3/Tek/Abce1/Ubr5/Mlxip/Ltv1/Cd36/Angpt1/Rsrc1/Jak2/Nutf2/Ddx5/Hnrnpa2b1/Uaca/Gbp4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
|GO:0051169 | -1.55|      0|    0.046|nuclear transport                                                                    |Kpna4/Ranbp2/Efcab7/Xpo1/Ptgs2/Rab23/Nup160/Nup155/Pkia/Anp32a/Nsun2/Nup98/Nup50/Txnip/Ei24/Ywhae/Appl1/Eif5a/Ipo7/Apod/Prkcq/Tsc1/Ctdspl2/Ythdc1/Nup88/Ipo11/Ipo4/Rab18/Ier3/Heatr3/Kpna3/Hsp90aa1/Kpna1/Stradb/Fyttd1/Sfn/Tmco6/Dnajc27/Camk4/Thoc7/Rbm15b/Epm2a/Rpain/Uhmk1/Dhx9/Ppm1a/Nupl2/Inpp4b/Tnpo3/Ifi27/Cse1l/Mapk1/Snupn/Nup153/Pik3r1/Stat3/Tek/Abce1/Ubr5/Mlxip/Ltv1/Cd36/Angpt1/Rsrc1/Jak2/Nutf2/Ddx5/Hnrnpa2b1/Uaca/Gbp4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
|GO:1903362 | -1.57|      0|    0.049|regulation of cellular protein catabolic process                                     |Vps35/Rgn/Pias1/Egf/Nfe2l1/Styx/Kcne2/Rnf19a/Pten/Rchy1/Nkd2/Chfr/Sirt6/Mtm1/Ube2v2/Sumo2/Atxn3/Ptk2b/Rnf144b/Wac/Mapk8/Glmn/Usp14/Rnf14/Dnajb2/Bbs7/Fbxl5/Rdx/Aqp11/Psmf1/Tlk2/Efna1/Ube2k/Alad/Xpo1/Csnk1d/Sgta/Cav1/Ubqln4/Ptk2/Faf1/Dda1/Psen2/Derl2/Det1/Arih1/Gclc/Rnft1/Rad23b/Cdk5rap3/Ctsc/Csnk2a1/Sdcbp/Hsp90aa1/Fmr1/Arih2/Usp13/Bag2/Usp25/Agtpbp1/Taf9/Trim39/Epm2a/Wnt1/Mylip/Uchl5/Psme1/Psme3/Bag5/Psme2/N4bp1/Nell1/Usp7/Rnf139/Mapk9/Psme2/Gabarapl2/Sumo1/Taf9/Clu/Hspa1a/Socs4/Nub1/Hamp                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
|GO:0048872 | -1.58|      0|    0.029|homeostasis of number of cells                                                       |Ripk3/Bpgm/Zfx/Hif1a/Il7/Rb1/Vhl/Inhba/Ezh1/Csf1/Ppp2r3c/Il2ra/Ccnb2/Pla2g10/Exoc6/Sp3/Rhag/Ikbkg/Rps24/Ncor1/Itgam/Ccr2/Cdk5rap3/Ccr7/Lyn/Prdx5/Bcl6/Id2/Wdr48/Prdx2/Ptpn2/Pde4b/Tmod3/Acvr2a/Sfxn1/Zbtb7a/Zc3h8/Sos1/Vps54/Slc40a1/Dnaja3/Cfh/Sox6/Epas1/Sp1/B2m/Fas/Cyld/Napepld/Mtch2/Stat5b/Alas2/Hoxa5/Stat3/Acvr1b/Kras/Ezh2/Ccl2/Cd24a/Bmi1/Jak2/Ppp3cb/Prdx1/Ccr4/Siva1/Stat1/Hamp/Hbb-b1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|GO:0042176 | -1.59|      0|    0.007|regulation of protein catabolic process                                              |Ophn1/Foxo1/Snx9/Arntl/Vps35/Rgn/Pias1/Sec22b/Egf/Nfe2l1/Styx/Kcne2/Rnf19a/Pten/Hace1/Rchy1/Nkd2/Chfr/Nsf/Sirt6/Mtm1/Ube2v2/Sumo2/Huwe1/Atxn3/Ptk2b/Rnf144b/Mdm4/Wac/Oaz1/Trim32/Mapk8/Glmn/Oaz2/Serpinb1b/Stx5a/Usp14/Mycbp2/Rnf14/Dnajb2/Bbs7/Fbxl5/Rdx/Aqp11/Psmf1/Fam83d/Hecw2/Tlk2/Efna1/Sorl1/Ube2k/Alad/Xpo1/Csnk1d/Atg4b/Rhbdd1/Vhl/Sgta/Cav1/Adam9/Ubqln4/Itch/Ptk2/Faf1/Dda1/Psen2/Derl2/Det1/Arih1/Cyp51/Gclc/Rnft1/Rad23b/Timp1/Psmd2/Asb5/Timp3/Sox17/Nedd4l/Cdk5rap3/Vps28/Ctsc/Csnk2a1/Ier3/Sdcbp/Nrd1/Hsp90aa1/Smad4/Fmr1/Arih2/Usp13/Snca/Bag2/Usp25/Agtpbp1/Dedd/Taf9/Trim39/Epm2a/Wnt1/Mylip/Asb11/Uchl5/Psme1/Snf8/Psme3/Snx33/Tiparp/Bag5/Psme2/N4bp1/Nell1/Samd9l/Fbxl20/Usp7/Rnf139/Mapk9/Smurf2/Psme2/Gabarapl2/Sumo1/Taf9/Clu/Hspa1a/Socs4/Nub1/Egln1/Serpine2/Hamp                                                                                                                                                                                                                                                                                                                                                                                                       |
|GO:0006417 | -1.59|      0|    0.028|regulation of translation                                                            |Tnrc6a/Aco1/Ilf3/Ep300/Itga2/Mettl5/Dhx36/Ireb2/Prkca/Per2/Mex3d/Ddx25/Eif2a/Hfm1/Qk/Eif2s3x/Tpr/Ncl/Fxr1/Ucn/Pa2g4/Ptbp2/Ptk2b/Elavl1/Eif4g3/Magoh/Sepsecs/Cdk4/Rps6kb1/Lsm14b/Xrn1/Nanos3/Tnrc6b/Barhl2/Eif2ak3/Pld1/Pstk/Ptcd3/Zcchc13/Lrpprc/Slc35a4/Cpeb4/Zc3h12d/Nck1/Qk/Ppp1r15b/Eif5a/Mrpl13/Eif5/Etf1/Tsc1/Dph3/Rock1/Rps27l/Rpusd4/Ythdf1/Enc1/Eif4g1/Fmr1/Otud6b/Nck2/Syncrip/Rock2/Fastkd2/C1qbp/Cnot7/Ddx3x/Pum2/Eif4e3/Cnot6/Cpeb1/Uhmk1/Paip2b/Dhx9/Fastkd3/Gspt1/Eif3e/Cnot1/Serp1/Larp1/Eif3d/Eif5b/Mapk1/Usp16/Pink1/Cnot8/Dapk3/Eif4e/Stat3/Hnrnpd/Prg3/Rnf139/Larp4/Lsm14a/Eif2ak2/Paip1/Rmnd1/Hbb-b1/S100a9                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|GO:0009617 | -1.61|      0|    0.001|response to bacterium                                                                |Pde4b/Klrk1/Dusp10/Defb40/Pglyrp1/Fcgr4/Tbk1/Ifi44/Mmp7/Prdx3/Gsdmd/Slamf8/Ly86/Fkbp5/Reg4/Cxcl13/Spn/Paf1/Ociad2/B2m/Slc10a2/Nr1h4/Peli1/Stat5b/Mapkapk2/Seh1l/Ido1/Optn/Mapk1/Cflar/Tnfsf8/Rarres2/Pglyrp4/Rps6ka3/Hmgb1/Tnfaip8/Plscr1/Slfn2/Ccl12/Ccl5/Ccl2/Chd7/Car3/Cd24a/Cd36/Eif2ak2/Amy1/Trem3/Cfd/Ifi204/Abcd2/Mgst1/Lpl/Jak2/Serpina3f/Defb39/Adipoq/Prg2/Cxcl11/Plac8/Macrod2/Cyp2e1/H2-K1/Cxcl10/Ly6a/Gdap10/Ifi205/Psmb9/Gzma/Defb9/Gbp5/Defb5/Stat1/Gbp3/Cxcl9/Hamp/Cd274/Saa3/Gbp6/Gbp2/S100a9/Gbp4/Irgm1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
|GO:0043161 | -1.61|      0|    0.007|proteasome-mediated ubiquitin-dependent protein catabolic process                    |Fbxo45/Styx/Rnf19a/Hace1/Rchy1/Nkd2/Rnf126/Fbxo31/Chfr/Fbxw11/Psma1/Sel1l/Sirt6/Mtm1/Hspa5/Ube2v2/Sumo2/Bfar/Huwe1/Atxn3/Rnf4/Rnf144b/Klhl15/Wac/Siah2/Mapk8/Ctnnb1/Glmn/Cdc27/Usp14/Herc2/Rnf14/Dnajb2/Ube2a/Bbs7/Fbxl5/Hecw2/Tlk2/Psmd7/Cul1/Ube2k/Anapc4/Erlin2/Cdc16/Xpo1/Fbxw5/Csnk1d/Dnajc10/Sgta/Cav1/Psmc6/Crbn/Ubqln4/Itch/Dda1/Psen2/Derl2/Det1/Arih1/Rbck1/Gclc/Ubr3/Cd2ap/Rad23b/Sirt1/Psmd2/Ube2g1/Fbxl3/Kat5/Sec61b/Hectd3/Nedd4l/Hsp90b1/Fbxl12/4931417E11Rik/Sdcbp/Cul2/Ubxn2b/Siah1b/Fbxl15/Arih2/Fbxl4/Ppp2r5c/Ubr2/Ubr3/Taf9/Trim39/Topors/Tdpoz1/Ubr1/Epm2a/Fbxo6/Ube2i/Psmb11/Uchl5/Peli1/Cul3/Bag5/Rbx1/N4bp1/Spop/Fbxl20/Psmb10/Rmnd5a/Ube2b/Usp7/Cul5/Anapc5/Ube2d3/Mapk9/Smurf2/Sumo1/Taf9/Clu/Hspa1a/Trim13/Socs4/Nub1/Psma3/Psmb9/Psmb8                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|GO:0044403 | -1.62|      0|    0.007|symbiotic process                                                                    |Ifitm3/Gtf2b/Apobec3/Pi4ka/Srpk1/Etf1/Inpp5k/Rab43/Ddx6/Gpx2/Dicer1/Mphosph8/Ripk1/Top2a/Ubp1/Fmr1/Trim56/Oasl2/Rock2/Rab9/Chd1/Nbn/Bst2/Eps15/Pglyrp1/Thoc7/Cnot7/Ddx3x/Insr/Ube2i/Sp1/Dhx9/Exoc7/Napepld/Zfp639/Snf8/Ppid/Larp1/Nucks1/Eif3d/Trim31/Pomc/Ccnk/Npc1/Pglyrp4/Mcts1/Plscr1/Ccl5/Eif3f/Exoc2/Aqp1/Pkn2/Eif2ak2/Trem3/Trim21/Chmp1b/Trim13/Ddx5/Zbp1/Stat1/Gbp3/Gbp6/Gbp2/Irgm1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
|GO:0002831 | -1.62|      0|    0.021|regulation of response to biotic stimulus                                            |Prdx2/Vsig4/Ptpn2/Ifi203/Klrk1/Dusp10/Klrb1c/Tbk1/Dnaja3/C1qbp/Cnot7/Pum2/Slamf8/Ly86/Clec2d/Il4/Spn/Dhx9/Irf7/Usp15/Ffar2/Trafd1/Stat5b/Dtx3l/Xcl1/Optn/Parp9/Pomc/Hmgb1/Plscr1/Ccl5/Cd24a/Otop1/Irf1/Lsm14a/Tap1/Trem3/Ifi204/Fgl2/Zbp1/Nmi/Samhd1/Ifi205/Gbp5/Stat1/Parp14/Cd274/Gbp4/Irgm1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
|GO:0031331 | -1.63|      0|    0.007|positive regulation of cellular catabolic process                                    |Bnip3/Foxo1/Snx9/Fyco1/Tnrc6a/Vps35/Rgn/Pias1/Egf/Dhx36/Sesn1/Mex3d/Kcne2/Rnf19a/Pten/Bcl2l11/Trp53inp1/Rchy1/Pnpt1/Nkd2/Atm/Qk/Chfr/Supv3l1/Ulk2/Sirt6/Ube2v2/Sumo2/Huwe1/Atxn3/Ptk2b/Rnf144b/Wac/Mapk8/Fabp1/Nod1/Hsf1/Sh3glb1/Nanos3/Tnrc6b/Rnf14/Dnajb2/Bbs7/Fbxl5/Rdx/Grsf1/Prkd1/Mul1/Hif1a/Acsl5/Csnk1d/Sgta/Zc3h12d/Cav1/Tlr2/Adam9/Ptpn1/Itch/Ptk2/Faf1/Dda1/Psen2/Gapdhs/Det1/Arih1/Qk/Plekhf1/Gclc/Hk2/Rnft1/Sirt1/Smcr8/Gnai3/Tmem59/Ikbkg/Kat5/Tsc1/Slc4a4/Rock1/Ctsc/Mfn2/Pnpla2/Stk11/Prkaa2/Pik3c2a/Nrd1/Hsp90aa1/Ythdf1/Fmr1/Pafah1b2/Arih2/Usp13/Snca/Rock2/Bag2/Agtpbp1/Map3k7/Tbk1/Pip4k2b/Cnot7/Insr/Epm2a/Il4/Cnot1/Snx33/Dtx3l/Larp1/Optn/Pink1/Rnf152/Scoc/Cnot8/Adrb2/Hmgb1/Hnrnpd/Calcoco2/Rnf139/Mapk9/Trim21/Abcd2/Sumo1/Clu/Hspa1a/Trim13/Bnip3l/Socs4/Nub1/Hamp/Irgm1                                                                                                                                                                                                                                                                                                                                                                                                |
|GO:0009896 | -1.63|      0|    0.002|positive regulation of catabolic process                                             |Foxo1/Snx9/Fyco1/Tnrc6a/Vps35/Rgn/Pias1/Sec22b/Egf/Dhx36/Sesn1/Mex3d/Kcne2/Rnf19a/Pten/Hace1/Bcl2l11/Trp53inp1/Rchy1/Pnpt1/Nkd2/Atm/Qk/Chfr/Nsf/Prkce/Supv3l1/Ulk2/Sirt6/Ube2v2/Sumo2/Huwe1/Atxn3/Ptk2b/Rnf144b/Wac/Oaz1/Trim32/Mapk8/Fabp1/Oaz2/Nod1/Hsf1/Sh3glb1/Stx5a/Nanos3/Tnrc6b/Rnf14/Dnajb2/Bbs7/Fbxl5/Rdx/Grsf1/Prkd1/Hecw2/Mul1/Sorl1/Hif1a/Acsl5/Csnk1d/Atg4b/Rhbdd1/Sgta/Zc3h12d/Cav1/Tlr2/Adam9/Ptpn1/Itch/Ptk2/Faf1/Dda1/Psen2/Gapdhs/Det1/Arih1/Qk/Plekhf1/Gclc/Hk2/Rnft1/Sirt1/Smcr8/Asb5/Gnai3/Tmem59/Ikbkg/Kat5/Tsc1/Slc4a4/Rock1/Sox17/Nedd4l/Vps28/Ctsc/Mfn2/Csnk2a1/Ier3/Pnpla2/Stk11/Prkaa2/Pik3c2a/Nrd1/Hsp90aa1/Ythdf1/Fmr1/Pafah1b2/Arih2/Usp13/Snca/Rock2/Bag2/Agtpbp1/Map3k7/Tbk1/Pip4k2b/Cnot7/Insr/Epm2a/Il4/Mylip/Asb11/Snf8/Cnot1/Snx33/Tiparp/Dtx3l/Larp1/Optn/Pink1/Rnf152/Scoc/Cnot8/Adrb2/Hmgb1/Hnrnpd/Calcoco2/Rnf139/Mapk9/Smurf2/Trim21/Abcd2/Sumo1/Clu/Hspa1a/Trim13/Bnip3l/Socs4/Nub1/Hamp/Irgm1                                                                                                                                                                                                                                                           |
|GO:0016567 | -1.63|      0|    0.001|protein ubiquitination                                                               |Per2/Rnf38/Rnf19a/Pten/Hace1/Med17/Rchy1/Rnf220/Cbx8/Rnf114/Rnf126/Chfr/Fbxw11/Prkce/Ube2q1/Rnf135/Bcl10/Trim44/Hspa5/Ube2v2/Bfar/Huwe1/Rnf4/Ube2l3/Rnf144b/Klhl18/Klhl15/Wac/Siah2/Trim32/Prmt3/Mapk8/Cnot4/Ctnnb1/Hdac8/Cdc27/Herc4/Uba6/Otub2/Mycbp2/Herc2/Keap1/Rnf14/Fbxo28/Trim37/Dnajb2/Ube2a/Mib1/Fbxl5/Rnf6/Trim41/Birc3/Senp2/Med30/Rbbp6/Hecw2/Mul1/Ccnc/Cul1/Traf5/Ube3b/Ube2k/Anapc4/Ube2r2/Pja2/Cdc16/Fbxw5/Rnf20/Traf7/Vhl/Cav1/Arrdc4/Crbn/Asb1/Rnf8/Rnf141/Itch/Dda1/Psen2/Wdsub1/Gtpbp4/Bcor/Det1/Arih1/Rbck1/Gclc/Ubr3/Cdk8/Rnft1/Sirt1/Ube2l6/Vcpip1/Cdc73/Asb5/Ube2g1/Fbxl3/Cul7/Nedd4l/Cdk5rap3/Vps28/Klhl2/Trim23/Ube2e1/Cbll1/Cul2/Rnf7/Mtbp/Hsp90aa1/Enc1/Fbxl15/Wdr48/Dcun1d1/Trim56/Arih2/Anapc13/Bag2/Hltf/Ivns1abp/Nsmce2/Kcmf1/Ubr2/Ubr3/Dnaja3/Ddx3x/Trim39/Topors/Trim24/Ubr1/Epm2a/Pdzrn3/Paf1/Pef1/Atg3/Mylip/Asb11/Med21/Peli2/Anapc11/Trim54/Peli1/Rfwd3/Aktip/Cul3/Klhl7/Bag5/Dtx3l/Suz12/Med6/Dcun1d3/Pink1/Rnf152/Trim31/Rbx1/Dnaja1/Uhrf1/Rnf213/N4bp1/Spop/Wdr70/Adrb2/Ubr5/Asb4/Ube2n/Ube2b/Cul5/Anapc5/Fbxo32/Ube2d3/Ubc/Rnf139/Mapk9/Smurf2/Trim21/Angpt1/Bmi1/Magel2/Asb12/Trim13/Nub1/Rnf113a2/Pcgf5/Ubd/Nmi/Dcun1d2/Rnf146/Ube2q2/Gbp4              |
|GO:0006511 | -1.66|      0|    0.000|ubiquitin-dependent protein catabolic process                                        |Rnf19a/Pten/Hace1/Rchy1/Rnf114/Nkd2/Rnf126/Fbxo31/Chfr/Fbxw11/Psma1/Sel1l/Sirt6/Mtm1/Hspa5/Ube2v2/Psmc2/Sumo2/Bfar/Huwe1/Ate1/Vps4a/Atxn3/Rnf4/Ube2l3/Ptk2b/Rnf144b/Klhl15/Wac/Siah2/Trim32/Mapk8/Cnot4/Ctnnb1/Glmn/Cdc27/Uba6/Usp33/Usp14/Usp37/Herc2/Keap1/Rnf14/Usp24/Dnajb2/Ube2a/Mib1/Bbs7/Usp47/Usp9x/Fbxl5/Rnf6/Psmd11/Psmf1/Rbbp6/Hecw2/Tlk2/Psmd7/Cul1/Ube3b/Ube2k/Anapc4/Erlin2/Ube2r2/Cdc16/Xpo1/Ubap1/Fbxw5/Rnf20/Csnk1d/Dnajc10/Sgta/Cav1/Psmc6/Crbn/Psmd13/Ubqln4/Rnf8/Itch/Ptk2/Dda1/Usp1/Psen2/Derl2/Tollip/Det1/Arih1/Rbck1/Gclc/Ubr3/Usp45/Cd2ap/Rad23b/Sirt1/Ube2l6/Psmd2/Ube2g1/Fbxl3/Cul7/Kat5/Sec61b/Hectd3/Nedd4l/Vps28/Hsp90b1/Fbxl12/Csnk2a1/4931417E11Rik/Sdcbp/Cul2/Ubxn2b/Rnf7/Siah1b/Ubl7/Fbxl15/Arih2/Usp13/Fbxl4/Usp25/Ppp2r5c/Agtpbp1/Ubr2/Ubr3/Taf9/Trim39/Topors/Tdpoz1/Ubr1/Epm2a/Fbxo6/Ube2i/Wnt1/Mylip/Psmb11/Usp15/Uchl5/Cyld/Anapc11/Zranb1/Peli1/Snf8/Cul3/Ntan1/Bag5/Cops3/Dtx3l/Fbxo8/Usp16/Fem1a/Rbx1/Uhrf1/Rnf213/N4bp1/Usp10/Spop/Vps37a/Fbxl20/Psmb10/Ube2n/Rmnd5a/Ube2b/Usp7/Cul5/Anapc5/Ube2d3/Rnf139/Mapk9/Smurf2/Sumo1/Kctd6/Taf9/Clu/Hspa1a/Trim13/Socs4/Nub1/Ubd/Psma3/Psmb9/Rnf146/Psmb8                                                      |
|GO:0010498 | -1.67|      0|    0.001|proteasomal protein catabolic process                                                |Kcne2/Rnf19a/Hace1/Rchy1/Nkd2/Rnf126/Fbxo31/Chfr/Fbxw11/Psma1/Sel1l/Sirt6/Mtm1/Hspa5/Ube2v2/Sumo2/Bfar/Huwe1/Ate1/Atxn3/Rnf4/Rnf144b/Nr1d1/Klhl15/Wac/Siah2/Mapk8/Ctnnb1/Glmn/Cdc27/Usp14/Herc2/Rnf14/Dnajb2/Ube2a/Bbs7/Fbxl5/Aqp11/Psmf1/Hecw2/Tlk2/Psmd7/Cul1/Ube2k/Anapc4/Alad/Erlin2/Cdc16/Xpo1/Fbxw5/Csnk1d/Rhbdd1/Vhl/Dnajc10/Sgta/Cav1/Psmc6/Ddi2/Crbn/Ubqln4/Itch/Dda1/Psen2/Derl2/Det1/Arih1/Rbck1/Gclc/Ubr3/Cd2ap/Rnft1/Rad23b/Sirt1/Psmd2/Ube2g1/Fbxl3/Kat5/Sec61b/Hectd3/Nedd4l/Hsp90b1/Fbxl12/4931417E11Rik/Sdcbp/Cul2/Ubxn2b/Siah1b/Enc1/Fmr1/Fbxl15/Arih2/Usp13/Fbxl4/Bag2/Usp25/Ppp2r5c/Ubr2/Ubr3/Taf9/Trim39/Topors/Tdpoz1/Ubr1/Epm2a/Fbxo6/Ube2i/Psmb11/Uchl5/Peli1/Psme1/Psme3/Cul3/Bag5/Psme2/Rbx1/N4bp1/Spop/Fbxl20/Psmb10/Rmnd5a/Ube2b/Usp7/Cul5/Anapc5/Psme4/Ube2d3/Rnf139/Mapk9/Smurf2/Psme2/Gabarapl2/Sumo1/Taf9/Clu/Hspa1a/Trim13/Socs4/Nub1/Psma3/Psmb9/Psmb8                                                                                                                                                                                                                                                                                                           |
|GO:0019941 | -1.67|      0|    0.000|modification-dependent protein catabolic process                                     |Rnf19a/Pten/Hace1/Rchy1/Rnf114/Nkd2/Rnf126/Fbxo31/Chfr/Fbxw11/Psma1/Sel1l/Sirt6/Mtm1/Hspa5/Ube2v2/Psmc2/Sumo2/Bfar/Huwe1/Ate1/Vps4a/Atxn3/Rnf4/Ube2l3/Ptk2b/Rnf144b/Klhl15/Wac/Siah2/Trim32/Mapk8/Cnot4/Ctnnb1/Glmn/Lonp1/Cdc27/Uba6/Usp33/Usp14/Usp37/Herc2/Keap1/Rnf14/Usp24/Dnajb2/Ube2a/Mib1/Bbs7/Usp47/Usp9x/Fbxl5/Rnf6/Psmd11/Psmf1/Rbbp6/Hecw2/Tlk2/Psmd7/Cul1/Ube3b/Ube2k/Anapc4/Erlin2/Ube2r2/Cdc16/Xpo1/Ubap1/Fbxw5/Rnf20/Csnk1d/Zmpste24/Dnajc10/Sgta/Cav1/Psmc6/Crbn/Psmd13/Ubqln4/Rnf8/Itch/Ptk2/Dda1/Usp1/Psen2/Derl2/Tollip/Det1/Arih1/Rbck1/Gclc/Ubr3/Usp45/Cd2ap/Rad23b/Sirt1/Ube2l6/Psmd2/Ube2g1/Fbxl3/Cul7/Kat5/Sec61b/Hectd3/Nedd4l/Vps28/Hsp90b1/Fbxl12/Csnk2a1/4931417E11Rik/Sdcbp/Cul2/Ubxn2b/Rnf7/Siah1b/Ubl7/Fbxl15/Arih2/Usp13/Fbxl4/Usp25/Ppp2r5c/Agtpbp1/Ubr2/Ubr3/Taf9/Trim39/Topors/Tdpoz1/Ubr1/Epm2a/Fbxo6/Ube2i/Wnt1/Mylip/Psmb11/Usp15/Uchl5/Cyld/Anapc11/Zranb1/Peli1/Snf8/Cul3/Ntan1/Bag5/Cops3/Dtx3l/Fbxo8/Usp16/Fem1a/Rbx1/Uhrf1/Rnf213/N4bp1/Usp10/Spop/Vps37a/Fbxl20/Psmb10/Ube2n/Rmnd5a/Ube2b/Usp7/Cul5/Anapc5/Ube2d3/Ubc/Rnf139/Mapk9/Smurf2/Sumo1/Kctd6/Taf9/Clu/Hspa1a/Trim13/Socs4/Nub1/Ubd/Psma3/Psmb9/Rnf146/Psmb8                                   |
|GO:0016071 | -1.69|      0|    0.000|mRNA metabolic process                                                               |Clns1a/Snw1/Dhx36/Prpf8/Ddx23/Cdc40/Ppp4r2/Zcchc8/Polr2i/Mex3d/Ppie/Mbnl1/Pnpt1/Atm/Qk/Ddx47/Tut1/Ncl/Pnn/Supv3l1/Rbm27/Ddx1/Fxr1/Brf1/Cwc22/Ptbp2/Zhx2/Elavl1/Magoh/Rbm5/Dhx15/Srpk2/Tbrg4/Setx/Isy1/Tcerg1/Cdc5l/Hsf1/Prpf6/Wbp4/Xrn1/Nanos3/Tnrc6b/Cstf2t/Zmat2/Ap3b1/Rngtt/Dnajb11/Ubl5/Mfap1b/Grsf1/Exosc9/Prpf40a/Rbbp6/Sltm/Cpsf2/Phf5a/Frg1/Zcrb1/Ppwd1/Rnf20/Rbm6/Zc3h12d/Dyrk1a/Rnmt/Larp7/Rbm25/Nsun2/Nup98/Snrnp35/Sf3b1/Ddx41/Kin/Qk/Prdx6/Luc7l2/Fus/Raver2/Cdc73/Smu1/Cwf19l2/Srpk1/Etf1/Prpf38a/Sf3a3/U2af1l4/Pus3/Ythdc1/Rock1/Zc3h14/Crnkl1/Ccnt1/Dicer1/Hipk3/Strap/Cbll1/Ddx46/Hnrnpk/Pcf11/Myd88/Cd2bp2/Zrsr1/Znrd1/Ythdf1/Rprd1b/Fmr1/Hnrnph1/Csde1/Syncrip/Xrn2/Exosc7/Rock2/Prpf18/Zbtb7a/Clp1/Rbm41/Bag4/Rprd2/Thoc7/Lsm10/C1qbp/Cnot7/Pum2/Rbm15b/Mrto4/Cpsf4l/Cnot6/Cwc22/Cpeb1/Dbr1/Sp1/Magohb/Sf3a2/Zfp830/Dhx9/Paf1/Hnrnpm/Fastkd3/Exosc5/Khdrbs3/Exosc3/Smn1/Gspt1/Ptcd2/Prmt5/Eif3e/Cnot1/Pcbp4/Mapkapk2/Larp1/Snrpb2/Ddx20/Tdrd3/Cnot8/Htatsf1/Stat3/Snrpn/Hnrnpd/Angel2/Rbmx2/Exosc8/Rnps1/Snrpa1/Lsm8/Rsrc1/Rbmx/Ddx5/Rnf113a2/Smndc1/Hnrnpa2b1/Snrpd1/Rbm7/Tsen15/Slbp                                                                                         |
|GO:0006397 | -1.69|      0|    0.002|mRNA processing                                                                      |Clns1a/Snw1/Dhx36/Prpf8/Ddx23/Cdc40/Ppp4r2/Zcchc8/Ppie/Mbnl1/Pnpt1/Qk/Ddx47/Tut1/Ncl/Pnn/Rbm27/Ddx1/Fxr1/Cwc22/Ptbp2/Magoh/Rbm5/Dhx15/Srpk2/Tbrg4/Setx/Isy1/Tcerg1/Cdc5l/Hsf1/Prpf6/Wbp4/Cstf2t/Zmat2/Rngtt/Ubl5/Mfap1b/Grsf1/Prpf40a/Rbbp6/Sltm/Cpsf2/Phf5a/Frg1/Zcrb1/Ppwd1/Rnf20/Rbm6/Dyrk1a/Rnmt/Larp7/Rbm25/Nup98/Snrnp35/Sf3b1/Ddx41/Kin/Qk/Prdx6/Luc7l2/Raver2/Cdc73/Smu1/Cwf19l2/Srpk1/Prpf38a/Sf3a3/U2af1l4/Ythdc1/Zc3h14/Crnkl1/Ccnt1/Strap/Ddx46/Hnrnpk/Pcf11/Cd2bp2/Zrsr1/Rprd1b/Fmr1/Hnrnph1/Syncrip/Xrn2/Prpf18/Zbtb7a/Clp1/Rbm41/Rprd2/Thoc7/Lsm10/C1qbp/Rbm15b/Cpsf4l/Cwc22/Cpeb1/Dbr1/Magohb/Sf3a2/Zfp830/Dhx9/Paf1/Hnrnpm/Khdrbs3/Smn1/Ptcd2/Prmt5/Pcbp4/Snrpb2/Ddx20/Tdrd3/Htatsf1/Snrpn/Rbmx2/Rnps1/Snrpa1/Lsm8/Rsrc1/Rbmx/Ddx5/Rnf113a2/Smndc1/Hnrnpa2b1/Snrpd1/Rbm7/Tsen15/Slbp                                                                                                                                                                                                                                                                                                                                                                                              |
|GO:0043632 | -1.70|      0|    0.000|modification-dependent macromolecule catabolic process                               |Rnf19a/Pten/Hace1/Rchy1/Pnpt1/Rnf114/Nkd2/Rnf126/Fbxo31/Chfr/Fbxw11/Psma1/Sel1l/Sirt6/Mtm1/Hspa5/Ube2v2/Psmc2/Sumo2/Bfar/Huwe1/Ate1/Vps4a/Atxn3/Rnf4/Ube2l3/Ptk2b/Rnf144b/Klhl15/Wac/Siah2/Trim32/Mapk8/Cnot4/Ctnnb1/Glmn/Lonp1/Cdc27/Uba6/Usp33/Usp14/Usp37/Herc2/Keap1/Rnf14/Usp24/Dnajb2/Ube2a/Mib1/Bbs7/Usp47/Usp9x/Fbxl5/Rnf6/Psmd11/Exosc9/Psmf1/Rbbp6/Hecw2/Tlk2/Psmd7/Cul1/Ube3b/Ube2k/Anapc4/Erlin2/Ube2r2/Cdc16/Xpo1/Ubap1/Fbxw5/Rnf20/Csnk1d/Zmpste24/Dnajc10/Sgta/Cav1/Psmc6/Crbn/Psmd13/Ubqln4/Rnf8/Itch/Ptk2/Dda1/Usp1/Psen2/Derl2/Tollip/Det1/Arih1/Rbck1/Gclc/Ubr3/Usp45/Cd2ap/Rad23b/Sirt1/Ube2l6/Psmd2/Ube2g1/Fbxl3/Cul7/Kat5/Sec61b/Hectd3/Nedd4l/Vps28/Hsp90b1/Fbxl12/Csnk2a1/4931417E11Rik/Sdcbp/Cul2/Ubxn2b/Rnf7/Siah1b/Ubl7/Fbxl15/Exosc7/Arih2/Usp13/Fbxl4/Usp25/Ppp2r5c/Agtpbp1/Ubr2/Ubr3/Taf9/Trim39/Topors/Tdpoz1/Ubr1/Epm2a/Fbxo6/Ube2i/Wnt1/Mylip/Psmb11/Usp15/Uchl5/Cyld/Anapc11/Zranb1/Exosc3/Peli1/Snf8/Cul3/Ntan1/Bag5/Cops3/Dtx3l/Fbxo8/Usp16/Fem1a/Rbx1/Uhrf1/Rnf213/N4bp1/Usp10/Spop/Vps37a/Fbxl20/Psmb10/Ube2n/Rmnd5a/Ube2b/Usp7/Cul5/Anapc5/Ube2d3/Exosc8/Ubc/Rnf139/Mapk9/Smurf2/Sumo1/Kctd6/Taf9/Clu/Hspa1a/Trim13/Socs4/Nub1/Ubd/Psma3/Psmb9/Rnf146/Psmb8 |
|GO:0008380 | -1.74|      0|    0.002|RNA splicing                                                                         |Clns1a/Snw1/Prpf8/Ddx23/Cdc40/Ppp4r2/Zcchc8/Ppie/Mbnl1/Qk/Ddx47/Ncl/Pnn/Ddx1/Fxr1/Cwc22/Ptbp2/Magoh/Rbm5/Dhx15/Srpk2/Setx/Isy1/Tcerg1/Cdc5l/Prpf6/Wbp4/Zmat2/Ubl5/Mfap1b/Grsf1/Prpf40a/Phf5a/Frg1/Zcrb1/Ppwd1/Rbm6/Dyrk1a/Larp7/Rbm25/Nup98/Snrnp35/Sf3b1/Ddx41/Qk/Prdx6/Clk4/Luc7l2/Clk1/Fus/Raver2/Smu1/Cwf19l2/Srpk1/Prpf38a/Sf3a3/U2af1l4/Ythdc1/Crnkl1/Strap/Ddx46/Hnrnpk/Cd2bp2/Zrsr1/Fmr1/Hnrnph1/Syncrip/Prpf18/Zbtb7a/Clp1/Ivns1abp/Rbm41/Thoc7/Lsm10/C1qbp/Rbm15b/Cwc22/Dbr1/Magohb/Sf3a2/Zfp830/Dhx9/Hnrnpm/Khdrbs3/Trpt1/Smn1/Prmt5/Rnf113a1/Pcbp4/Snrpb2/Ddx20/Htatsf1/Pik3r1/Snrpn/Rbmx2/Rnps1/Snrpa1/Lsm8/Rsrc1/Rbmx/Ddx5/Rnf113a2/Smndc1/Hnrnpa2b1/Snrpd1/Rbm7/Tsen15                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
|GO:0000375 | -1.77|      0|    0.007|RNA splicing, via transesterification reactions                                      |Clns1a/Snw1/Prpf8/Ddx23/Cdc40/Ppie/Mbnl1/Qk/Ncl/Fxr1/Cwc22/Ptbp2/Magoh/Rbm5/Srpk2/Setx/Isy1/Cdc5l/Prpf6/Wbp4/Zmat2/Ubl5/Mfap1b/Prpf40a/Phf5a/Zcrb1/Rbm6/Dyrk1a/Larp7/Rbm25/Nup98/Snrnp35/Sf3b1/Ddx41/Qk/Prdx6/Luc7l2/Raver2/Smu1/Cwf19l2/Srpk1/Prpf38a/Sf3a3/U2af1l4/Ythdc1/Crnkl1/Strap/Ddx46/Hnrnpk/Zrsr1/Fmr1/Prpf18/Zbtb7a/Rbm41/C1qbp/Rbm15b/Cwc22/Dbr1/Magohb/Sf3a2/Dhx9/Hnrnpm/Khdrbs3/Smn1/Prmt5/Pcbp4/Snrpb2/Ddx20/Htatsf1/Snrpn/Rbmx2/Rnps1/Snrpa1/Lsm8/Rsrc1/Rbmx/Ddx5/Rnf113a2/Hnrnpa2b1/Snrpd1/Rbm7                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|GO:0000377 | -1.77|      0|    0.007|RNA splicing, via transesterification reactions with bulged adenosine as nucleophile |Clns1a/Snw1/Prpf8/Ddx23/Cdc40/Ppie/Mbnl1/Qk/Ncl/Fxr1/Cwc22/Ptbp2/Magoh/Rbm5/Srpk2/Setx/Isy1/Cdc5l/Prpf6/Wbp4/Zmat2/Ubl5/Mfap1b/Prpf40a/Phf5a/Zcrb1/Rbm6/Dyrk1a/Larp7/Rbm25/Nup98/Snrnp35/Sf3b1/Ddx41/Qk/Prdx6/Luc7l2/Raver2/Smu1/Cwf19l2/Srpk1/Prpf38a/Sf3a3/U2af1l4/Ythdc1/Crnkl1/Strap/Ddx46/Hnrnpk/Zrsr1/Fmr1/Prpf18/Zbtb7a/Rbm41/C1qbp/Rbm15b/Cwc22/Dbr1/Magohb/Sf3a2/Dhx9/Hnrnpm/Khdrbs3/Smn1/Prmt5/Pcbp4/Snrpb2/Ddx20/Htatsf1/Snrpn/Rbmx2/Rnps1/Snrpa1/Lsm8/Rsrc1/Rbmx/Ddx5/Rnf113a2/Hnrnpa2b1/Snrpd1/Rbm7                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|GO:0000398 | -1.77|      0|    0.007|mRNA splicing, via spliceosome                                                       |Clns1a/Snw1/Prpf8/Ddx23/Cdc40/Ppie/Mbnl1/Qk/Ncl/Fxr1/Cwc22/Ptbp2/Magoh/Rbm5/Srpk2/Setx/Isy1/Cdc5l/Prpf6/Wbp4/Zmat2/Ubl5/Mfap1b/Prpf40a/Phf5a/Zcrb1/Rbm6/Dyrk1a/Larp7/Rbm25/Nup98/Snrnp35/Sf3b1/Ddx41/Qk/Prdx6/Luc7l2/Raver2/Smu1/Cwf19l2/Srpk1/Prpf38a/Sf3a3/U2af1l4/Ythdc1/Crnkl1/Strap/Ddx46/Hnrnpk/Zrsr1/Fmr1/Prpf18/Zbtb7a/Rbm41/C1qbp/Rbm15b/Cwc22/Dbr1/Magohb/Sf3a2/Dhx9/Hnrnpm/Khdrbs3/Smn1/Prmt5/Pcbp4/Snrpb2/Ddx20/Htatsf1/Snrpn/Rbmx2/Rnps1/Snrpa1/Lsm8/Rsrc1/Rbmx/Ddx5/Rnf113a2/Hnrnpa2b1/Snrpd1/Rbm7                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|GO:0060759 | -1.80|      0|    0.030|regulation of response to cytokine stimulus                                          |Mul1/Hpx/Ube2k/Padi2/Il7/Vrk2/Cav1/Tlr2/Ptprc/Csf1/Il1f5/Cd300lf/Casp1/Ripk1/Ptpn2/Hipk1/Tbk1/Dnaja3/Cnot7/Taf9/Dhx9/Irf7/Cyld/Nr1h4/Pafah1b1/Parp9/Ifih1/Cd24a/Otop1/Lsm14a/Casp4/Angpt1/Taf9/Adipoq/Rnf113a2/Zbp1/Samhd1/Parp14/Irgm1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
|GO:0019395 | -1.80|      0|    0.046|fatty acid oxidation                                                                 |Acadl/Cpt2/Fabp1/Adipor1/Pex13/Hacl1/Pex7/Acsl5/Ppargc1a/Acad11/Adipor2/Abcd3/Auh/Acadvl/Etfdh/Etfb/Pdk4/Crot/Fabp3/Acox1/Lonp2/Echdc2/Acacb/Crat/Phyh/Hadha/Hsd17b4/Acat1/Echs1/Slc25a17/Acat3/Echdc1/Acadm/Cd36/Abcd2/C1qtnf9/Adipoq                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|GO:2001252 | -1.89|      0|    0.001|positive regulation of chromosome organization                                       |Setdb1/Wbp2/Jdp2/Terc/Men1/Xrcc5/Daxx/Akap8l/Lrrk2/Tnks/Ruvbl2/Jarid2/Rif1/Ep300/Nos1/Snw1/Dhx36/Nipbl/Dnmt1/Cct3/Atrx/Atm/Sdr16c5/Tpr/Sirt6/Cct2/Lig4/Mecp2/Gcg/Cct4/Dmrtc2/Ctnnb1/Pphln1/Mtf2/Cdc27/Smc5/Mre11a/Slk/Nek2/Gnl3/Prkd1/Pot1a/Kat2a/Anapc4/Cdc16/Rnf20/Rb1/Wdr5/Rad21/Hmbox1/Sirt1/Fen1/Rps6ka5/Prkcq/Terf1/Brd7/Ncor1/Gtf2h2/Mphosph8/Suv39h1/Map2k7/Bcl6/Smad4/Fmr1/Terf2/Nbn/Nsmce2/Paf1/Anapc11/Cul3/Mapk1/Pink1/Ing2/Nek7/Mier1/Hnrnpd/Lpin1/Ube2n/Anapc5/Kat2b/Hnrnpa2b1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
|GO:0033866 | -2.03|      0|    0.043|nucleoside bisphosphate biosynthetic process                                         |Ppcdc/Acsl5/Pdhx/Pank1/Pank4/Pank2/Acot7/Acsl1/Pdk4/Pdk2/Acacb/Dld/Snca/Acat1/Pdhb/Ppcs/Pank3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
|GO:0034030 | -2.03|      0|    0.043|ribonucleoside bisphosphate biosynthetic process                                     |Ppcdc/Acsl5/Pdhx/Pank1/Pank4/Pank2/Acot7/Acsl1/Pdk4/Pdk2/Acacb/Dld/Snca/Acat1/Pdhb/Ppcs/Pank3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
|GO:0034033 | -2.03|      0|    0.043|purine nucleoside bisphosphate biosynthetic process                                  |Ppcdc/Acsl5/Pdhx/Pank1/Pank4/Pank2/Acot7/Acsl1/Pdk4/Pdk2/Acacb/Dld/Snca/Acat1/Pdhb/Ppcs/Pank3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
|GO:0043124 | -2.06|      0|    0.018|negative regulation of I-kappaB kinase/NF-kappaB signaling                           |Sirt1/Ripk1/Tnip3/Dnaja3/Trim39/Zmynd11/Ppm1a/Nr1h4/Optn/Rora/Usp10/Casp8/Adipoq/Stat1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|GO:0044406 | -2.09|      0|    0.026|adhesion of symbiont to host                                                         |Gbp3/Gbp6/Gbp2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
|GO:0060333 | -2.16|      0|    0.014|interferon-gamma-mediated signaling pathway                                          |Ptpn2/Dnaja3/Parp9/Otop1/Irf1/Jak2/Nmi/Stat1/Parp14/Irgm1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
|GO:0034341 | -2.21|      0|    0.000|response to interferon-gamma                                                         |Stxbp4/Syncrip/Snca/Ptpn2/Bst2/Ccl11/Dnaja3/Vamp3/Stx8/Ccl8/Vamp4/Xcl1/Parp9/Actr2/Ccl7/Dapk3/Ccl12/Ccl5/Ccl2/Otop1/Irf1/Trim21/Jak2/Ubd/Nmi/Gbp5/Stat1/Gbp3/Parp14/Gbp6/Gbp2/Gbp4/Irgm1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
|GO:0071346 | -2.28|      0|    0.000|cellular response to interferon-gamma                                                |Stxbp4/Syncrip/Ptpn2/Ccl11/Dnaja3/Vamp3/Stx8/Ccl8/Vamp4/Xcl1/Parp9/Actr2/Ccl7/Dapk3/Ccl12/Ccl5/Ccl2/Otop1/Irf1/Jak2/Nmi/Gbp5/Stat1/Gbp3/Parp14/Gbp6/Gbp2/Gbp4/Irgm1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
|GO:0035456 | -2.33|      0|    0.001|response to interferon-beta                                                          |Plscr1/Irf1/Ifi204/Ifi205/Stat1/Gbp3/Gbp6/Gbp2/Irgm1/Igtp                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
|GO:0035458 | -2.36|      0|    0.000|cellular response to interferon-beta                                                 |Irf1/Ifi204/Ifi205/Stat1/Gbp3/Gbp6/Gbp2/Irgm1/Igtp                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |

## Empirical Muscle ATF4 Gene Set


```r
list(ATF4q0.25 = tT2 %>% filter(adj.P.Val <0.25) %>% pull(Gene.symbol),
     ATF4p0.001 = tT2 %>% filter(P.Value <0.001) %>% pull(Gene.symbol),
     ATF4p0.01 = tT2 %>% filter(P.Value <0.01) %>% pull(Gene.symbol),
     ATF4p0.05 = tT2 %>% filter(P.Value <0.05) %>% pull(Gene.symbol)) -> atf4.gene.lists
library(fgsea)

mtsc.ranks <- pull(mtsc.data, log2FoldChange)
names(mtsc.ranks) <- mtsc.data$external_gene_name

mtsc.atf4.pathways <- fgsea(pathways = atf4.gene.lists, 
                  stats    = sort(mtsc.ranks,decreasing=T),
                  minSize  = 1,
                  maxSize  = 500)

mtsc.atf4.pathways %>% 
  arrange(-NES) %>%
  kable(caption="GSEA results for mTSC differential expression results compared to empirically determined ATF4 dependent changes in muscle")
```



Table: GSEA results for mTSC differential expression results compared to empirically determined ATF4 dependent changes in muscle

|pathway    |  pval|  padj| log2err|    ES|  NES| size|leadingEdge                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
|:----------|-----:|-----:|-------:|-----:|----:|----:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|ATF4p0.001 | 0.000| 0.000|   0.593| 0.855| 2.06|   15|Uchl1 , Sln   , Krt18 , Ankrd1, Myl4  , Tnnt2 , Slc7a5, Tnni1 , Asns  , Igf2  , Chrna1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
|ATF4p0.01  | 0.000| 0.000|   0.675| 0.612| 2.06|   86|Uchl1  , Sln    , Krt18  , Ankrd1 , Myl4   , Tnnt2  , Zfp697 , Casq2  , Slc7a5 , Ncam1  , Tnni1  , Tnnc1  , Asns   , Chrnd  , Ptpn5  , Cyb5r3 , Maged2 , Nes    , Epor   , Slc7a3 , Olfr550, Hspb7  , Slpi   , Kcnn3  , Tmem37 , Slc7a11, Igf2   , Csrp3  , Gadd45a, Chrna1 , Dhrs9                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
|ATF4q0.25  | 0.001| 0.001|   0.477| 0.862| 1.82|    8|Uchl1 , Ankrd1, Myl4  , Tnnt2 , Tnni1 , Igf2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|ATF4p0.05  | 0.000| 0.000|   0.557| 0.401| 1.53|  410|Uchl1        , Sln          , Nefm         , Psat1        , Krt18        , Ankrd1       , Pcbd1        , Grin3a       , Myl4         , Tnnt2        , Mustn1       , Txndc2       , Fam171b      , Zfp697       , Casq2        , C1qtnf3      , Slc7a5       , Tubb3        , Skint7       , Rab15        , Olfr1368     , Ncam1        , Tnni1        , Gck          , Nt5c2        , Gpr3         , Tnnc1        , Asns         , Chrnd        , Hspa2        , Carhsp1      , Ptpn5        , Cyb5r3       , Cd160        , Grb10        , Tal2         , Bbox1        , Kcnv2        , Maged2       , Nes          , Dclk1        , Fgf17        , Myo1h        , Uchl4        , Epor         , Slc7a3       , Nr2e3        , Olfr550      , Hspb7        , Slpi         , Pln          , 6430548M08Rik, Htatip2      , Kcnn3        , Tmem159      , Tmem37       , Slc7a11      , Lrrn1        , Cib1         , Igf2         , Xirp1        , Csrp3        , Iah1         , Gamt         , Ch25h        , Rassf3       , Gadd45a      , Lypd1        , Pomc         , Chrna1       , Dhrs9        , Rnf113a2     , Apol9b       , Slc16a12     , Rnd2         , Cdkn1a       , Myog         , Rassf5       , Dcun1d3      , Rgs9bp       , Zcchc3       , Napepld      , Cth          , Fetub        , Pdia6        , Cpxm2        , Hspa4l       , Ift122       , Igtp         , Usp11        , Slitrk6      , Peg12        , Mgst3        , Lyzl6        , Trem3        , Gars         , Meox1        , Tubb6        , Snap47       , Crybb3       , Ormdl2       , Ttf2         , Ahcy         , 2310002L09Rik, Apoe         , Zdhhc23      , Dusp6        , Cebpa        , Uck2         , Bin3         , Aven         , Dennd2d      , Rcan3        , Tceal5       , Angpt1       , Cdca4        , Saa3         , Col12a1      , Fgd2         , Pnmt         , Lysmd3       , Hmgb3 |

```r
plotEnrichment(atf4.gene.lists[["ATF4p0.001"]],
               !(is.na(sort(mtsc.ranks, decreasing=T))))

plotGseaTable(atf4.gene.lists, !(is.na(sort(mtsc.ranks,decreasing=T))), mtsc.atf4.pathways, 
              gseaParam=0.5)
```

![](figures/gene-set-atf4-1.png)<!-- -->

## Pathway Analysis of Dual TSC/ATF4 Genes

The co-regulated genes in TSC and ATF4 muscles were:

*Upregulated*: Tnni1, Myl4, Tnnt2, Ankrd1, Igf2, Uchl1, Apoe, Asns, Chrna1, Slc7a5, Krt18, Sln, Nes, Csrp3, Tmem37, Ncam1, Casq2, Kcnn3, Tnnc1, Chrnd, Cyb5r3, Cth, Zfp697, Gadd45a, Maged2, Hspb7, Txndc2, Carhsp1, Psat1, Chrnb1, Myog, Dclk1, Lypd1, Grb10, Gars, Cib1, Kcnv2, Pcbd1, Uck2, Gamt, Mustn1, Rnd2, Rassf3, Dusp6, C1qtnf3, Fam171b, Bin3, Rab15, Hspa2, Tal2, Lrrn1, Xirp1, Cdca4, Nt5c2, Snap47, Iah1, Ift122, Tubb6, Hspa4l, Nefm, Tubb3, Usp11, 6430548M08Rik, Meox1, Ttf2, Nol3, Htatip2, Gck, Tceal5

*Downregulated*:Asb14, Rbmx, Tubd1, Ppp1r1a, Serpine2, Bcap29, Ppp3cb, Ostn, Dcun1d2, Tmem56, Gbp4, March7, Ankrd46, Rgs5, Hnrnpa2b1, Eif3f, Smyd2, Rpap3, Tek, Hook1, Sdr39u1, Lpin1, Mapk9, Strbp, Chd7, Pptc7, Cdkal1, Sec62, Rpl21, Car14, Atp8a1, Usp10, Fmo2, Gpd2, Cd24a, Scn4b, Trim54, Fam118b, Zfp62, Gbp2, Cep350, Art3, Stat1, 4430402I18Rik, Zfp458, Eif1, Hist2h2be, Nub1, Cmbl, Fam161a, Chad, Kcnf1, Bmi1, Ube2b, Ddx5, Ahctf1, Ccdc91



```r
enrich.bp <- enrichGO(c(sig.overlap.up,sig.overlap.down), 
             ont ="BP", 
             keyType = "SYMBOL", 
             OrgDb = org.Mm.eg.db)
library(enrichplot)
#upsetplot(enrich.bp)
#dotplot(enrich.bp, orderBy="NES")
#enrich.bp <- pairwise_termsim(enrich.bp, method="JC")
# emapplot(enrich.bp,
#          cex_label_category=0.5,
#          showCategory = 50,
#          color="NES",
#          min_edge=0.1)

enrich.bp %>% as.data.frame %>% dplyr::select(GeneRatio,pvalue,p.adjust,Description) %>%  filter(p.adjust<0.05) %>% kable(caption="Significant GO-BP pathways for co-regulated genes")
```



Table: Significant GO-BP pathways for co-regulated genes

|           |GeneRatio | pvalue| p.adjust|Description                                                                             |
|:----------|:---------|------:|--------:|:---------------------------------------------------------------------------------------|
|GO:0006941 |10/121    |  0.000|    0.000|striated muscle contraction                                                             |
|GO:0003012 |13/121    |  0.000|    0.000|muscle system process                                                                   |
|GO:0014706 |13/121    |  0.000|    0.000|striated muscle tissue development                                                      |
|GO:0060537 |13/121    |  0.000|    0.001|muscle tissue development                                                               |
|GO:0003009 |5/121     |  0.000|    0.001|skeletal muscle contraction                                                             |
|GO:0007517 |12/121    |  0.000|    0.001|muscle organ development                                                                |
|GO:0006936 |10/121    |  0.000|    0.001|muscle contraction                                                                      |
|GO:0060048 |7/121     |  0.000|    0.001|cardiac muscle contraction                                                              |
|GO:0055002 |8/121     |  0.000|    0.002|striated muscle cell development                                                        |
|GO:0050879 |5/121     |  0.000|    0.002|multicellular organismal movement                                                       |
|GO:0050881 |5/121     |  0.000|    0.002|musculoskeletal movement                                                                |
|GO:0055001 |8/121     |  0.000|    0.002|muscle cell development                                                                 |
|GO:0070296 |4/121     |  0.000|    0.003|sarcoplasmic reticulum calcium ion transport                                            |
|GO:0055008 |5/121     |  0.000|    0.005|cardiac muscle tissue morphogenesis                                                     |
|GO:0051146 |9/121     |  0.000|    0.006|striated muscle cell differentiation                                                    |
|GO:0042692 |10/121    |  0.000|    0.006|muscle cell differentiation                                                             |
|GO:0048738 |8/121     |  0.000|    0.008|cardiac muscle tissue development                                                       |
|GO:0042493 |8/121     |  0.000|    0.008|response to drug                                                                        |
|GO:0060415 |5/121     |  0.000|    0.009|muscle tissue morphogenesis                                                             |
|GO:0032869 |7/121     |  0.000|    0.010|cellular response to insulin stimulus                                                   |
|GO:0042770 |5/121     |  0.000|    0.011|signal transduction in response to DNA damage                                           |
|GO:0045214 |4/121     |  0.000|    0.011|sarcomere organization                                                                  |
|GO:0048644 |5/121     |  0.000|    0.011|muscle organ morphogenesis                                                              |
|GO:0010880 |3/121     |  0.000|    0.012|regulation of release of sequestered calcium ion into cytosol by sarcoplasmic reticulum |
|GO:0010677 |4/121     |  0.000|    0.012|negative regulation of cellular carbohydrate metabolic process                          |
|GO:0060047 |7/121     |  0.000|    0.014|heart contraction                                                                       |
|GO:0055010 |4/121     |  0.000|    0.014|ventricular cardiac muscle tissue morphogenesis                                         |
|GO:0014808 |3/121     |  0.000|    0.014|release of sequestered calcium ion into cytosol by sarcoplasmic reticulum               |
|GO:0045912 |4/121     |  0.000|    0.015|negative regulation of carbohydrate metabolic process                                   |
|GO:1903514 |3/121     |  0.000|    0.015|release of sequestered calcium ion into cytosol by endoplasmic reticulum                |
|GO:0003015 |7/121     |  0.000|    0.015|heart process                                                                           |
|GO:0032868 |7/121     |  0.000|    0.016|response to insulin                                                                     |
|GO:0043462 |4/121     |  0.000|    0.017|regulation of ATPase activity                                                           |
|GO:0071356 |6/121     |  0.000|    0.017|cellular response to tumor necrosis factor                                              |
|GO:0007274 |3/121     |  0.000|    0.018|neuromuscular synaptic transmission                                                     |
|GO:0090257 |7/121     |  0.000|    0.018|regulation of muscle system process                                                     |
|GO:0003229 |4/121     |  0.000|    0.018|ventricular cardiac muscle tissue development                                           |
|GO:0030330 |4/121     |  0.000|    0.018|DNA damage response, signal transduction by p53 class mediator                          |
|GO:0071375 |7/121     |  0.000|    0.021|cellular response to peptide hormone stimulus                                           |
|GO:0030239 |4/121     |  0.000|    0.024|myofibril assembly                                                                      |
|GO:0035094 |3/121     |  0.000|    0.025|response to nicotine                                                                    |
|GO:0034612 |6/121     |  0.000|    0.025|response to tumor necrosis factor                                                       |
|GO:0005979 |3/121     |  0.000|    0.026|regulation of glycogen biosynthetic process                                             |
|GO:0010962 |3/121     |  0.000|    0.026|regulation of glucan biosynthetic process                                               |
|GO:0007519 |6/121     |  0.001|    0.028|skeletal muscle tissue development                                                      |
|GO:0062012 |8/121     |  0.001|    0.028|regulation of small molecule metabolic process                                          |
|GO:0003007 |7/121     |  0.001|    0.030|heart morphogenesis                                                                     |
|GO:0048747 |4/121     |  0.001|    0.030|muscle fiber development                                                                |
|GO:0060538 |6/121     |  0.001|    0.032|skeletal muscle organ development                                                       |
|GO:0003208 |4/121     |  0.001|    0.032|cardiac ventricle morphogenesis                                                         |
|GO:0005977 |4/121     |  0.001|    0.032|glycogen metabolic process                                                              |
|GO:0006073 |4/121     |  0.001|    0.032|cellular glucan metabolic process                                                       |
|GO:0044042 |4/121     |  0.001|    0.032|glucan metabolic process                                                                |
|GO:0043516 |3/121     |  0.001|    0.033|regulation of DNA damage response, signal transduction by p53 class mediator            |
|GO:0042326 |9/121     |  0.001|    0.033|negative regulation of phosphorylation                                                  |
|GO:0002931 |3/121     |  0.001|    0.037|response to ischemia                                                                    |
|GO:0070873 |3/121     |  0.001|    0.037|regulation of glycogen metabolic process                                                |
|GO:1901653 |7/121     |  0.001|    0.038|cellular response to peptide                                                            |
|GO:0003206 |5/121     |  0.001|    0.039|cardiac chamber morphogenesis                                                           |
|GO:0034765 |9/121     |  0.001|    0.039|regulation of ion transmembrane transport                                               |
|GO:0032781 |3/121     |  0.001|    0.040|positive regulation of ATPase activity                                                  |
|GO:0032885 |3/121     |  0.001|    0.040|regulation of polysaccharide biosynthetic process                                       |
|GO:0006112 |4/121     |  0.001|    0.040|energy reserve metabolic process                                                        |
|GO:0014883 |2/121     |  0.001|    0.040|transition between fast and slow fiber                                                  |
|GO:0030049 |2/121     |  0.001|    0.040|muscle filament sliding                                                                 |
|GO:2001138 |2/121     |  0.001|    0.040|regulation of phospholipid transport                                                    |
|GO:2001140 |2/121     |  0.001|    0.040|positive regulation of phospholipid transport                                           |
|GO:0048638 |8/121     |  0.001|    0.045|regulation of developmental growth                                                      |
|GO:0045445 |4/121     |  0.001|    0.046|myoblast differentiation                                                                |
|GO:0014889 |2/121     |  0.001|    0.046|muscle atrophy                                                                          |
|GO:0048630 |2/121     |  0.001|    0.046|skeletal muscle tissue growth                                                           |
|GO:0043535 |4/121     |  0.001|    0.046|regulation of blood vessel endothelial cell migration                                   |
|GO:0070252 |4/121     |  0.001|    0.046|actin-mediated cell contraction                                                         |
|GO:0005978 |3/121     |  0.002|    0.047|glycogen biosynthetic process                                                           |
|GO:0009250 |3/121     |  0.002|    0.047|glucan biosynthetic process                                                             |
|GO:0010675 |5/121     |  0.002|    0.048|regulation of cellular carbohydrate metabolic process                                   |
|GO:1903169 |5/121     |  0.002|    0.048|regulation of calcium ion transmembrane transport                                       |
|GO:0014874 |2/121     |  0.002|    0.048|response to stimulus involved in regulation of muscle adaptation                        |
|GO:0030949 |2/121     |  0.002|    0.048|positive regulation of vascular endothelial growth factor receptor signaling pathway    |
|GO:0033275 |2/121     |  0.002|    0.048|actin-myosin filament sliding                                                           |
|GO:0032881 |3/121     |  0.002|    0.048|regulation of polysaccharide metabolic process                                          |
|GO:0010959 |8/121     |  0.002|    0.048|regulation of metal ion transport                                                       |
|GO:0010594 |5/121     |  0.002|    0.048|regulation of endothelial cell migration                                                |
|GO:0043502 |4/121     |  0.002|    0.048|regulation of muscle adaptation                                                         |
|GO:0062014 |4/121     |  0.002|    0.048|negative regulation of small molecule metabolic process                                 |

# Specific Genes

## Sarcolipin


```r
sln.id <- tT2 %>% filter(Gene.symbol=='Sln') %>% rownames()
sln.exprs <- exprs(gset)[sln.id,]

sln.data <- data.frame(sln.exprs, gs)

library(forcats)
sln.data %>%
  group_by(gs) %>%
  summarize(avg=mean(2^sln.exprs),
            error=se(2^sln.exprs)) %>%
  mutate(norm.mean = avg/avg[1],
         norm.error = error/avg[1]) %>%
  ggplot(aes(x=gs,
             y=norm.mean, 
             ymin=norm.mean-norm.error,
             ymax=norm.mean+norm.error)) +
  geom_bar(stat='identity', position='dodge') +
  geom_errorbar(width=0.5) +
  scale_x_discrete(labels=c("Control" = "Empty Vector", "ATF4.OE" = "ATF4")) +
  theme_classic() +
  theme(text=element_text(size=20)) +
  labs(y='Relative mRNA Expression',
       x='',
       title='Effects ATF4 Overexpression on Sarcolipin')
```

![](figures/atf4-sln-barplot-1.png)<!-- -->

## Ketolytic Genes

```{atf4-ketolytic-barplot}
genes.to.test <- c('Slc16a1','Bdh1','Hmgcs2')
kt.id <- tT2 %>% filter(Gene.symbol%in%genes.to.test) %>% rownames()
kt.exprs <- exprs(gset)[kt.id,] 
rownames(kt.exprs) <- genes.to.test

library(tibble)
kt.data <- as.data.frame(kt.exprs) %>%
  rownames_to_column(var="Gene") %>%
  pivot_longer(-Gene, names_to='Sample') %>%
  mutate(Treatment=rep(as.character(gs),3))


library(forcats)
kt.data %>%
  group_by(Gene,Treatment) %>%
  summarize(avg=mean(2^sln.exprs),
            error=se(2^sln.exprs),
            .groups='keep') %>% #not working here
  mutate(norm.mean = avg/avg[Treatment=='Control'],
         norm.error = error/avg[Treatment=='Control']) %>%
  ggplot(aes(x=Treatment,
             y=avg
             ymin=norm.mean-norm.error,
             ymax=norm.mean+norm.error)) +
  geom_bar(stat='identity', position='dodge') +
  geom_errorbar(width=0.5) +
  scale_x_discrete(labels=c("Control" = "Empty Vector", "ATF4.OE" = "ATF4")) +
  theme_classic() +
  theme(text=element_text(size=16)) +
  labs(y='Relative mRNA Expression',
       x='',
       title='Effects ATF4 Overexpression on Ketolysis Genes')
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
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] forcats_0.5.1          fgsea_1.16.0           enrichplot_1.10.2     
##  [4] org.Mm.eg.db_3.12.0    AnnotationDbi_1.52.0   IRanges_2.24.1        
##  [7] S4Vectors_0.28.1       clusterProfiler_3.18.1 venneuler_1.1-0       
## [10] rJava_1.0-6            ggplot2_3.3.5          maptools_1.1-2        
## [13] sp_1.4-6               umap_0.2.7.0           limma_3.46.0          
## [16] GEOquery_2.58.0        Biobase_2.50.0         BiocGenerics_0.36.1   
## [19] broom_0.7.11           dplyr_1.0.7            tidyr_1.1.4           
## [22] knitr_1.37            
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_2.0-2    ellipsis_0.3.2      qvalue_2.22.0      
##   [4] rstudioapi_0.13     farver_2.1.0        graphlayouts_0.8.0 
##   [7] ggrepel_0.9.1       bit64_4.0.5         RSpectra_0.16-0    
##  [10] fansi_1.0.0         scatterpie_0.1.7    xml2_1.3.3         
##  [13] splines_4.0.2       cachem_1.0.6        GOSemSim_2.16.1    
##  [16] polyclip_1.10-0     jsonlite_1.7.2      GO.db_3.12.1       
##  [19] png_0.1-7           ggforce_0.3.3       BiocManager_1.30.16
##  [22] readr_2.1.1         compiler_4.0.2      rvcheck_0.2.1      
##  [25] backports_1.4.1     assertthat_0.2.1    Matrix_1.4-0       
##  [28] fastmap_1.1.0       cli_3.1.0           tweenr_1.0.2       
##  [31] htmltools_0.5.2     tools_4.0.2         igraph_1.2.6       
##  [34] gtable_0.3.0        glue_1.6.0          reshape2_1.4.4     
##  [37] DO.db_2.9           fastmatch_1.1-3     Rcpp_1.0.7         
##  [40] jquerylib_0.1.4     vctrs_0.3.8         nlme_3.1-153       
##  [43] ggraph_2.0.5        xfun_0.29           stringr_1.4.0      
##  [46] lifecycle_1.0.1     DOSE_3.16.0         MASS_7.3-54        
##  [49] scales_1.1.1        tidygraph_1.2.0     vroom_1.5.7        
##  [52] hms_1.1.1           ggupset_0.3.0       RColorBrewer_1.1-2 
##  [55] yaml_2.2.1          curl_4.3.2          memoise_2.0.1      
##  [58] reticulate_1.22     gridExtra_2.3       downloader_0.4     
##  [61] ggfun_0.0.4         yulab.utils_0.0.4   sass_0.4.0         
##  [64] stringi_1.7.6       RSQLite_2.2.9       highr_0.9          
##  [67] BiocParallel_1.24.1 rlang_0.4.12        pkgconfig_2.0.3    
##  [70] evaluate_0.14       lattice_0.20-45     purrr_0.3.4        
##  [73] labeling_0.4.2      cowplot_1.1.1       shadowtext_0.1.1   
##  [76] bit_4.0.4           tidyselect_1.1.1    plyr_1.8.6         
##  [79] magrittr_2.0.1      R6_2.5.1            generics_0.1.1     
##  [82] DBI_1.1.2           pillar_1.6.4        foreign_0.8-81     
##  [85] withr_2.4.3         mgcv_1.8-38         tibble_3.1.6       
##  [88] crayon_1.4.2        utf8_1.2.2          tzdb_0.2.0         
##  [91] rmarkdown_2.11      viridis_0.6.2       grid_4.0.2         
##  [94] data.table_1.14.2   blob_1.2.2          digest_0.6.29      
##  [97] openssl_1.4.6       munsell_0.5.0       viridisLite_0.4.0  
## [100] bslib_0.3.1         askpass_1.1
```
