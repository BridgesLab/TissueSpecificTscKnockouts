Running R Scripts
===================

Once the counts tables are generated for both gene counts and exon counts, all the other analyses are run using the scripts used here.  All analysis file are R scripts as Rmarkdown (Rmd) files.  These files can be run inside RStudio https://www.rstudio.com/ and will generate the md and html files which include the processed data.  For more information on using R and using these files see http://cran.us.r-project.org/.  The scripts should be run in this order.

## Running the Differential Expression Analyses

We use DESeq2 to do the differential expression analysis.  This takes the htseq-counts file (**../data/raw/htseq-counts**) along with the sample mapping file (**../data/raw/sample_mapping.csv**) and generates the results file (**../data/processed/Binary DESeq Results.csv**).  It also generates the normalized counts file used in other scripts.  This file is found at **../data/processed/Binary Normalized Counts.csv**.  This script also generates the diagnostic dispersion, MA and PCA plots as well as the heatmap.

## Generating the RPKM Counts Table

The normalized counts are not normalized to transcript length.  This normalization is done in the **RPKM_counts_table.Rmd** script.  It generates the output file **../data/processed/RPKM_count.csv** which is now normalized RPKM values.

## Generating barplots of the data

The **barplots.Rmd** script is used to generate the barplots from the RPKM data.  These barplots are automatically saved in the scripts/figures folder.

## Performing GSEA Analysis

This is done separately using a standalone program, see the data/processed folder for more details.  To generate the input files for the preranked GSEA analysis the script **GSEA_inputs.Rmd** takes the DESeq results file and puts it in the correct format.  To parse the results, a python script called **MakeGSEATable.py** takes the gsea output and created the summary files.