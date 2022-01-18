Some of these analyses are done externally on a high performance computing cluster.  The data is then taken out of the cluster and processed locally.  

Processing done on the high performance computing cluster
-----------------------------------------------------------

The externally run scripts are found in the **external_scripts** folder.  These extract the xsq files, aligns them to the mouse genome and then generates the counts tables.


### File Preparaton

Prepared csfasta and quals files using the XSQ_Tools package on the XSQ files.

First downloaded the xsq files from the MRC server
Used convertFromXSQ script with xsqconvert version -- 1.5

```bash
qsub -d /mnt/sdc1/lifescopeusers/dbridge/bridges ./xsq_convert_L01.sh
qsub -d /mnt/sdc1/lifescopeusers/dbridge/bridges ./xsq_convert_L02.sh
qsub -d /mnt/sdc1/lifescopeusers/dbridge/bridges ./xsq_convert_L03.sh
```

### Alignments

Ran all the alignments using the colorspace options of bowtie 1 and against the Mus_musculus.GRCm38.75 genome using the script tophat_run.sh
Each run had to be manually submitted and directed to a free node using qsub -l nodes=XXXXX.  The software versions were:

* Used Samtools version 0.1.19-44428cd
* Used TopHat v2.0.10
* Used bowtie version 1.0.0


### Counts Table Generation

Generated exon and transcript counts tables using HTSeq version 0.5.4p5:

* HTSeq-counts-table.sh
* HTSeq-counts-table-exons.sh

Locally run analyses
----------------------

The two counts tables and the sample mapping files are the raw data for the rest of the differential expression analyses.

### Raw Data

All the raw data is in the **data/raw** folder.  This data is automatically obtained from our internal LIMS system and can only be updated from that point.  The script files use this raw data to do the analysis.  Data which has been processed is saved in the **data/processed** folder.

### Script files

This code base includes the raw data and reproducible R code for these analyses within the **scripts** folder.  A README file in that folder describes these analyses.  All analysis file are R scripts as Rmarkdown (Rmd) files.  These files can be run inside RStudio https://www.rstudio.com/ and will generate the md and html files which include the processed data.  For more information on using R and using these files see http://cran.us.r-project.org/.
