# RNAseq Analysis Log

## Design Considerations

Used Salmon for transcript mapping followed by txnip and DESeq2 for quantification based on the suggested protocol from the DESeq2 manual https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#tximport. Recommendation is to use --gcBias flag.

## Installations and versioning

Salmon v 1.3.0
Mouse Genome GRCm38.p6 (GCA_000001635.8)


### Installing Salmon

Used conda to install and activate a salmon environment

conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n salmon salmon
conda activate salmon

## Index for Mouse Genome

from http://uswest.ensembl.org/info/data/ftp/index.html

ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/
ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/
ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/cds/ #for the coding sequences

To get the transcriptome files
curl ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz --output Mus_musculus.GRCm38.cds.all.fa.gz


salmon index -t Mus_musculus.GRCm38.cds.all.fa.gz -i mouse_index

#note that the index is not in the git repositiry and has to be remade if you download these data

this file is set as salmon-quant-bridges.sh or salmon-quant-gregg.sh
#!/bin/bash
for fn in 1415-NEH/NovaA-284/1415-NEH/Sample_1415-NEH-{1..11};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
f1=`ls ${fn}/*R1_001.fastq.gz`
f2=`ls ${fn}/*R2_001.fastq.gz`
echo "Processing file `basename ${f1}`"
echo "Processing file `basename ${f2}`"
salmon quant -i mouse_index -l A \
         -1 ${f1} \
         -2 ${f2} \
         -p 4 --validateMappings --gcBias -o quants/NEH-${samp}_quant
done 