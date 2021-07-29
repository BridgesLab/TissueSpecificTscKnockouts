#!/bin/bash

#This script cycles through all accepted_hits.bam files in a directory then runs it through htseq-count
#The name of the file is taken in and used for the output files
GTF="reference-annotations/Mus_musculus.GRCm38.75.gtf"
TOPHAT_OUT="tophat_out"
OUTPUT_DIR="htseq-counts"

#clear output directory
rm -r $OUTPUT_DIR
mkdir $OUTPUT_DIR

echo "Using reference genome file in $GTF"
echo "Searching for accepted_hits.bam files in $TOPHAT_OUT"


for directory in `ls $TOPHAT_OUT` ; do
    #then write a script which pipes the samtools view output (which makes the sam file) then makes a file called counts
    echo "Processing $directory"
    samtools view $TOPHAT_OUT/$directory/accepted_hits.bam |  htseq-count --stranded=no  --mode=intersection-nonempty - $GTF > $OUTPUT_DIR/$directory-counts.txt
done
