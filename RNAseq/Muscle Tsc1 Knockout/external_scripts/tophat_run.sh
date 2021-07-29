#!/bin/bash
#to use this script, the following links need to be established:
#reference-annotations folder is symlinked

#this is the GTF file for the assembly
GTF="reference-annotations/Mus_musculus.GRCm38.75.gtf"
#this is the name (without the .fa) of the index built with bowtie2-build
REFERENCE="reference-annotations/Mus_musculus.GRCm38.75.dna.toplevel"

#this is where the csfasta and qual files are located
FASTA="fasta_files/Libraries"
LOCATION="/mnt/sdc1/lifescopeusers/dbridge/bridges"
echo "This alignment is using bowtie with the $REFERENCE genome"
echo "The alignments and assemblies are made using $GTF"

#remove existing alignments and assemblies
rm -r tophat_out

mkdir tophat_out


for d in $FASTA/*
do
  sample=`basename $d`
  #run tophat alignment
  #this initializes the output directories to avoid a filesystem detection problem in glusterfs
  echo "stat tophat_out" > $sample.sh
  #this uses multiple processors (-p 11) and to first map to known sequences (-G) before matching other sequences.  The lanes are combined by commas.
  echo "tophat2 -p 11  --bowtie1 --color --qual -G $GTF -o tophat_out/$sample $REFERENCE $FASTA/${sample}/F3/reads/L01_${sample}_F3.csfasta,$FASTA/${sample}/F3/reads/L02_${sample}_F3.csfasta,$FASTA/${sample}/F3/reads/L03_${sample}_F3.csfasta  $FASTA/${sample}/F3/reads/L01_${sample}_F3.QV.qual,$FASTA/${sample}/F3/reads/L02_${sample}_F3.QV.qual,$FASTA/${sample}/F3/reads/L03_${sample}_F3.QV.qual" >> $sample.sh
  qsub -d $LOCATION -h $sample.sh
  echo "Submitting sample $sample"
  #rm $sample.sh
done
