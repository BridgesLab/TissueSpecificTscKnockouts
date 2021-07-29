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
