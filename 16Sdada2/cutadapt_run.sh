#!/bin/bash/python

#Author: Ioannis Kampouris
#Purpose: Use fastqc and multiqc in the R Server (UBUNTU OS)
#Dependencies a path for fastqc and 
# Primer Sequences According to NOVOGENE
# FW: CCTAYGGGRBGCASCAG RE:GGACTACNNGGGTATCTAAT

cd '2020_DiControl_AP4-2Mais'

filename='s_2020_sample_list.txt'
echo "Start 2020 samples" 
mkdir trimmed_cutadapt
mkdir
while read i;
   do 
   SAMPLE=$(echo ${i} | sed "s/raw_1\.fq\.gz//")
   echo "$SAMPLE"
   do cutadapt -g CCTAYGGGRBGCASCAG -G GGACTACNNGGGTATCTAAT \ 
   --discard-untrimmed -o ${i}_2020_trimmed_cutadapt_1.fq -p ${i}_2020_trimmed_cutadapt_2.fq \
   ${i}raw_1.fq ${i}raw_2.fq  > ${i}cutadapt.log

done 