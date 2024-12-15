#!/bin/bash
#Author: Ioannis Kampouris
#Purpose: Use cutadapt to remove primers in the R Server (UBUNTU OS)
#Dependencies a path for fastqc and 
# Primer Sequences According to NOVOGENE
# FW: CCTAYGGGRBGCASCAG RE:GGACTACNNGGGTATCTAAT

cd '2020'
find -name  "*raw_1*"|sed 's/^..//'  > s_2020_sample_list.txt #Change accordingly to your path 

filename='s_2020_sample_list.txt'
echo "Start 2020 samples" 
mkdir -p trimmed_cutadapt/first_2020
mkdir reports

while read i;
	do 
   SAMPLE=$(echo ${i} | sed "s/raw_1\.fq\.gz//")
   echo "$SAMPLE"
   cutadapt -g CCTAYGGGRBGCASCAG -G GGACTACNNGGGTATCTAAT  --discard-untrimmed   \
   -o ${SAMPLE}.2020.trimmed_cutadapt_1.fq -p         \
     ${SAMPLE}.2020.trimmed_cutadapt_2.fq    ${SAMPLE}raw_1.fq.gz ${SAMPLE}raw_2.fq.gz \
	 > ${SAMPLE}.log
    
	
	
done <"$filename"

find . -name '*.fq' -exec mv {} trimmed_cutadapt/first_2020 \;
find . -name '*.log' -exec mv {} reports \;

cd ..

cd '2021'
find -name  "*raw_1*"|sed 's/^..//'  > s_2021_sample_list.txt #Change accordingly to your path 

filename='s_2021_sample_list.txt'
echo "Start 2021 samples" 
mkdir -p trimmed_cutadapt/second_2021
mkdir reports

while read i;
	do 
   SAMPLE=$(echo ${i} | sed "s/raw_1\.fq\.gz//")
   echo "$SAMPLE"
   cutadapt -g CCTAYGGGRBGCASCAG -G GGACTACNNGGGTATCTAAT  --discard-untrimmed    \
   -o ${SAMPLE}.2021.trimmed_cutadapt_1.fq -p         \
     ${SAMPLE}.2021.trimmed_cutadapt_2.fq    ${SAMPLE}raw_1.fq.gz ${SAMPLE}raw_2.fq.gz \
	  > ${SAMPLE}.log
    
	
	
done <"$filename"

cd '2022'
find -name  "*_1*"|sed 's/^..//'  > s_2021_sample_list.txt #Change accordingly to your path 

filename='s_2021_sample_list.txt'
echo "Start 2021 samples" 
mkdir -p trimmed_cutadapt/second_2021
mkdir reports

while read i;
	do 
   SAMPLE=$(echo ${i} | sed "s/raw_1\.fq\.gz//")
   echo "$SAMPLE"
   cutadapt -g CCTAYGGGRBGCASCAG -G GGACTACNNGGGTATCTAAT    \
   -o ${SAMPLE}.2021.trimmed_cutadapt_1.fq -p         \
     ${SAMPLE}.2021.trimmed_cutadapt_2.fq    ${SAMPLE}raw_1.fq.gz ${SAMPLE}raw_2.fq.gz \
	  > ${SAMPLE}.log
    
	
	
done <"$filename"



find . -name '*.fq' -exec mv {} trimmed_cutadapt/second_2021 \;
find . -name '*.log' -exec mv {} reports \;
