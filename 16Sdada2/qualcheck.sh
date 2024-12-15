#!/bin/bash
#Author: Ioannis Kampouris
#Purpose: Use fastqc and multiqc in the R Server (UBUNTU OS)
#Dependencies a path for fastqc and 

export PATH=$PATH:"/home/ioannis.kampouris/fastqc/FastQC" #Change accordingly to your path 
export PATH=$PATH:"/home/ioannis.kampouris/fastqc/MultiQC" #Change accordingly to your path 

cd '2020_DiControl_AP4-2Mais'
find -name  "*raw_1*"|sed 's/^..//'  > s_2020_sample_list.txt #Change accordingly to your path 

mkdir QC

filename='s_2020_sample_list.txt'
echo "Start 2020 samples" 
mkdir 
while read i;
   do 
   SAMPLE=$(echo ${i} | sed "s/raw_1\.fq\.gz//")
      echo "$SAMPLE"
   fastqc  ${SAMPLE}raw_2.fq.gz -o QC #You can remove the "RAW" if you want your files with NOVOGENE cleaning
    fastqc  ${SAMPLE}raw_1.fq.gz -o QC
     
done <"$filename"
 
 mkdir QC/FW
 mkdir QC/RE

mv QC/*raw_1* QC/FW
mv QC/*raw_2* QC/RE
cd QC/FW
python3 -m multiqc  .
cd ..
cd RE
python3 -m multiqc  .

cd ..
cd ..
cd ..

cd '2021_DiControl_AP4_WdhMais'
find -name  "*raw_1*"|sed 's/^..//'  > s_2021_sample_list.txt #Change accordingly to your path 
mkdir QC

filename='s_2021_sample_list.txt'
echo "Start 2021 samples" 
 
while read i;
   do 
   SAMPLE=$(echo ${i} | sed "s/raw_1\.fq\.gz//")
      echo "$SAMPLE"
   fastqc  ${SAMPLE}raw_2.fq.gz -o QC #You can remove the "RAW" if you want your files with NOVOGENE cleaning
    fastqc  ${SAMPLE}raw_1.fq.gz -o QC
     
done <"$filename"
 
 mkdir QC/FW
 mkdir QC/RE

mv QC/*raw_1* QC/FW
mv QC/*raw_2* QC/RE
cd QC/FW
python3 -m multiqc  .
cd ..
cd RE
python3 -m multiqc  .