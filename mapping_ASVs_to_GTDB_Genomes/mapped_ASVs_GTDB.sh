#/bin/bash
#Author: Ioannis Kampouris
#Purpose: To perform alignment for finding sequences from GTDB and align them to ASVs.
#Dependences use usearch
#For this you need to convert the saved fasta files of Comamonadaceae ASVs to genomes
 
echo "Start proscessing"
#Enter the directory

i=$1
M=$(echo ${i} | sed "s/fna.gz/fna/")
echo "$M" 
gunzip "$i" 
usearch -usearch_local  "$M" -db /vol/IDK2/ms_model_ASVs/ASVs_MS_Models.udb -id 0.97  -strand both  -evalue  1e-10 -threads 27  \
  -userout /vol/IDK2/ms_model_ASVs/Genomes_aligned/$M"_Genome.txt -userfields query+target+evalue+id+alnlen+qs+ts 
 
gzip "$M" 				    	               
echo "$i" >/vol/IDK2/ms_model_ASVs/Genomes_aligned/name 
paste -d ";" /vol/IDK2/ms_model_ASVs/Genomes_aligned/name /vol/IDK2/ms_model_ASVs/Genomes_aligned/$M"_Genome.txt > /vol/IDK2/ms_model_ASVs/Genome2.txt
					    	                 

