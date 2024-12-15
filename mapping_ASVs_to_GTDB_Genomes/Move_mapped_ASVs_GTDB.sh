#!/bin/bash
#Author: Ioannis Kampouris
#Purpose: To move aligned genomes to our ASVs.
#Dependences use BASH
 
echo "Start proscessing"

mkdir s_Genomes
#Create an alignment file

filename='/vol/IDK2/ms_model_ASVs/Genomes_aligned/selected_list.csv'


while read i;
  do 

 	M=$(echo ${i} | sed "s/fna.gz/fna/g")
	S=$(echo ${i} | sed "s/\\//_/g")
	F=$(echo ${i} | sed "s/fna/fna\\.gz/g")



	cp  /vol/IDK2/release214/"$i" /vol/IDK2/release214/s_Genomes/"$S"
	cp  /vol/IDK2/release214/"$F" /vol/IDK2/release214/s_Genomes/"$S".gz
	# /vol/IDK2/prokka/bin/prokka    --outdir s_Genomes/"$S" "$S"  --cpus 24 --force --prefix "$S"

	echo "$F"
	#echo "/vol/IDK2/release214/${M}"
	#gzip "/vol/IDK2/release214/${M}"
done<"$filename"
					    	                       
echo "All complete"

					    	                       							