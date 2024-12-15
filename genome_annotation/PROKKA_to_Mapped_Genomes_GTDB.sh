#!/bin/bash
#Author: Ioannis Kampouris
#Purpose: To annotate aligned genomes to our ASVs.
#Dependences use PROKKA
 
echo "Start proscessing"


#Create an alignment file
cd /vol/IDK2/release214/s_Genomes
mkdir annotations

for i in *.fna;
  do 


   
    /vol/IDK2/prokka/bin/prokka   "$i"  --outdir annotations/"$i"   --cpus 24 --force 

	echo "$F"

done

echo "All complete"

					    	                       							