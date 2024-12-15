library(ggplot2) #ggplot2: Fancy graphs
library(ggpubr) #ggplot expansion for more fancy graphs
library(readxl)
library(readr)
library(dplyr)# Data handling
library(vegan)# Multivariate stat tool
library(tidyr) # Data handling
library(RColorBrewer) # Colours for fancy graphs
library(tibble)# Data handling: rownames/columnanmes
library(phyloseq) # tool for 16S amplicon analysis
library(reshape2)  # Data handling
library(dunn.test)# non-parametric post-hoc test
library(ape) # tool for phylogenetic trees
library(phangorn) #tool for phylogenetic trees
library(stringr)# string manipulation and text modification
library(plyr)
library(ARTool)
library(forcats)
library(tidyverse)
library(pairwiseAdonis)

#Script for alpha and beta-diversity for Rhizosphere and 
# Soil 16S rRNA amplicon: Dataset Bernburg 2020
# Assume that you have your DNA concentrations or the Copies per qpcr
ASVtable_working <- read.csv("dada2_out/ASV_Table_LTE_2020_2021.csv")%>%select(ASV, everything())
DNA_NG=readxl::read_excel("dada2_out/Di_Control_DNA_2020_2021.xlsx", 
sheet = "Dcont")
#Load the available 16S rRNA gene copies measured for several samples
Copies_NG=readxl::read_excel("dada2_out/Di_Control_DNA_2020_2021.xlsx", 
                          sheet = "Dcont2")%>%na.omit() %>%aggregate(Copies~Sample, mean)
Copies_DNA= full_join(DNA_NG, Copies_NG, by="Sample")
#In this dataset MW is  concentration and einwage the weight. 
set.seed(17011990)

#Function to normalize the ASV table 
normalise_100 <- function(x){100*(x/sum(x))}

nrow(ASVtable_working%>%
       filter(., str_detect(tax.Kingdom, "Archaea")))

#33 Archaea in our ASV table
ASVtax_table <- ASVtable_working[,1:9] %>%as.data.frame()
ASVtax_table <- ASVtax_table %>% mutate(Taxonomy=paste0(tax.Kingdom,";",tax.Phylum,";",tax.Class,";",tax.Order,";",tax.Family,";",tax.Genus))%>%
  mutate(Taxonomy =gsub("[[:space:]]", "_", Taxonomy))%>%
  mutate(Taxonomy=str_replace_all(Taxonomy, ".ncultured.*.", "Unclassified"))%>%
  mutate(Taxonomy=str_replace_all(Taxonomy, "uncultured*", "Unclassified"))%>%
  mutate(Taxonomy=str_replace(Taxonomy, "metagenome", "Unclassified"))

ASVtax_table <- select(ASVtax_table,ASV,Sequence ,Taxonomy)%>%
  separate(Taxonomy,
           into = c("Kingdom","Phylum","Class","Order","Family","Genus"), sep=";")
ASVtax_table[ASVtax_table=="NA"]="Unclassified"
ASVtax_table[ASVtax_table=="uncultured"]<-"Unclassified"
ASVtax_table[ASVtax_table=="metagenome"]<-"Unclassified"
rownames(ASVtax_table)=ASVtax_table$ASV
#Assign metadata via sample names
metadata_assign=data.frame(Sample=
colnames(ASVtable_working[,16:ncol(ASVtable_working)]), 
ID=colnames(ASVtable_working[,16:ncol(ASVtable_working)]))%>%
  mutate(Sample=str_replace(Sample,".trimmed",""))%>%
  mutate(Sample=str_replace(Sample,".trimmed",""))%>%
  separate(Sample, into = c("Sample", "Year"), sep = "..20" )%>%
  mutate(Year= paste0("20", Year))%>%
  
  mutate(Type= str_replace(Sample, ".*.T0.*.","FS_T0"))%>%
  
  mutate(Type= str_replace(Type, "FS.*.","Root-Associated_Soil"))%>%
  mutate(Type= str_replace(Type, "RH.*","Rhizosphere"))%>%
  mutate(Tillage= str_replace(Sample, ".*CT.*.","CT"))%>%
  mutate(Tillage= str_replace(Tillage, ".*MP.*.","MP"))%>%
  mutate(Int= str_replace(Sample, ".*In.*.","Intensive"))%>%
  mutate(Int= str_replace(Int, ".*Ex.*.","Extensive"))%>%
  mutate(BMc= str_replace(Sample, ".*.t.I.*","BMc"))%>%
  mutate(BMc= str_replace(BMc, ".*.x.I.*","BMc"))%>%
  mutate(BMc= str_replace(BMc, ".*.n.I.*","BMc"))%>%
  mutate(Int= str_replace(Int, ".*Ex.*.","Extensive"))%>%
  mutate(Block=str_sub(Sample, -1, -1)  )%>% 
  mutate(Root_Window=str_replace(Sample, ".*.RW.*.","RW"))%>%filter(!Root_Window=="RW")
metadata_assign$BMc[metadata_assign$BMc!="BMc"]<-"Ctrl"

metadata_assign$Year=as.numeric(metadata_assign$Year)
metadata_assign$Block=as.numeric(metadata_assign$Block)


metadata_DNA= full_join(Copies_DNA%>%select(-Sample), metadata_assign, by=colnames(metadata_assign%>%select(-ID,-Sample,-Root_Window)))%>%
  filter(Type=="Rhizosphere")

#Load the Taxonomy Table of the Rarefied ASV table. 
rownames(ASVtax_table)=ASVtax_table$ASV
#Assign metadata via sample names 
#If you have any other metadata file please use it.

ASVtable=dplyr::select(ASVtable_working, ASV,
c(metadata_DNA$ID))%>%column_to_rownames(var="ASV")


#I will train my model on the with the Rhizosphere, because I lack for the water content of rhizosphere

ASVtable2= ASVtable
ASVtable2=apply(ASVtable2+1,2,function (x) 100*(x/sum(x)))%>%na.omit()


correlations=function(x){
file2= cor.test(x, 100*(metadata_DNA$MW)/metadata_DNA$Weight, method = "spearman")
file3=data.frame( p=file2$p.value, rho=file2$estimate)
return(file3)
}
#Apply the spearman correlations in parrallel for all ASVs
#Use your sample name order from metadata (DNA/Gene copies concentration)  
#to ensure you have the correct order

spearman_correlations= apply(ASVtable[,metadata_DNA$ID],1, function(x) correlations(x))
#Create the data.frame from the list.
spearman_correlations_df <- as.data.frame(do.call(rbind, spearman_correlations))
#Find the contaminatns
contaminants=filter(spearman_correlations_df, p<0.05&rho<0)
#Check your contaminants
View(contaminants)
rownames(metadata_DNA)<-metadata_DNA$ID
#After selecting the dataset for training we will work with the whole datase
ASVtable3= dplyr::select(ASVtable_working, ASV,
                         c(metadata_DNA$ID))%>%column_to_rownames(var="ASV")
ASVtable3=apply(ASVtable3,2, function (x) 100*(x)/sum(x))
#Get the Relative Abundance of contaminants
contaminants_RA= data.frame(contaminants= colSums(ASVtable3[rownames(contaminants),]))
distribution_of_contaminants= ggplot(contaminants_RA, aes(x=contaminants)) + geom_bar(stat = "density")+
  ggtitle(paste0("Mean: ", round(mean(contaminants_RA$contaminants),2), " Standard Deviation: ",
                 round(sd(contaminants_RA$contaminants),2))) + xlab("Contaminants RA (%)")
distribution_of_contaminants

contaminant_free=filter(ASVtable_working, !ASV%in%c(rownames(contaminants)))

ASV_counts=data.frame(Counts=rowSums(contaminant_free[,metadata_assign$ID]), ASV=contaminant_free$ASV)%>%
  filter(Counts>=10)
contaminant_free2=filter(ASVtable_working%>%select(colnames(ASVtable_working[,1:9]), metadata_assign$ID), ASV%in%c(ASV_counts$ASV))
contamination_info=data.frame( number_of_removed=nrow(contaminant_free2),
                               number_of_ASVs=nrow(ASVtable_working))
contamination_info2=as.data.frame( cbind(colSums(ASVtable_working[,metadata_assign$ID]),
      colSums(contaminant_free2[,metadata_assign$ID])))
colnames(contamination_info2)<-c("Merged", "Contaminant_free")

write.csv(contamination_info, file = "dada2_out/number_of_ASVs_contamination_passed_LTE_2020_2021.csv")
write.csv(contamination_info2, file = "dada2_out/number_of_reads_contamination_passed__LTE_2020_2021.csv")
write.csv(contaminant_free2, file ="dada2_out/plastid_clean_low_read_clean_contaminant_free_ASV_Table_LTE_2020_2021.csv", row.names =T)

file1=read.csv("dada2_out/plastid_clean_low_read_clean_contaminant_free_ASV_Table_LTE_2020_2021.csv")
metadata_assign_RH=filter(metadata_assign, Type=="Rhizosphere")
file2=file1[c("ASV", metadata_assign_RH$ID)]
file2[,2:ncol(file2)]%>%sum()
file2[,2:ncol(file2)]%>%apply(.,2, function(x) sum(x))%>%mean()
file2[,2:ncol(file2)]%>%apply(.,2, function(x) sum(x))%>%sd()
file3= as.data.frame( file2[,2:ncol(file2)]%>%apply(.,1, function(x) 100*sum(x>0)/length(x)))
colnames(file3)="Prev"
file4=filter(file3, Prev>0)
ASVtable_working[,metadata_assign_RH$ID]%>%sum()

file5= read.csv("dada2_out/ASV_Table_LTE_2020_2021_track_final.csv")%>%filter(X%in%gsub("trimmed","trimmed_cutadapt_1.fq",metadata_assign_RH$ID))
sum(file5$Reads.In)

file6
rownames(ASVtable_working)=ASVtable_working$ASV
file6=ASVtable_working[,metadata_assign_RH$ID, drop=F]
rownames(file6)


coverage1= as.data.frame(  apply(file6, 2, function(x) 100*(1-sum(x==1)/sum(x))))
mean(coverage1[,1])
sd(coverage1[,1])
