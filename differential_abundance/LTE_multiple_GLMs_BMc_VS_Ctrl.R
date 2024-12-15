library(ggplot2) #ggplot2: Fancy graphs
library(ggpubr) #ggplot expansion for more fancy graphs
library(readxl) # upload files  
library(reshape2)
library(readr) # read the files
library(dplyr)# Data handling
library(vegan)# Multivariate stat tool
library(tidyr) # Data handling
library(RColorBrewer) # Colours for fancy graphs
library(tibble)# Data handling: rownames/columnanmes
library(stringr) # Manipulate strings
library(randomForest)
set.seed(172590)


#------------------------------------------------------------------------------
# First Part Analyze Bacteria in the rhizosphere
ASV_table= read.csv( file=paste0( "/home/ioannis.kampouris",
  "/DiControl_JKI_Phase3_repos/",
"LTE_2_years/16S/decontam_Rarefied_1000_times_ASV_table_2020_2021.csv"),
row.names = 1)
#Metadata to sample name
metadata=data.frame(Sample=rownames(ASV_table), 
ID=rownames(ASV_table))%>%
  mutate(Sample=str_replace(Sample,".trimmed",""))%>%
  mutate(Sample=str_replace(Sample,".trimmed",""))%>%
  separate(Sample, into = c("Sample", "Year"), sep = "..20" )%>%
  mutate(Time_point= str_replace(Sample, ".*.T0.*.","T0"))%>%
  
  mutate(Soil_Type= str_replace(Sample, ".*.T0.*.","FS_T0"))%>%
mutate(Soil_Type= str_replace(Soil_Type, "FS.*.","FS"))%>%
mutate(Soil_Type= str_replace(Soil_Type, "RH.*","RH"))%>%
mutate(Tillage= str_replace(Sample, ".*CT.*.","CT"))%>%
  mutate(Tillage= str_replace(Tillage, ".*MP.*.","MP"))%>%
  mutate(Int= str_replace(Sample, ".*In.*.","Intensive"))%>%
  mutate(Int= str_replace(Int, ".*Ex.*.","Extensive"))%>%
  mutate(BMc= str_replace(Sample, ".*.t.I.*","BMc"))%>%
  mutate(BMc= str_replace(BMc, ".*.x.I.*","BMc"))%>%
  mutate(BMc= str_replace(BMc, ".*.n.I.*","BMc"))%>%
  mutate(Int= str_replace(Int, ".*Ex.*.","Extensive"))%>%
  mutate(Block=str_sub(Sample, -1, -1)  )%>% 
  mutate(Root_Window=str_replace(Sample, ".*.RW.*.","RW"))
metadata$BMc[metadata$BMc!="BMc"]<-"Ctrl"
metadata$Block=as.numeric(metadata$Block)
metadata$Root_Window[metadata$Root_Window!="RW"]<-"Normal-Sample"

Taxonomy=read.csv(
"/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/LTE_2_years/16S/ASV_taxonomy.csv", row.names = "X")
Taxonomy=mutate(Taxonomy, tax.Order=str_replace_all(tax.Order, " ", "_"))%>%mutate(., tax.Order=str_replace_all(tax.Order, "-", "_"))
Taxonomy=mutate(Taxonomy, tax.Family=str_replace_all(tax.Family, " ", "_"))%>%mutate(., tax.Family=str_replace_all(tax.Family, "-", "_"))%>%
  mutate(.,tax.Family=str_replace_all(tax.Family, "\\[", ""))%>%
  mutate(., tax.Family=str_replace_all(tax.Family, "]_", "_"))
Taxonomy=Taxonomy%>%mutate(., tax.Genus=str_replace_all(tax.Genus, " ", "_"))%>%
  mutate(.,tax.Genus=str_replace_all(tax.Genus, "\\[", ""))%>%
  mutate(., tax.Genus=str_replace_all(tax.Genus, "]_", "_"))%>%mutate(., tax.Genus=str_replace_all(tax.Genus, "-", "_"))

Taxonomy[is.na(Taxonomy)]="Unclassified"
Taxonomy[Taxonomy=="uncultured"]<-"Unclassified"
Taxonomy[Taxonomy=="metagenome"]<-"Unclassified"

for(i in 1:nrow(Taxonomy)){
  gn =paste0(Taxonomy[i,7]) 
  if (Taxonomy[i,7]=="Unclassified"){
    gn =paste0("Unclassified_", Taxonomy[i,6]) 
    if (Taxonomy[i,6]=="Unclassified"){
      gn =paste0("Unclassified_", Taxonomy[i,5]) 
      if (Taxonomy[i,5]=="Unclassified"){
        gn =paste0("Unclassified_", Taxonomy[i,4]) 
        if (Taxonomy[i,4]=="Unclassified"){
          gn =paste0("Unclassified_", Taxonomy[i,3]) 
          
          if (Taxonomy[i,3]=="Unclassified"){
            gn =paste0("Unclassified_", Taxonomy[i,2]) 
            
          }
        }
      }
    }
  }
  Taxonomy[i,8]=paste0(gn)
 print(gn)
}
Taxonomy$Acronym=paste0( Taxonomy$tax.Genus, "_", Taxonomy$ASV)

# Generate count tables over the different taxa.

metadata_RH_2020=filter(metadata, Soil_Type=="RH"&Time_point!="T0"&Root_Window=="Normal-Sample"&Year==20&ID!="RH.MP.Int.C1..2020.trimmed")
metadata_RH_2021=filter(metadata, Soil_Type=="RH"&Time_point!="T0"&Root_Window=="Normal-Sample"&Year==21&ID!="RH.MP.Int.C1..2020.trimmed")
metadata_RH=filter(metadata, Soil_Type=="RH"&Time_point!="T0"&Root_Window=="Normal-Sample"&ID!="RH.MP.Int.C1..2020.trimmed")

ggplot(cbind(metadata_RH, ASV_table[metadata_RH$Sample,]),
        aes(x=ASV9,y=ASV694))+geom_point()+stat_cor(method = "pearson")+
  facet_grid(~BMc~Year)


summed_ASVs= cbind(t(ASV_table),Taxonomy[rownames(t(ASV_table)),
c("ASV","Acronym")])%>%select(-ASV)%>%
  aggregate(.~Acronym, sum)
rownames(summed_ASVs)=summed_ASVs$Acronym
summed_ASVs_check_2020= summed_ASVs[,metadata_RH_2020$ID]%>%t()%>% 
  as.data.frame()%>% rownames_to_column(var="Sample")%>%
  gather(-Sample, key="ASV", value = "RA")%>% filter(RA>1)%>% 
  group_by(ASV)%>% dplyr::summarise(count=n())%>%
  filter(count>=3)

summed_ASVs_check_2021= summed_ASVs[,metadata_RH_2021$ID]%>%t()%>% 
as.data.frame()%>% rownames_to_column(var="Sample")%>%
gather(-Sample, key="ASV", value = "RA")%>% filter(RA>1)%>% group_by(ASV)%>% 
  dplyr::summarise(count=n())%>%
filter(count>=3)

LTE2020_M_central_2021_growth <- 
read_excel("LTE_2_years/LTE2020_M_central_2021_growth.xlsx")%>%as.data.frame()
LTE2020_M_central_2021_growth$Shoot_DW=compositions::clr(LTE2020_M_central_2021_growth$Shoot_DW)
LTE2020_M_central_2021_growth$Name=paste0(LTE2020_M_central_2021_growth$TILLAGE,
                                          LTE2020_M_central_2021_growth$FERTILIZATION,LTE2020_M_central_2021_growth$BENEFICIALS,
                                          LTE2020_M_central_2021_growth$BLOCK, LTE2020_M_central_2021_growth$YEAR)
rownames( LTE2020_M_central_2021_growth)=LTE2020_M_central_2021_growth$Name

metadata_RH_2020$Name=paste0(metadata_RH_2020$Tillage,metadata_RH_2020$Int,metadata_RH_2020$BMc,metadata_RH_2020$Block, paste0("20", metadata_RH_2020$Year))
metadata_RH_2021$Name=paste0(metadata_RH_2021$Tillage,metadata_RH_2021$Int,metadata_RH_2021$BMc,metadata_RH_2021$Block, paste0("20", metadata_RH_2021$Year))

multiple_glm= function(x){
   
  file1$Treatment2=file1$BMc
  file1$Treatment2[file1$Treatment2=="BMc"]<-1
  file1$Treatment2[file1$Treatment2=="Ctrl"]<-0
  file1$Treatment2=as.numeric(  file1$Treatment2)


  Table_glm=glm(file1$Treatment2~scale(x), family=binomial(link="logit"))
  file2=anova(Table_glm, test = "Chisq")
  FC1=dcast(data.frame(T=file1$BMc, A=x), T~., value.var = "A",  mean)
  FC2=log2((FC1[1,2]+1)/(FC1[2,2]+1))
  file3=data.frame(coeff=(Table_glm$coefficients[2]), p=file2$`Pr(>Chi)`[2], FC=FC2, Mean=mean(x))
  
  return(file3)
}


comparison_list=list()
for (i in unique(paste0( metadata_RH$Tillage,"_", metadata_RH$Int,"_", metadata_RH$Year ))){
  print(i)
  file1=mutate(metadata_RH, Group=paste0( metadata_RH$Tillage,"_", metadata_RH$Int, "_", metadata_RH$Year))
  file1=filter(file1, Group==paste0(i) )
  input1=(as.data.frame(ASV_table[file1$ID,]))
  input2=as.data.frame( apply(input1, 2, function (x) sum(x>0)))
  colnames(input2)<-"Count"
  input2=filter(input2, Count>3)
  input3=as.data.frame(( input1[file1$ID,rownames(input2)]))
  multiple_glms_2020_ASVs= apply(input3,2, function(x) multiple_glm(x))
  multiple_glms_ASVs_df <- as.data.frame(do.call(rbind, multiple_glms_2020_ASVs))
  multiple_glms_ASVs_df[is.na(multiple_glms_ASVs_df)]<-1
  multiple_glms_ASVs_df$padj=p.adjust(multiple_glms_ASVs_df$p, method="BH")
  multiple_glms_ASVs_df=filter(multiple_glms_ASVs_df)%>%mutate(Int=unique(file1$Int),
  Tillage=unique(file1$Tillage), Year=unique(file1$Year))
  comparison_list=rbind(multiple_glms_ASVs_df%>%rownames_to_column(var="ASV"),comparison_list)
  
}
Taxonomy2=filter(Taxonomy, ASV%in%c(unique(comparison_list$ASV)))

write_csv(as.data.frame(comparison_list)%>%full_join(Taxonomy2),
file=paste0( "/home/ioannis.kampouris",
"/DiControl_JKI_Phase3_repos/",
"LTE_2_years/16S/RH_ASV_Responders"))

comparison_list2=filter(as.data.frame(comparison_list)%>%full_join(Taxonomy2), padj<0.05)
library(ggplot2)
bacterial_ASVs= ggplot(comparison_list2, aes(x=(100*Mean/20000), y=FC, 
fill=paste0(tax.Phylum), shape=paste0(Tillage,"-", Int)))+
  facet_wrap(~paste0("20", Year), scales = "fixed")+
  geom_point( size=6, colour="black")+scale_shape_manual(values = c(21,22,23,24), name="")+
  ylab("Log2FC(BMc/Ctrl)")+xlab("Mean Relative Abundance (%)") + 
   scale_fill_brewer(palette = "Paired",  name="Phylum")  +
  theme_bw() +
  
  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=45, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
  theme(plot.title =  element_text( size=45, face="bold", colour = "black"))+
  
  theme( strip.text = element_text( size=25, face = "bold", colour="white")) +
  theme( strip.background = element_rect( fill="skyblue4")) +
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+ 
  guides(fill = guide_legend(override.aes = list(shape = 21, size=10)))+ 
  guides(colour = guide_legend(override.aes = list(linewidth = 0)))+guides(shape = guide_legend(override.aes = list(size=14)))+
  ggtitle("Differential abundance of bacterial ASVs")
  
bacterial_ASVs

#Select Comamonadeceae sequences of differentially abundant ASVs for further analysis
sequences=unique(paste0(d$Sequence,"_",d$ASV))
comamonadaceae= Biostrings::DNAStringSet(c(gsub("_.*","", sequences) ))
names(comamonadaceae)=gsub("*.*._","", sequences)
Biostrings::writeXStringSet(comamonadaceae,
paste0(  "/home/ioannis.kampouris",
"/DiControl_JKI_Phase3_repos/",
"LTE_2_years/16S/comamonadaceae.fasta")  , append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

blast_results=read.table(paste0(
  "/home/ioannis.kampouris",
  "/DiControl_JKI_Phase3_repos/",
  "LTE_2_years/16S/output_blast_results_comamonadaceae.txt"
))

# Second part: Moving to Fungi
ASV_table= read.delim( file=paste0( "/home/ioannis.kampouris",
"/DiControl_JKI_Phase3_repos/",
"LTE_2_years/18S/fungal_ASV_Rarefied_Table_1000_times.txt"), sep = "\t", row.names = 1)%>%t()

Taxonomy= read.delim( file=paste0( "/home/ioannis.kampouris",
"/DiControl_JKI_Phase3_repos/",
"LTE_2_years/18S/ITS2_20202021_TAX_reduced.txt"), sep = "\t")


Taxonomy$ASVname=paste0("ASV", 1:nrow(Taxonomy))
metadata=data.frame(Sample=colnames(ASV_table), ID=colnames(ASV_table))%>%
  mutate(Year=str_replace(Sample,".*.20.*.","20"))%>%
  mutate(Year=str_replace(Year,".*.21*.*.","21"))%>%
  mutate(Time_point= str_replace(Sample, ".*.T0.*.","T0"))%>%
  
  mutate(Soil_Type= str_replace(Sample, ".*.T0.*.","FS_T0"))%>%
  
  mutate(Soil_Type= str_replace(Soil_Type, "FS.*.","FS"))%>%
  mutate(Soil_Type= str_replace(Soil_Type, "RH.*","RH"))%>%
  mutate(Tillage= str_replace(Sample, ".*CT.*.","CT"))%>%
  mutate(Tillage= str_replace(Tillage, ".*MP.*.","MP"))%>%
  mutate(Int= str_replace(Sample, ".*In.*.","Intensive"))%>%
  mutate(Int= str_replace(Int, ".*Ex.*.","Extensive"))%>%
  mutate(BMc= str_replace(Sample, ".*BM*.*","BMc"))%>%
  mutate(Int= str_replace(Int, ".*Ex.*.","Extensive"))%>%
  mutate(Block=str_sub(Sample, -1, -1)  )%>% 
  mutate(Root_Window=str_replace(Sample, ".*.RW.*.","RW"))
metadata$BMc[metadata$BMc!="BMc"]<-"Ctrl"

metadata$Block=as.numeric(metadata$Block)
metadata$Root_Window[metadata$Root_Window!="RW"]<-"Normal-Sample"



Taxonomy[is.na(Taxonomy)]="Unclassified"
Taxonomy[Taxonomy=="uncultured"]<-"Unclassified"
Taxonomy[Taxonomy=="metagenome"]<-"Unclassified"
Taxonomy$ASV=rownames(Taxonomy)
Taxonomy=select(Taxonomy, ASV, everything())
for(i in 1:nrow(Taxonomy)){
  gn =paste0(Taxonomy[i,7]) 
  if (Taxonomy[i,7]=="Unclassified"){
    gn =paste0("Unclassified_", Taxonomy[i,6]) 
    if (Taxonomy[i,6]=="Unclassified"){
      gn =paste0("Unclassified_", Taxonomy[i,5]) 
      if (Taxonomy[i,5]=="Unclassified"){
        gn =paste0("Unclassified_", Taxonomy[i,4]) 
        if (Taxonomy[i,4]=="Unclassified"){
          gn =paste0("Unclassified_", Taxonomy[i,3]) 
          
          if (Taxonomy[i,3]=="Unclassified"){
            gn =paste0("Unclassified_", Taxonomy[i,2]) 
            
          }
        }
      }
    }
  }
  Taxonomy[i,8]=paste0(gn)
  print(gn)
}
Taxonomy$ASV2=paste0(Taxonomy$ASV, rownames(Taxonomy))
Taxonomy$Acronym=paste0( Taxonomy$genus, "_", Taxonomy$ASV2)


# Generate count tables over the different taxa.

metadata_RH_2020=filter(metadata, Soil_Type=="RH"&Root_Window=="Normal-Sample"&Year==20&ID!="RH.MP.Int.C1..2020.trimmed")
metadata_RH_2021=filter(metadata, Soil_Type=="RH"&Root_Window=="Normal-Sample"&Year==21&ID!="RH.MP.Int.C1..2020.trimmed")
metadata_RH=filter(metadata, Soil_Type=="RH"&Root_Window=="Normal-Sample"&ID!="RH.MP.Int.C1..2020.trimmed")


FS_summed_ASVs= cbind((ASV_table),Taxonomy[rownames((ASV_table)),c("ASV","Acronym")])%>%select(-ASV)%>%
  aggregate(.~Acronym, sum)
rownames(FS_summed_ASVs)=FS_summed_ASVs$Acronym
FS_summed_ASVs_check_2020= FS_summed_ASVs[,metadata_RH_2020$ID]%>%t()%>% 
  as.data.frame()%>% rownames_to_column(var="Sample")%>%
  gather(-Sample, key="ASV", value = "RA")%>% filter(RA>1)%>% 
  group_by(ASV)%>% dplyr::summarise(count=n())%>%
  filter(count>=3)


FS_summed_ASVs_check_2021= FS_summed_ASVs[,metadata_RH_2021$ID]%>%t()%>% 
  as.data.frame()%>% rownames_to_column(var="Sample")%>%
  gather(-Sample, key="ASV", value = "RA")%>% filter(RA>1)%>% group_by(ASV)%>% dplyr::summarise(count=n())%>%
  filter(count>=3)


multiple_glm= function(x){
  
  file1$Treatment2=file1$BMc
  file1$Treatment2[file1$Treatment2=="BMc"]<-1
  file1$Treatment2[file1$Treatment2=="Ctrl"]<-0
  file1$Treatment2=as.numeric(  file1$Treatment2)
  
  
  Table_glm=glm(file1$Treatment2~scale(x), family=binomial(link="logit"))
  file2=anova(Table_glm, test = "Chisq")
  FC1=dcast(data.frame(T=file1$BMc, A=x), T~., value.var = "A",  mean)
  FC2=log2((FC1[1,2]+1)/(FC1[2,2]+1))
  file3=data.frame(coeff=(Table_glm$coefficients[2]), p=file2$`Pr(>Chi)`[2], FC=FC2, Mean=mean(x))
  
  return(file3)
}


comparison_list=list()
for (i in unique(paste0( metadata_RH$Tillage,"_", metadata_RH$Int,"_", metadata_RH$Year ))){
  print(i)
  file1=mutate(metadata_RH, Group=paste0( metadata_RH$Tillage,"_", metadata_RH$Int, "_", metadata_RH$Year))
  file1=filter(file1, Group==paste0(i) )
  input1=(as.data.frame(t(ASV_table)[file1$ID,]))
  input2=as.data.frame( apply(input1, 2, function (x) sum(x>0)))
  colnames(input2)<-"Count"
  input2=filter(input2, Count>3)
  input3=as.data.frame(( input1[file1$ID,rownames(input2)]))
  multiple_glms_2020_ASVs= apply(input3,2, function(x) multiple_glm(x))
  multiple_glms_ASVs_df <- as.data.frame(do.call(rbind, multiple_glms_2020_ASVs))
  multiple_glms_ASVs_df[is.na(multiple_glms_ASVs_df)]<-1
  multiple_glms_ASVs_df$padj=p.adjust(multiple_glms_ASVs_df$p, method="BH")
  multiple_glms_ASVs_df=filter(multiple_glms_ASVs_df)%>%mutate(Int=unique(file1$Int),
                                                               Tillage=unique(file1$Tillage), Year=unique(file1$Year))
  comparison_list=rbind(multiple_glms_ASVs_df%>%rownames_to_column(var="ASV"),comparison_list)
  
}

Taxonomy2=filter(Taxonomy, ASV%in%c(unique(comparison_list$ASV)))


write_csv(as.data.frame(comparison_list)%>%full_join(Taxonomy2),file=paste0( "/home/ioannis.kampouris",
"/DiControl_JKI_Phase3_repos/",
"LTE_2_years/18S/Fungal_ASV_Responders"))


comparison_list2=filter(as.data.frame(comparison_list)%>%full_join(Taxonomy2), padj<0.05)
library(ggplot2)

comparison_list2$phylum=gsub("Fungi_phy_Incertae_sedis","Incertae_sedis",
                             comparison_list2$phylum)
fungal_ASVs=
ggplot(comparison_list2, aes(x=(100*Mean/48690), y=FC, fill=paste0(phylum), shape=paste0(Tillage,"-", Int)))+
  facet_wrap(~paste0("20", Year), scales = "fixed")+
  geom_point( size=6,  colour="black")+scale_shape_manual(values = c(21,22,23,24), name="")+
  ylab("Log2FC(BMc/Ctrl)")+xlab("Mean Relative Abundance (%)") + 
  scale_fill_brewer(palette = "Paired", name="Phylum")  +
  theme_bw() +
  
  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=45, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
  theme(plot.title =  element_text( size=45, face="bold", colour = "black"))+
  
  theme( strip.text = element_text( size=25, face = "bold", colour="white")) +
  theme( strip.background = element_rect( fill="skyblue4")) +
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
                                                                   size = 2, linetype = "solid"))+ 
  guides(fill = guide_legend(override.aes = list(shape = 21, size=10)))+ 
  guides(colour = guide_legend(override.aes = list(linewidth = 0)))+
  guides(shape = guide_legend(override.aes = list(size=14))) +
scale_fill_brewer(palette = "Paired", name="")+
  ggtitle("Differential abundance of fungal ASVs")


combo1= ggarrange(bacterial_ASVs, 
          fungal_ASVs, 
          font.label = list(size = 45, color = "black", face = "bold", family = NULL), ncol=1,
          labels=c("A)","B)"),align = "hv")

combo2= ggarrange(classes_plot_bacteria, 
                  classes_plot_fungi, 
                  font.label = list(size = 45, color = "black", face = "bold", family = NULL), ncol=1,
                  labels=c("C)","D)"),align = "hv")

png(filename = "figure4MS4.png", width = 14000, height = 8000,
    res = 300)
print(ggarrange(combo1,
combo2))
dev.off()
nrow(filter(comparison_list2, FC>0&Year==20))
nrow(filter(comparison_list2, FC<0&Year==20))


nrow(filter(comparison_list2, FC>0&Year==21))
nrow(filter(comparison_list2, FC<0&Year==21))


ASV_table= read.csv( file=paste0( "/home/ioannis.kampouris",
                                  "/DiControl_JKI_Phase3_repos/",
                                  "LTE_2_years/16S/decontam_Rarefied_1000_times_ASV_table_2020_2021.csv"), row.names = 1)



metadata=data.frame(Sample=rownames(ASV_table), ID=rownames(ASV_table))%>%
  mutate(Sample=str_replace(Sample,".trimmed",""))%>%
  mutate(Sample=str_replace(Sample,".trimmed",""))%>%
  separate(Sample, into = c("Sample", "Year"), sep = "..20" )%>%
  mutate(Time_point= str_replace(Sample, ".*.T0.*.","T0"))%>%
  
  mutate(Soil_Type= str_replace(Sample, ".*.T0.*.","FS_T0"))%>%
  
  mutate(Soil_Type= str_replace(Soil_Type, "FS.*.","FS"))%>%
  mutate(Soil_Type= str_replace(Soil_Type, "RH.*","RH"))%>%
  mutate(Tillage= str_replace(Sample, ".*CT.*.","CT"))%>%
  mutate(Tillage= str_replace(Tillage, ".*MP.*.","MP"))%>%
  mutate(Int= str_replace(Sample, ".*In.*.","Intensive"))%>%
  mutate(Int= str_replace(Int, ".*Ex.*.","Extensive"))%>%
  mutate(BMc= str_replace(Sample, ".*.t.I.*","BMc"))%>%
  mutate(BMc= str_replace(BMc, ".*.x.I.*","BMc"))%>%
  mutate(BMc= str_replace(BMc, ".*.n.I.*","BMc"))%>%
  mutate(Int= str_replace(Int, ".*Ex.*.","Extensive"))%>%
  mutate(Block=str_sub(Sample, -1, -1)  )%>% 
  mutate(Root_Window=str_replace(Sample, ".*.RW.*.","RW"))
metadata$BMc[metadata$BMc!="BMc"]<-"Ctrl"
metadata$Block=as.numeric(metadata$Block)
metadata$Root_Window[metadata$Root_Window!="RW"]<-"Normal-Sample"

Taxonomy=read.csv("/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/LTE_2_years/16S/ASV_taxonomy.csv", row.names = "X")
Taxonomy=mutate(Taxonomy, tax.Order=str_replace_all(tax.Order, " ", "_"))%>%mutate(., tax.Order=str_replace_all(tax.Order, "-", "_"))
Taxonomy=mutate(Taxonomy, tax.Family=str_replace_all(tax.Family, " ", "_"))%>%mutate(., tax.Family=str_replace_all(tax.Family, "-", "_"))%>%
  mutate(.,tax.Family=str_replace_all(tax.Family, "\\[", ""))%>%
  mutate(., tax.Family=str_replace_all(tax.Family, "]_", "_"))
Taxonomy=Taxonomy%>%mutate(., tax.Genus=str_replace_all(tax.Genus, " ", "_"))%>%
  mutate(.,tax.Genus=str_replace_all(tax.Genus, "\\[", ""))%>%
  mutate(., tax.Genus=str_replace_all(tax.Genus, "]_", "_"))%>%mutate(., tax.Genus=str_replace_all(tax.Genus, "-", "_"))


Taxonomy[is.na(Taxonomy)]="Unclassified"
Taxonomy[Taxonomy=="uncultured"]<-"Unclassified"
Taxonomy[Taxonomy=="metagenome"]<-"Unclassified"

for(i in 1:nrow(Taxonomy)){
  gn =paste0(Taxonomy[i,7]) 
  if (Taxonomy[i,7]=="Unclassified"){
    gn =paste0("Unclassified_", Taxonomy[i,6]) 
    if (Taxonomy[i,6]=="Unclassified"){
      gn =paste0("Unclassified_", Taxonomy[i,5]) 
      if (Taxonomy[i,5]=="Unclassified"){
        gn =paste0("Unclassified_", Taxonomy[i,4]) 
        if (Taxonomy[i,4]=="Unclassified"){
          gn =paste0("Unclassified_", Taxonomy[i,3]) 
          
          if (Taxonomy[i,3]=="Unclassified"){
            gn =paste0("Unclassified_", Taxonomy[i,2]) 
            
          }
        }
      }
    }
  }
  Taxonomy[i,8]=paste0(gn)
  print(gn)
}
Taxonomy$Acronym=paste0( Taxonomy$tax.Genus, "_", Taxonomy$ASV)

