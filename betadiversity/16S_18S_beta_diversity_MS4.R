library(ggplot2) #ggplot2: Fancy graphs
library(ggpubr) #ggplot expansion for more fancy graphs
library(readxl) # upload files  
library(readr) # read the files
library(dplyr)# Data handling
library(vegan)# Multivariate stat tool
library(tidyr) # Data handling
library(RColorBrewer) # Colours for fancy graphs
library(tibble)# Data handling: rownames/columnanmes
library(stringr) # Manipulate strings

#Load the 16S data 
ASV_table= read.csv( file=paste0( "/home/ioannis.kampouris",
"/DiControl_JKI_Phase3_repos/",
"LTE_2_years/16S/decontam_Rarefied_1000_times_ASV_table_2020_2021.csv"),  
sep = ",", row.names = 1)
#Load the metadata 
#Use the samples names for creating the metadata
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
  mutate(Block=str_sub(Sample, -2, -2)  )%>% 
  mutate(Root_Window=str_replace(Sample, ".*.RW.*.","RW"))
metadata$BMc[metadata$BMc!="BMc"]<-"Ctrl"
metadata$Root_Window[metadata$Root_Window!="RW"]<-"Normal-Sample"


#PERMANOVA and ANOSIM for Rhizospheric Soil with Bray-Curtis Distance
metadata_RH=filter(metadata, Soil_Type=="RH"&Root_Window=="Normal-Sample"&ID!="RH.MP.Int.C1..2020.trimmed")
metadata_RH_2020=filter(metadata, Soil_Type=="RH"&Root_Window=="Normal-Sample"&Year==20&ID!="RH.MP.Int.C1..2020.trimmed")
metadata_RH_2021=filter(metadata, Soil_Type=="RH"&Root_Window=="Normal-Sample"&Year==21&ID!="RH.MP.Int.C1..2020.trimmed")

#Calculate the distances
bray_RH_distance_20 = 
  vegdist((ASV_table[metadata_RH_2020$ID,]), method = "bray")

bray_RH_distance_21 = 
  vegdist((ASV_table[metadata_RH_2021$ID,]), method = "bray")

bray_permanova_RH_20= adonis2(bray_RH_distance_20~Tillage*Int*BMc,
                              data =metadata_RH_2020, 
                              permutations =9999 )
bray_permanova_RH_20

bray_permanova_RH_21= adonis2(bray_RH_distance_21~Tillage*Int*BMc, 
                              data =metadata_RH_2021, 
                              permutations =9999 )
bray_permanova_RH_21
permutest(betadisper(bray_RH_distance_20, metadata_RH_2020$BMc),
          permutations = 9999)
permutest(betadisper(bray_RH_distance_20, metadata_RH_2020$Tillage),
          permutations = 9999)
permutest(betadisper(bray_RH_distance_20, metadata_RH_2020$Int),
          permutations = 9999)
permutest(betadisper(bray_RH_distance_21, metadata_RH_2021$BMc),
          permutations = 9999)
permutest(betadisper(bray_RH_distance_21, metadata_RH_2021$Tillage),
          permutations = 9999)
permutest(betadisper(bray_RH_distance_21, metadata_RH_2021$Int),
          permutations = 9999)


#MDS for Rhizosphere with Bray-Curtis Distance
bray_meta_MDS_soil_RH_20= cmdscale(bray_RH_distance_20,eig = T)
bray_meta_MDS_soil_RH_20.scrs <- as.data.frame(scores(bray_meta_MDS_soil_RH_20,display = "sites"))%>% 
  rownames_to_column(var="ID") %>% full_join(., metadata_RH_2020, by="ID")


bray_meta_MDS_soil_RH_21= cmdscale(bray_RH_distance_21,eig = T)
bray_meta_MDS_soil_RH_21.scrs <- as.data.frame(scores(bray_meta_MDS_soil_RH_21,
                                                      display = "sites"))%>% 
  rownames_to_column(var="ID") %>% full_join(., metadata_RH_2021, by="ID")


bray_meta_MDS_soil_RH_21.scrs$tillage_fertilization=paste0(bray_meta_MDS_soil_RH_21.scrs$Tillage, "-",bray_meta_MDS_soil_RH_21.scrs$Int)

bray_meta_MDS_soil_RH_21.scrs$inoculation <- factor(bray_meta_MDS_soil_RH_21.scrs$BMc, levels = c("Ctrl", "BMc"))
bray_meta_MDS_soil_RH_21.scrs$tillage_fertilization_inoculation=paste0(bray_meta_MDS_soil_RH_21.scrs$Tillage, "-",
bray_meta_MDS_soil_RH_21.scrs$Int, "-",bray_meta_MDS_soil_RH_21.scrs$BMc)

bray_meta_MDS_soil_RH_21.scrs$tillage_fertilization_inoculation <- factor(bray_meta_MDS_soil_RH_21.scrs$tillage_fertilization_inoculation, 
                                                 levels = c("MP-Intensive-Ctrl", "MP-Intensive-BMc", "MP-Extensive-Ctrl", "MP-Extensive-BMc",
                                                            "CT-Intensive-Ctrl", "CT-Intensive-BMc", "CT-Extensive-Ctrl", "CT-Extensive-BMc"))
mdsplot <- ggplot(bray_meta_MDS_soil_RH_21.scrs, aes(x = Dim1, y = Dim2, fill = tillage_fertilization, linetype=Tillage, shape = inoculation)) +
  geom_point(size =6)+
#  stat_ellipse(aes(group = paste0(inoculation, Tillage)), level = 0.95, show.legend = F, size = 0.5) +
  theme(legend.text = element_text(colour = "black", size = 22, face = "bold"), 
        legend.title = element_text(colour = "black", size = 20, face = "bold"),
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line("white"),
        legend.box.background = element_rect(colour = "white"), 
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.line.x.top = element_blank(),
        axis.line.x.bottom = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        axis.line.y.right = element_blank(),
        axis.text = element_text(colour = "black", size = 20))

mdsplot <- mdsplot + scale_shape_manual(values = c(21,24),
                                        name = "",
                                        breaks = c("Ctrl", "BMc"),
                                        labels = c("Ctrl", "BMc"))

mdsplot <- mdsplot + scale_fill_manual(values = c("#8CC0DE", "#B8F1B0","#F38181", "#FCE38A"),
                                       name = "",
                                       breaks = c("MP-Intensive", "MP-Extensive", "CT-Intensive", "CT-Extensive"),
                                       labels = c("MP-Int", "MP-Ext","CT-Int", "CT-Ext"))

mdsplot <- mdsplot + scale_colour_manual(values = c("black","black","black", "black","black","black","black", "black"))
mdsplot <- mdsplot #+ scale_x_continuous(labels = unicode_minus) + scale_y_continuous(labels = unicode_minus)
mdsplot_2021 <- mdsplot + guides(fill = guide_legend(override.aes = list(shape = 21)))
mdsplot_2021=mdsplot_2021+ xlab(paste("MDS1","[",round(100*bray_meta_MDS_soil_RH_21$eig[1]/sum(bray_meta_MDS_soil_RH_21$eig),2),"%]"))+
  ylab(paste("MDS2","[",round(100*bray_meta_MDS_soil_RH_21$eig[2]/sum(bray_meta_MDS_soil_RH_21$eig),2),"%]"))+
  theme(axis.title = element_text(size=30))+theme(legend.position = "right")


png(filename = "Loreen_MDS_2021.png", width=550, height=450)
print(mdsplot_2021)
dev.off()


bray_meta_MDS_soil_RH_20.scrs$tillage_fertilization=paste0(bray_meta_MDS_soil_RH_20.scrs$Tillage, "-",bray_meta_MDS_soil_RH_20.scrs$Int)

bray_meta_MDS_soil_RH_20.scrs$inoculation <- factor(bray_meta_MDS_soil_RH_20.scrs$BMc, levels = c("Ctrl", "BMc"))
bray_meta_MDS_soil_RH_20.scrs$tillage_fertilization_inoculation=paste0(bray_meta_MDS_soil_RH_20.scrs$Tillage, "-",
                                                                       bray_meta_MDS_soil_RH_20.scrs$Int, "-",bray_meta_MDS_soil_RH_20.scrs$BMc)

bray_meta_MDS_soil_RH_20.scrs$tillage_fertilization_inoculation <- factor(bray_meta_MDS_soil_RH_20.scrs$tillage_fertilization_inoculation, 
                                                                          levels = c("MP-Intensive-Ctrl", "MP-Intensive-BMc", "MP-Extensive-Ctrl", "MP-Extensive-BMc",
                                                                                     "CT-Intensive-Ctrl", "CT-Intensive-BMc", "CT-Extensive-Ctrl", "CT-Extensive-BMc"))
mdsplot <- ggplot(bray_meta_MDS_soil_RH_20.scrs, aes(x = Dim1, y = Dim2, 
                                                     shape = tillage_fertilization, linetype=Tillage, fill = inoculation)) +
  geom_point(size =6)+
#  stat_ellipse(aes(group = paste0(inoculation, Tillage)), level = 0.95, show.legend = F, size = 0.5) +
  theme(legend.text = element_text(colour = "black", size = 22, face = "bold"), 
        legend.title = element_text(colour = "black", size = 20, face = "bold"),
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line("white"),
        legend.box.background = element_rect(colour = "white"), 
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.line.x.top = element_blank(),
        axis.line.x.bottom = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        axis.line.y.right = element_blank(),
        axis.text = element_text(colour = "black", size = 20))+
 scale_shape_manual(values = c(21,22,23,24),
                       name = "",
                       breaks = c("MP-Intensive", "MP-Extensive", "CT-Intensive", "CT-Extensive"),
                       labels = c("MP-Int", "MP-Ext","CT-Int", "CT-Ext")
  )+ scale_fill_manual(values = c( "#146eb4","#ff9900"),name = "",
                       breaks = c("Ctrl", "BMc"),
                       labels = c("Ctrl", "BMc")
  )

mdsplot <- mdsplot + scale_colour_manual(values = c("black","black","black", "black","black","black","black", "black"))
mdsplot <- mdsplot #+ scale_x_continuous(labels = unicode_minus) + scale_y_continuous(labels = unicode_minus)
mdsplot_2020 <- mdsplot + guides(fill = guide_legend(override.aes = list(shape = 21)))

mdsplot_2020=mdsplot_2020 + xlab(paste("MDS1","[",round(100*bray_meta_MDS_soil_RH_20$eig[1]/sum(bray_meta_MDS_soil_RH_20$eig),2),"%]"))+
   ylab(paste("MDS2","[",round(100*bray_meta_MDS_soil_RH_20$eig[2]/sum(bray_meta_MDS_soil_RH_20$eig),2),"%]"))+
  theme(axis.title = element_text(size=30))+theme(legend.position = "right")
png(filename = "Loreen_MDS_2020.png", width=550, height=450)
print(mdsplot_2020)
dev.off()



#``````````````````````````````````````````````````````````````````````````
#Fungi
#``````````````````````````````````````````````````````````````````````````````
#Fungal Table
f_ASV_table= read.delim( file=paste0( "/home/ioannis.kampouris",
  "/DiControl_JKI_Phase3_repos/",
  "LTE_2_years/18S/fungal_ASV_Rarefied_Table_1000_times.txt"), sep = "\t", row.names = 1)
f_metadata=data.frame(ID=rownames(f_ASV_table), Sample=rownames(f_ASV_table))%>%
  separate(.,Sample,into = c("Soil_Type","Year","Int","Tillage"), sep="_")
f_metadata$BMc=gsub("BM.","BMc",gsub("C.","Ctrl", gsub("*.*_","", f_metadata$ID)))
f_metadata$Tillage=paste0(f_metadata$Tillage,"ensive")

#PERMANOVA and ANOSIM for Rhizospheric Soil with Bray-Curtis Distance
f_metadata=filter(f_metadata, Soil_Type=="RH")
f_metadata_RH_2020=filter(f_metadata, Soil_Type=="RH"&Year==20)
f_metadata_RH_2021=filter(f_metadata, Soil_Type=="RH"&Year==21)
#Fungal bray distance
f_bray_RH_distance_20 = 
  vegdist((f_ASV_table[f_metadata_RH_2020$ID,]), method = "bray")

f_bray_RH_distance_21 = 
  vegdist((f_ASV_table[f_metadata_RH_2021$ID,]), method = "bray")

f_bray_permanova_RH_20= adonis2(f_bray_RH_distance_20~Tillage*Int*BMc,
                              data =f_metadata_RH_2020, 
                              permutations =9999 )
f_bray_permanova_RH_20

f_bray_permanova_RH_21= adonis2(f_bray_RH_distance_21~Tillage*Int*BMc, 
                              data =f_metadata_RH_2021, 
                              permutations =9999 )
f_bray_permanova_RH_21

permutest(betadisper(f_bray_RH_distance_20, f_metadata_RH_2020$BMc),
          permutations = 9999)
permutest(betadisper(f_bray_RH_distance_20, f_metadata_RH_2020$Tillage),
          permutations = 9999)
permutest(betadisper(f_bray_RH_distance_20, f_metadata_RH_2020$Int),
          permutations = 9999)
permutest(betadisper(f_bray_RH_distance_21, f_metadata_RH_2021$BMc),
          permutations = 9999)
permutest(betadisper(f_bray_RH_distance_21, f_metadata_RH_2021$Tillage),
          permutations = 9999)
permutest(betadisper(f_bray_RH_distance_21, f_metadata_RH_2021$Int),
          permutations = 9999)


#NMDS for Rhizosphere with Bray-Curtis Distance


f_bray_meta_MDS_soil_RH_21= cmdscale(f_bray_RH_distance_21,k=2, eig = T)
f_bray_meta_MDS_soil_RH_21.scrs <- as.data.frame(scores(f_bray_meta_MDS_soil_RH_21,
                                                      display = "sites"))%>% 
  rownames_to_column(var="ID") %>% full_join(., f_metadata_RH_2021, by="ID")


f_bray_meta_MDS_soil_RH_21.scrs$tillage_fertilization=paste0(f_bray_meta_MDS_soil_RH_21.scrs$Int, "-",f_bray_meta_MDS_soil_RH_21.scrs$Tillage)

f_bray_meta_MDS_soil_RH_21.scrs$inoculation <- factor(f_bray_meta_MDS_soil_RH_21.scrs$BMc, levels = c("Ctrl", "BMc"))
f_bray_meta_MDS_soil_RH_21.scrs$tillage_fertilization_inoculation=paste0(f_bray_meta_MDS_soil_RH_21.scrs$Tillage, "-",
                                                                       f_bray_meta_MDS_soil_RH_21.scrs$Int, "-",f_bray_meta_MDS_soil_RH_21.scrs$BMc)

f_mdsplot <- ggplot(f_bray_meta_MDS_soil_RH_21.scrs, 
aes(x = Dim1, y = Dim2, shape = tillage_fertilization, linetype=inoculation, fill = inoculation)) +
  geom_point(size =6)+
  # stat_ellipse(aes(group = paste0(inoculation, Tillage)), level = 0.95, show.legend = F, size = 0.5) +
  theme(legend.text = element_text(colour = "black", size = 12, face = "bold"), 
        legend.title = element_text(colour = "black", size = 10, face = "bold"),
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line("white"),
        legend.box.background = element_rect(colour = "white"), 
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.line.x.top = element_blank(),
        axis.line.x.bottom = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        axis.line.y.right = element_blank(),
        axis.text = element_text(colour = "black", size = 10))+ scale_shape_manual(values = c(21,22,23,24),
                                                                                   name = "")+ scale_fill_manual(values = c( "#146eb4","#ff9900"),name = "",
                             breaks = c("Ctrl", "BMc"),
                             labels = c("Ctrl", "BMc")
        ) +
  scale_colour_manual(values = c("black","black","black", "black","black","black","black", "black"))+ guides(fill = guide_legend(override.aes = list(shape = 21)))+ xlab(paste("MDS1","[",round(100*f_bray_meta_MDS_soil_RH_21$eig[1]/sum(f_bray_meta_MDS_soil_RH_21$eig),2),"%]"))+
  ylab(paste("MDS2","[",round(100*f_bray_meta_MDS_soil_RH_21$eig[2]/sum(f_bray_meta_MDS_soil_RH_21$eig),2),"%]"))
                    
                    
f_mdsplot_2021 <- f_mdsplot + guides(fill = guide_legend(override.aes = list(shape = 21)))

png(filename = "f_Loreen_MDS_2021.png", width=550, height=450)
print(f_mdsplot_2021)
dev.off()

f_bray_meta_MDS_soil_RH_20=  cmdscale(f_bray_RH_distance_20,k=2,eig=T)
f_bray_meta_MDS_soil_RH_20.scrs <- as.data.frame(scores(f_bray_meta_MDS_soil_RH_20,display = "sites"))%>% 
  rownames_to_column(var="ID") %>% full_join(., f_metadata_RH_2020, by="ID")

f_bray_meta_MDS_soil_RH_20.scrs$tillage_fertilization=paste0(f_bray_meta_MDS_soil_RH_20.scrs$Int,"-",f_bray_meta_MDS_soil_RH_20.scrs$Tillage)

f_bray_meta_MDS_soil_RH_20.scrs$inoculation <- factor(f_bray_meta_MDS_soil_RH_20.scrs$BMc, levels = c("Ctrl", "BMc"))
f_bray_meta_MDS_soil_RH_20.scrs$tillage_fertilization_inoculation=paste0(f_bray_meta_MDS_soil_RH_20.scrs$Tillage, "-",
                                                                       f_bray_meta_MDS_soil_RH_20.scrs$Int, "-",f_bray_meta_MDS_soil_RH_20.scrs$BMc)

f_mdsplot <- ggplot(f_bray_meta_MDS_soil_RH_20.scrs, 
                    aes(x = Dim1, y = Dim2, shape = tillage_fertilization, linetype=inoculation, fill = inoculation)) +
  geom_point(size =6)+
  # stat_ellipse(aes(group = paste0(inoculation, Tillage)), level = 0.95, show.legend = F, size = 0.5) +
  theme(legend.text = element_text(colour = "black", size = 12, face = "bold"), 
        legend.title = element_text(colour = "black", size = 10, face = "bold"),
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line("white"),
        legend.box.background = element_rect(colour = "white"), 
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.line.x.top = element_blank(),
        axis.line.x.bottom = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        axis.line.y.right = element_blank(),
        axis.text = element_text(colour = "black", size = 10))+ scale_shape_manual(values = c(21,22,23,24), name="")+ scale_fill_manual(values = c( "#146eb4","#ff9900"),name = "",
                             breaks = c("Ctrl", "BMc"),
                             labels = c("Ctrl", "BMc")
        ) +
  scale_colour_manual(values = c("black","black","black", "black","black","black","black", "black"))+ guides(fill = guide_legend(override.aes = list(shape = 21)))+ 
  xlab(paste("MDS1","[",round(100*f_bray_meta_MDS_soil_RH_20$eig[1]/sum(f_bray_meta_MDS_soil_RH_20$eig),2),"%]"))+
  ylab(paste("MDS2","[",round(100*f_bray_meta_MDS_soil_RH_20$eig[2]/sum(f_bray_meta_MDS_soil_RH_20$eig),2),"%]"))

                    


f_mdsplot_2020 <- f_mdsplot 
f_mdsplot_2021=f_mdsplot_2021+labs(title = "Fungi:", 
subtitle = expression( paste("BMc Inoculation: R"^{2}, " = 32.8%****")))+
  theme(plot.subtitle = element_text(size=13, face="bold"))+
  theme(plot.title = element_text(size=13, face="bold"))

f_mdsplot_2021

f_mdsplot_2020=f_mdsplot_2020+labs(title="Fungi:" ,
subtitle=  expression(paste("BMc Inoculation: R"^{2}, " = 28.1%****")))+
  theme(plot.subtitle = element_text(size=13, face="bold"))+
  theme(plot.title = element_text(size=13, face="bold"))
f_mdsplot_2020

plots=ggarrange(both[[1]], f_mdsplot_2020,both[[2]],f_mdsplot_2021, align = "hv", common.legend = T, legend = "right",
                labels = c("", "2020  ", "", "2021  "), font.label = list(size=15),label.x =-0.1, label.y = 1.01 )
png(filename = "plot16SITSMS4.png", width = 3000, height = 2500, res=300)
print(plots)
dev.off()

png(filename = "f_Loreen_MDS_2020.png", width=550, height=450)
print(f_mdsplot_2020)
dev.off()
write.csv(as.matrix(bray_permanova_RH_20) , paste0( "/home/ioannis.kampouris",
                  "/DiControl_JKI_Phase3_repos/",
                  "LTE_2_years/Bac_PERMANOVA_2020.csv"))

write.csv(as.matrix(bray_permanova_RH_21) , paste0( "/home/ioannis.kampouris",
                                                    "/DiControl_JKI_Phase3_repos/",
                                                    "LTE_2_years/Bac_PERMANOVA_2021.csv"))


write.csv(as.matrix(f_bray_permanova_RH_20) , paste0( "/home/ioannis.kampouris",
                                                    "/DiControl_JKI_Phase3_repos/",
                                                    "LTE_2_years/Fungal_PERMANOVA_2020.csv"))

write.csv(as.matrix(f_bray_permanova_RH_21) , paste0( "/home/ioannis.kampouris",
                                                    "/DiControl_JKI_Phase3_repos/",
                                                    "LTE_2_years/Fungal_PERMANOVA_2021.csv"))
