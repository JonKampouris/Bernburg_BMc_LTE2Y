  library(ggplot2) #ggplot2: Fancy graphs
  library(ggpubr) #ggplot expansion for more fancy graphs
  library(readxl) # upload files  
  library(reshape2)# Chnange the shapes and apply casting
  library(readr) # read the files
  library(dplyr)# Data handling
  library(vegan)# Multivariate stat tool
  library(tidyr) # Data handling
  library(RColorBrewer) # Colours for fancy graphs
  library(tibble)# Data handling: rownames/columnanmes
  library(stringr) # Manipulate strings
  library(forcats)
  library(randomForest) # Random Forest model to train over samples

  set.seed(172590)
  
  source(paste0("/home/ioannis.kampouris/",
  "DiControl_JKI_Phase3_repos/",
  "MS3_correlation_network/functions_load.R"))
  
ASV_table= read.csv( file=paste0( "/home/ioannis.kampouris",
  "/DiControl_JKI_Phase3_repos/",
"LTE_2_years/16S/decontam_Rarefied_1000_times_ASV_table_2020_2021.csv"), row.names = 1)

#Create the metadata from the sample names.
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
metadata_RH=filter(metadata, Soil_Type=="RH"&Root_Window=="Normal-Sample"&ID!="RH.MP.Int.C1..2020.trimmed")
  
  #Creating the tables for the Random Forest Model
aggregation_table=cbind(t(ASV_table),Taxonomy[rownames(t(ASV_table)),])
#Combine every taxonomic rank
aggregation_table=aggregation_table %>%
  mutate(tax.Genus=paste(tax.Class,tax.Order,tax.Family,tax.Genus))
#Select the samples
aggregation_table2=as.data.frame( t(ASV_table[metadata_RH$ID,]))#%>%select(taxa=Acronym, everything())
input_for_the_rf_regression2=cbind(metadata_RH,(ASV_table[metadata_RH$ID,]))
LTE2020_M_central_2021_growth <- read_excel("LTE_2_years/LTE2020_M_central_2021_growth.xlsx", 
sheet = "Sheet1")%>%mutate(Fe=`Fe[mg/kg]`)%>%
  select(Int, Year, BMc, Tillage, Block,Fe)

input_for_the_rf_regression2$Year=as.double(input_for_the_rf_regression2$Year)
input_for_the_rf_regression3=  full_join(LTE2020_M_central_2021_growth,
input_for_the_rf_regression2, by=c("Int", "Year", "Block", "Tillage","BMc"))%>%na.omit()
rownames(input_for_the_rf_regression3)=input_for_the_rf_regression3$Sample
Fe1=input_for_the_rf_regression3
input_for_the_rf_regression4=input_for_the_rf_regression3[,12:ncol(input_for_the_rf_regression3)]
input_for_the_rf_regression4=apply(input_for_the_rf_regression4, 2,function(x) as.numeric(x))
var= as.data.frame( apply(input_for_the_rf_regression4,2, function (x) sum(x>0)))
colnames(var)="var"
var=filter(var, var>4)
input_for_the_rf_regression4=input_for_the_rf_regression4[, rownames(var)]

#Run a model against Iron concentrations.
linear_regressions=function(x){
  file1=lm(log10(x+0.1)~Fe1$Fe)
  file2=anova(file1)
  file3=data.frame(coeff=file1$coefficients[2], p=file2[1,5])
  return(file3)
}

input_for_the_rf_regression4=as.data.frame(input_for_the_rf_regression4)
linear_regressions_rf= apply( input_for_the_rf_regression4,2, function (x) linear_regressions(x))
linear_regressions_rf <- as.data.frame(do.call(rbind, linear_regressions_rf))
linear_regressions_rf[is.na(linear_regressions_rf)]<-1
linear_regressions_rf$padj=p.adjust(linear_regressions_rf$p, method="BH")
linear_regressions_rf2=(linear_regressions_rf)%>%rownames_to_column(var="Taxon")
linear_regressions_rf3=filter(linear_regressions_rf2, coeff>0&padj<0.05)
write.csv(

full_join(
  linear_regressions_rf3, mutate(Taxonomy, Taxon=ASV))%>%filter(padj!="NA"),
file =  "iron_associated_bacterial_ASVs_2020_2021.csv")

rf_bacteria= randomForest(scale(log10( Fe1$Fe))~., 
data=scale(log10(input_for_the_rf_regression4[,linear_regressions_rf3$Taxon]+1)),ntree=1000000, mtry=5, num.cores =20) 
rf_bacteria

important_bac_Fe=cbind(Taxonomy[rownames(rf_bacteria$importance),], rf_bacteria$importance,
                        linear_regressions_rf3)
important_bac_Fe_per=as.data.frame(table(important_bac_Fe$tax.Family))
important_bac_Fe_per2=important_bac_Fe_per[order(important_bac_Fe_per$Freq, decreasing = T),]
important_bac_Fe_per2=important_bac_Fe_per2[1:10,]
important_bac_Fe_per2=full_join(important_bac_Fe, important_bac_Fe_per2%>%mutate(tax.Family=Var1))%>%na.omit()

Bact_R2=100-100*(sum((rf_bacteria$y-rf_bacteria$predicted)^2)/sum((rf_bacteria$y-mean(rf_bacteria$y))^2))
no_of_ASVs= ggplot(important_bac_Fe_per2,
       aes(x=Freq,y=fct_reorder( tax.Family, Freq), fill=tax.Phylum))+geom_point(size=20, shape=21)+
  xlab("No. of Iron-Associated ASVs")+ylab("")+theme_pubr(legend = "right")+
  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=45, face = "bold", colour="black")) + 
  theme( plot.title =  element_text( size=35, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+scale_fill_brewer(palette = "Set1", name="")
no_of_ASVs=no_of_ASVs+ ggtitle(paste0("Cross validation with RF: Expl. Var=",round(Bact_R2,2), "%"))
no_of_ASVs

important_bac_Fe$coeff2=important_bac_Fe$coeff
important_bac_Fe$coeff2[important_bac_Fe$coeff2>0]<-"PS_Iron_Associated"
important_bac_Fe$coeff2[important_bac_Fe$coeff2!="PS_Iron_Associated"]<-"NG_Iron_Associated"
important_bac_Fe2=filter(important_bac_Fe, coeff2=="PS_Iron_Associated")


ASV_table_iron=select(ASV_table, c(important_bac_Fe2$ASV))
ASV_table_iron= cbind( ASV_table_iron[metadata_RH$ID,],metadata_RH)%>%
  gather(-c( colnames(metadata_RH)), key = "ASV",value = "RA")
ASV_table_iron$BMc=factor(ASV_table_iron$BMc, levels = c("Ctrl", "BMc"))  
ASV_table_iron= dcast(ASV_table_iron, BMc+Year+Block+Tillage+Int~., value.var = "RA", sum)
ASV_table_iron3=ASV_table_iron%>%full_join(Fe1%>%mutate(Year=as.character(Year)))


ASV_table_iron3$G=paste0(ASV_table_iron3$Tillage,"-", ASV_table_iron3$Int)
ASV_table_iron3$RA=100*ASV_table_iron3$`.`/20000




models_iron=list()
for(i in unique(ASV_table_iron3$G)){
  file1=filter(ASV_table_iron3, G==paste0(i))
  file2= lm(file1$Fe~log(RA), data = file1)
  file3=data.frame( sim=0.01:100)%>%mutate(G=paste0(i))
  (file2)
  file3$sim2=(file2$coefficients[1]+file2$coefficients[2]*log(file3$sim))
  file3$sim2_se2=(file2$coefficients[1]+(file2$coefficients[2]+coef(summary(file2))[, "Std. Error"][2])*(log(file3$sim)))
  file3$sim2_se3=(file2$coefficients[1]+(file2$coefficients[2]-coef(summary(file2))[, "Std. Error"][2])*(log(file3$sim)))
  models_iron=rbind(file3,models_iron)
  print(shapiro.test(resid(file2)))
  print(anova(file2))
}
models_iron_bacteria= ggplot(filter( models_iron), aes(y=sim, x=sim2, colour=G))+geom_smooth(stat="identity", size=2)+
geom_line(stat="identity", linetype=2,size=2, aes(y=sim, x=sim2_se2))+
  geom_line(stat="identity", linetype=2,size=2, aes(y=sim, x=sim2_se3))+
  scale_color_brewer(palette = "Set1", name="Group")+  theme_bw() +
  
  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=45, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
size = 2, linetype = "solid"))+ ylab("Relative Abundance (%)")+xlab(expression(bold( "Predicted Fe concentration (mg SDM kg"^{-1}~")")))+
  theme(legend.position = "right")
models_iron_bacteria

ggplot(ASV_table_iron, 
       aes(x=paste0("Year 20", Year),y=100*`.`/20000, fill=BMc))+theme_pubr()+
geom_boxplot()+stat_compare_means(method = "wilcox", label = "p.signif", size=20)+
theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
theme( axis.title.y = element_text( size=45, face = "bold", colour="black")) + 
theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
coord_flip()+xlab("RA(%) of Iron-Associated ASVs")+ylab("")+
  scale_fill_manual(name="",values=c("#146eb4", "#ff9900"))

file1=ASV_table_iron 
file1=select(file1,  RA=., everything())

file1$G=paste0(file1$BMc,"_",file1$Year)
dunns= pairwise.wilcox.test(file1$RA, file1$G, method="BH")
comparisons=as.data.frame(dunns$p.value)%>% rownames_to_column(var="row")%>%
  gather(-row, key=col, value = "p")%>%na.omit()
orderlist= dcast( file1,G~. ,value.var = "RA",  mean)%>%select(., mean=., G)
orderlist2= dcast( file1,G~. ,value.var = "RA",  max)%>%select(., max=., G)

comparisons=full_join(comparisons, orderlist%>%mutate(row=G), by="row")%>%na.omit()

comparisons=comparisons[order(comparisons$mean, decreasing = T),]

names1=paste0(comparisons$row,"-", comparisons$col)
p=comparisons$p
names(p)=names1
letters=  multcompView::multcompLetters(p,  compare = "<", reversed =F)$Letters%>%as.data.frame()

colnames(letters)="L1"
letters$G=rownames(letters)

file1=full_join(file1,letters)%>%full_join(orderlist[,c(1, ncol(orderlist))], by="G")%>%
  full_join(orderlist2)
new_with_letters=(file1)

new_with_letters=as.data.frame(new_with_letters)
new_with_letters$L1[is.na(new_with_letters$L1)]=""

ironplot2=ggplot(new_with_letters, 
                 aes(y=paste0("Year 20", Year), fill=BMc))+theme_pubr()+
  geom_boxplot(aes(x=100*(RA/20000)), position = position_dodge(width = 1))+geom_text(stat="summary",
                                                aes(x=(100*(max/20000)), label=L1),
                                                    position = position_dodge(width = 1),
                                                size=10,hjust=-0.5)+
  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=45, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
xlab("RA(%) of Iron-Associated ASVs")+ylab("")+
  scale_fill_manual(name="",values=c("#146eb4", "#ff9900"))
ironplot2


png("LTE_2_years/16S/Iron_Associated_Bacteria2.png", res=100,   width = 2000, height=500)
print(ironplot2) 
dev.off()

important_bac_Fe2=important_bac_Fe2[order(important_bac_Fe2$IncNodePurity, decreasing = T),]
important_bac_Fe3=important_bac_Fe2[1:30,,drop=F]

ggplot(important_bac_Fe3,
aes(y=fct_reorder(paste0(tax.Family,";", Acronym),coeff), x=100*IncNodePurity/max(important_bac_Fe3$IncNodePurity),
                      fill=paste0( tax.Phylum) ))+
  geom_point(shape=21, size=10, colour="black")+theme_pubr(legend = "right")+ylab("")+
  xlab("Importance based on RF model")+
  theme( axis.text.y = element_text( size=15, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=25, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=15, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
  scale_fill_brewer(palette = "Spectral", name="Phyla") + 
  theme(legend.background =  element_rect(linewidth = 1, colour="black"))
freq1= as.data.frame( table(important_bac_Fe3$tax.Genus))
ggplot(freq1, aes(y=Var1, x=Freq))+geom_bar(stat="summary") + xlab("Frequency of Taxa assoicated with iron")


linear_regressions_rf4=full_join(mutate(linear_regressions_rf3,ASV=Taxon), Taxonomy)%>%na.omit()
Comanadaceae_Iron= filter(important_bac_Fe2, tax.Family=="Comamonadaceae")%>%full_join(
  cbind( ASV_table[metadata_RH$ID,important_bac_Fe2$ASV],metadata_RH)%>%
  gather(-c( colnames(metadata_RH)), key = "ASV",value = "RA"), by="ASV")%>%na.omit()
Comanadaceae_Iron$Group=paste0(Comanadaceae_Iron$BMc,"-", Comanadaceae_Iron$Year)
Comanadaceae_Iron2=dcast(Comanadaceae_Iron, tax.Family+Year+BMc+Group+Block~., value.var = "RA", sum)
Comanadaceae_Iron2$G2=paste0(Comanadaceae_Iron2$BMc,Comanadaceae_Iron2$Year)
Comanadaceae_Iron2$log=log10(Comanadaceae_Iron2$.)
tests=rstatix::pairwise_t_test( Comanadaceae_Iron2,
                       log~G2, p.adjust.method = "BH")

comparisons=as.data.frame(tests[,c("group1",
"group2","p.adj")])
orderlist= dcast( Comanadaceae_Iron2,G2~. ,value.var = ".",  mean)%>%select(., mean=., G2)
orderlist2= dcast( Comanadaceae_Iron2,G2~. ,value.var = ".",  max)%>%select(., max=., G2)

comparisons=full_join(comparisons, orderlist%>%mutate(group1=G2), by="group1")%>%na.omit()

comparisons=comparisons[order(comparisons$mean, decreasing = T),]

names1=paste0(comparisons$group1,"-", comparisons$group2)
p=comparisons$p.adj
names(p)=names1
letters=  multcompView::multcompLetters(p,  compare = "<", reversed =F)$Letters%>%as.data.frame()

colnames(letters)="L1"
letters$G2=rownames(letters)


#########################
#Calculate means and SDs
#############################
Comanadaceae_Iron3=dcast(Comanadaceae_Iron, tax.Family+Year+BMc+Group+Block+tax.Genus~., value.var = "RA", sum)%>%
  select(., RA=`.`, everything())%>%mutate(RA=100*RA/20000)%>%
  dcast(., tax.Family+Year+BMc+tax.Genus~., value.var = "RA", mean)%>%mutate(G2=paste0(BMc,Year))%>%
  full_join(letters)%>%full_join(orderlist)

Comanadaceae_Iron4=dcast(Comanadaceae_Iron, tax.Family+Year+BMc+Group+Block~., value.var = "RA", sum)%>%
  select(., RA=`.`, everything())%>%mutate(RA=100*RA/20000)%>%
  dcast(., tax.Family+Year+BMc~., value.var = "RA", sd)%>% select(., sd=`.`, everything())

Comanadaceae_Iron5=dcast(Comanadaceae_Iron, tax.Family+Year+BMc+Group+Block~., value.var = "RA", sum)%>%
  select(., RA=`.`, everything())%>%mutate(RA=100*RA/20000)%>%
  dcast(., tax.Family+Year+BMc~., value.var = "RA", mean)%>%
  select(., Mean=`.`, everything())%>%
  full_join(Comanadaceae_Iron4)

maxes= dcast(Comanadaceae_Iron, tax.Family+Year+BMc+Group+Block+tax.Genus~., value.var = "RA", sum)%>%
  select(., RA=`.`, everything())%>%
  dcast(., tax.Family+Year+BMc+Group+tax.Genus+Block~., value.var = "RA", sum)%>%
  dcast(., tax.Family+Year+BMc+Group+Block~., value.var = ".", sum)%>%
  select(., RA=`.`, everything())%>%mutate(RA=100*RA/20000)%>%
  dcast(., tax.Family+Year+BMc+Group~., value.var = "RA", max)
full_join(Comanadaceae_Iron4,Comanadaceae_Iron3)
Comanadaceae_Iron4$min=Comanadaceae_Iron4$.-Comanadaceae_Iron4$sd
Comanadaceae_Iron4$max=Comanadaceae_Iron4$.+Comanadaceae_Iron4$sd

Comanadaceae_Iron6=full_join(Comanadaceae_Iron3,Comanadaceae_Iron5)%>%full_join(letters)
families_plot1= ggplot(Comanadaceae_Iron6, aes(y=factor( paste0(BMc,"_20", Year),
levels = c("Ctrl_2020","BMc_2020", "Ctrl_2021","BMc_2021")),x=.))+
  geom_bar(stat = "summary", aes( fill=tax.Genus), alpha=0.8)+geom_pointrange(aes(x=Mean, xmin=Mean-sd,
                                                                     xmax=Mean+sd),size=1, linewidth=2.5, colour="red")+
  theme_bw()+
  scale_fill_viridis_d(name="")+geom_text(aes(x=Mean+0.4+sd, label=L1), stat="summary", size=10)+

  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=15, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
  xlab("RA(%)")+ylab("")
families_plot1

plotsalplha= cowplot::plot_grid(no_of_ASVs,ironplot2,families_plot1, align = "hv", ncol = 1, labels = c("A)","B)","C)"),
                                label_size = 70)
png(filename = "test1.png",height = 2000, width = 2000)
print(plotsalplha)
dev.off()

multiple_glm= function(x){
  if(sum(x>0)>3){
  file1$Treatment2=file1$BMc
  file1$Treatment2[file1$Treatment2=="BMc"]<-1
  file1$Treatment2[file1$Treatment2=="Ctrl"]<-0
  file1$Treatment2=as.numeric(  file1$Treatment2)
  file1[,"ASV"]=x
  means= dcast(file1, BMc~., value.var = "ASV",mean)
  FC=log2((means[1,2]+1)/(means[2,2]+1))
  
  Table_glm=glm(file1$Treatment2~scale(x), family=binomial(link="logit"))
  file2=anova(Table_glm, test = "Chisq")
  file3=data.frame(coeff=(Table_glm$coefficients[2]), mean_RA=mean(x),FC=FC, p=file2$`Pr(>Chi)`[2])
  return(file3)
  }
  }
metadata_RH$Group=paste0(metadata_RH$Tillage,";",metadata_RH$Int, metadata_RH$Year )

diff_abundant_Bacterial_ASVs=list()
for(i in unique(metadata_RH$Group)){
  print(i)
file1=filter(metadata_RH, Group==paste0(i))
  multiple_glms= apply(as.data.frame(ASV_table[file1$ID,]),2, function(x) multiple_glm(x))
multiple_glms_df <- as.data.frame(do.call(rbind, multiple_glms))
multiple_glms_df[is.na(multiple_glms_df)]<-1
multiple_glms_df$padj=p.adjust(multiple_glms_df$p, method="BH")
multiple_glms_df=filter(multiple_glms_df)%>%mutate(Group=paste0(i))
multiple_glms_df=rownames_to_column(multiple_glms_df, var="ASV")
diff_abundant_Bacterial_ASVs=rbind(diff_abundant_Bacterial_ASVs,multiple_glms_df)
}
diff_abundant_Bacterial_ASVs2=full_join( important_bac_Fe[,c("ASV","coeff2")],
diff_abundant_Bacterial_ASVs, by="ASV")%>%full_join(Taxonomy, by="ASV")
diff_abundant_Bacterial_ASVs2$coeff2[is.na(diff_abundant_Bacterial_ASVs2$coeff2)]="Not-Associated with Fe"
diff_abundant_Bacterial_ASVs3=filter(diff_abundant_Bacterial_ASVs2,padj<0.05&coeff2=="PS_Iron_Associated")
diff_abundant_Bacterial_ASVs3$Year=gsub(".*.2","Year 202", diff_abundant_Bacterial_ASVs3$Group)
diff_abundant_Bacterial_ASVs3$G=gsub("2.*.","", diff_abundant_Bacterial_ASVs3$Group)
ASV_table_iron_bac_FC= ggplot(diff_abundant_Bacterial_ASVs3, aes(x=100*(mean_RA)/20000, y=FC,shape=G,fill=tax.Phylum))+
  geom_point(size=8)+facet_wrap(~Year,nrow = 1)+scale_shape_manual(values = c(21,22,23,24), name="")+
  scale_fill_brewer(palette = "Spectral", name="Phylum") +theme_bw()+
    theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
    theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
    theme( axis.title.y = element_text( size=45, face = "bold", colour="black")) +
   theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) + 
  theme( strip.text = element_text( size=45, face = "bold", colour="white")) + 
  theme( strip.background = element_rect(fill="skyblue4", colour="black",linewidth = 4)) + 
  
    theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
    theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
    theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
          size = 2, linetype = "solid"))+ ylab("log2FC(BMc/Control)")+xlab("RA (%)")+
    theme(legend.position = "right")+guides(fill=guide_legend(override.aes = list(shape=21, size=12)))+
  geom_hline(yintercept =0, linetype=2 )

png("LTE_2_years/16S/Iron_Associated_Bacteria_FC.png", res=100,   width = 2000, height=1000)
print(ASV_table_iron_bac_FC) 
dev.off()

diff_abundant_Bacterial_ASVs4=filter(diff_abundant_Bacterial_ASVs2,coeff2=="PS_Iron_Associated")

write.csv(diff_abundant_Bacterial_ASVs4,"LTE_2_years/16S/Iron_Associated_Bacteria_FC.csv")

#``````````````````````````````````````````````````````````````````````````
#Fungi
#``````````````````````````````````````````````````````````````````````````````
ASV_table= read.delim( file=paste0( "/home/ioannis.kampouris",
"/DiControl_JKI_Phase3_repos/",
"LTE_2_years/18S/fungal_ASV_Rarefied_Table_1000_times.txt"), sep = "\t", row.names = 1)%>%t()

Taxonomy= read.delim( file=paste0( "/home/ioannis.kampouris",
                                   "/DiControl_JKI_Phase3_repos/",
                                   "LTE_2_years/18S/ITS2_20202021_TAX_reduced.txt"), sep = "\t")
#Create numbered ASVs similar to the bacteria data
Taxonomy$ASV3=paste0("ASV", 1:nrow(Taxonomy))

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


aggregation_table2=as.data.frame( (ASV_table[,metadata_RH$ID]))#%>%select(taxa=Acronym, everything())
input_for_the_rf_regression2=cbind(metadata_RH,t(ASV_table[,metadata_RH$ID]))
LTE2020_M_central_2021_growth <- read_excel("LTE_2_years/LTE2020_M_central_2021_growth.xlsx", 
 sheet = "Sheet1")%>%select(Int, Year, BMc, Tillage, Block,Fe=`Fe[mg/kg]`)


input_for_the_rf_regression2$Year=as.double(input_for_the_rf_regression2$Year)
input_for_the_rf_regression3=  full_join(LTE2020_M_central_2021_growth,
                                         input_for_the_rf_regression2,
by=c("Int", "Year", "Block", "Tillage","BMc"))%>%na.omit()%>%as.data.frame()
rownames(input_for_the_rf_regression3)=input_for_the_rf_regression3$Sample
Fe1=input_for_the_rf_regression3
input_for_the_rf_regression4=input_for_the_rf_regression3[,12:ncol(input_for_the_rf_regression3)]
input_for_the_rf_regression4=apply(input_for_the_rf_regression4, 2,function(x) as.numeric(x))
var= as.data.frame( apply(input_for_the_rf_regression4,2, function (x) sum(x>0)))
colnames(var)="var"
var=filter(var, var>4)
input_for_the_rf_regression4=input_for_the_rf_regression4[, rownames(var)]


linear_regressions=function(x){
  file1=lm(log10(x+0.1)~Fe1$Fe)
  file2=anova(file1)
  file3=data.frame(coeff=file1$coefficients[2], p=file2[1,5])
  return(file3)
}

linear_regressions_rf= apply( input_for_the_rf_regression4,2, function (x) linear_regressions(x))
linear_regressions_rf <- as.data.frame(do.call(rbind, linear_regressions_rf))

linear_regressions_rf[is.na(linear_regressions_rf)]<-1
linear_regressions_rf$padj=p.adjust(linear_regressions_rf$p, method="BH")
linear_regressions_rf2=(linear_regressions_rf)%>%rownames_to_column(var="Taxon")
linear_regressions_rf3=filter(linear_regressions_rf2, coeff>0&padj<0.05)
input_for_the_rf_regression4=as.data.frame(input_for_the_rf_regression4)
linear_regressions_rf= apply( input_for_the_rf_regression4,2, function (x) linear_regressions(x))
linear_regressions_rf <- as.data.frame(do.call(rbind, linear_regressions_rf))
linear_regressions_rf[is.na(linear_regressions_rf)]<-1
linear_regressions_rf$padj=p.adjust(linear_regressions_rf$p, method="BH")
linear_regressions_rf2=(linear_regressions_rf)%>%rownames_to_column(var="Taxon")
linear_regressions_rf3=filter(linear_regressions_rf2, coeff>0&padj<0.05)
write.csv(
  
  full_join(
    linear_regressions_rf3, mutate(Taxonomy, Taxon=ASV))%>%filter(padj!="NA"),
  file =  "iron_associated_fungal_ASVs_2020_2021.csv")


rf_fungi= randomForest(scale(log10( Fe1$Fe))~., 
                          data=scale(log10(input_for_the_rf_regression4[,linear_regressions_rf3$Taxon]+1)),ntree=1000000, mtry=5, num.cores =20) 
rf_fungi
important_fungi_Fe=cbind(Taxonomy[(linear_regressions_rf3$Taxon),], rf_fungi$importance,
                       linear_regressions_rf3)


important_fungi_Fe$coeff2=important_fungi_Fe$coeff
important_fungi_Fe_per=as.data.frame(table(important_fungi_Fe$family))
important_fungi_Fe_per2=important_fungi_Fe_per[order(important_fungi_Fe_per$Freq, decreasing = T),]
important_fungi_Fe_per2=important_fungi_Fe_per2[1:10,]
important_fungi_Fe_per2=important_fungi_Fe_per2
important_fungi_Fe_per2=full_join(important_fungi_Fe, important_fungi_Fe_per2%>%mutate(family=Var1))%>%na.omit()

important_fungi_Fe_per2$family=gsub("_fam_*.*", "",important_fungi_Fe_per2$family)


Fungi_R2=100-100*(sum((rf_fungi$y-rf_fungi$predicted)^2)/sum((rf_fungi$y-mean(rf_fungi$y))^2))

no_of_fungal_ASVs=
ggplot(important_fungi_Fe_per2,
       aes(x=Freq,y=fct_reorder( family, Freq), fill=phylum))+geom_point(size=20, shape=21)+
  xlab("No. of Iron-Associated ASVs")+ylab("")+theme_pubr(legend = "right")+
  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=45, face = "bold", colour="black")) + 
    theme( plot.title =  element_text( size=35, face = "bold", colour="black")) + 

  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+scale_fill_brewer(palette = "Set1", name="")
no_of_fungal_ASVs=no_of_fungal_ASVs+ ggtitle(paste0("Cross validation with RF: Expl. Var=",round(Fungi_R2,2), "%"))

no_of_fungal_ASVs




important_fungi_Fe$coeff2[important_fungi_Fe$coeff2>0]<-"PS_Iron_Associated"
important_fungi_Fe$coeff2[important_fungi_Fe$coeff2!="PS_Iron_Associated"]<-"NG_Iron_Associated"
important_fungi_Fe2=filter(important_fungi_Fe, coeff2=="PS_Iron_Associated")
ASV_table_iron=select(t(ASV_table)%>%as.data.frame(), c(important_fungi_Fe2$ASV))
ASV_table_iron= cbind( ASV_table_iron[metadata_RH$ID,],metadata_RH)%>%
  gather(-c( colnames(metadata_RH)), key = "ASV",value = "RA")



ASV_table_iron$BMc=factor(ASV_table_iron$BMc, levels = c("Ctrl", "BMc"))  
ASV_table_iron= dcast(ASV_table_iron, BMc+Year+Block+Tillage+Int~., value.var = "RA", sum)
ASV_table_iron3=ASV_table_iron%>%full_join(Fe1%>%mutate(Year=as.character(Year)))
ASV_table_iron3$G=paste0(ASV_table_iron3$Tillage,"-", ASV_table_iron3$Int)
ASV_table_iron3$RA=100*ASV_table_iron3$`.`/48690 
models_iron=list()
for(i in unique(ASV_table_iron3$G)){
  file1=filter(ASV_table_iron3, G==paste0(i))
  file2= lm(file1$Fe~log(RA), data = file1)
  file3=data.frame( sim=0.01:100)%>%mutate(G=paste0(i))
  (file2)
  file3$sim2=(file2$coefficients[1]+file2$coefficients[2]*log(file3$sim))
  file3$sim2_se2=(file2$coefficients[1]+(file2$coefficients[2]+coef(summary(file2))[, "Std. Error"][2])*(log(file3$sim)))
  file3$sim2_se3=(file2$coefficients[1]+(file2$coefficients[2]-coef(summary(file2))[, "Std. Error"][2])*(log(file3$sim)))
  print(file2$coefficients[2])
  models_iron=rbind(file3,models_iron)
  print(shapiro.test(resid(file2)))
}
models_iron_fungi=ggplot(models_iron, aes(y=sim, x=sim2, colour=G))+geom_smooth(stat="identity", size=2)+
  geom_line(stat="identity", linetype=2,size=2, aes(y=sim, x=sim2_se2))+
  geom_line(stat="identity", linetype=2,size=2, aes(y=sim, x=sim2_se3))+
  scale_color_brewer(palette = "Set1", name="Group")+  theme_bw() +
  
  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=45, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
  theme(legend.position="right") + theme( axis.line = element_line(colour = "black", 
size = 2, linetype = "solid"))+ ylab("Relative Abundance (%)")+xlab("Predicted Fe concentration (ng/SDM g)")+
  theme(legend.position = "right")+xlab(expression(bold( "Predicted Fe concentration (mg SDM kg"^{-1}~")")))
models_iron_fungi

plots=ggarrange(models_iron_bacteria, models_iron_fungi,
                labels = c("A)", "B)"), font.label = list(size=45), common.legend = T, legend = "right", ncol = 1 )
plots

png(filename = "FigureS4.png", width = 9000, height = 9000, res=200)
print(plots)
dev.off()



ASV_table_iron_plot_f=ggplot(ASV_table_iron, 
                           aes(x=paste0("Year 20", Year),y=100*`.`/48690, fill=BMc))+theme_pubr()+
  geom_boxplot()+stat_compare_means(method = "wilcox", label = "p.signif", size=20)+
  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=45, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
  coord_flip()+ylab("RA(%) of Iron-Associated ASVs")+xlab("")+
  scale_fill_manual(name="",values=c("#146eb4", "#ff9900"))
ASV_table_iron_plot_f


linear_regressions_rf4=full_join(mutate(linear_regressions_rf3,ASV=Taxon), Taxonomy[,1:8] )%>%na.omit()
Hypocreaceae_Iron= filter(linear_regressions_rf4, family=="Hypocreaceae")%>%full_join(
  cbind( t(ASV_table)[metadata_RH$ID,linear_regressions_rf4$ASV],metadata_RH)%>%
    gather(-c( colnames(metadata_RH)), key = "ASV",value = "RA"), by="ASV")%>%na.omit()
Hypocreaceae_Iron2=dcast(Hypocreaceae_Iron, family+Year+BMc+Tillage+Int+Block~., value.var = "RA", sum)
Hypocreaceae_Iron2$G2=paste0(Hypocreaceae_Iron2$BMc,Hypocreaceae_Iron2$Year)
Hypocreaceae_Iron2$log=log10(Hypocreaceae_Iron2$.)
tests=rstatix::  pairwise_t_test( Hypocreaceae_Iron2,
                       log~G2, p.adjust.method = "BH")

comparisons=as.data.frame(tests[,c("group1",
                                   "group2","p.adj")])
orderlist= dcast( Hypocreaceae_Iron2,G2~. ,value.var = ".",  mean)%>%select(., mean=., G2)
orderlist2= dcast( Hypocreaceae_Iron2,G2~. ,value.var = ".",  max)%>%select(., max=., G2)

comparisons=full_join(comparisons, orderlist%>%mutate(group1=G2), by="group1")%>%na.omit()

comparisons=comparisons[order(comparisons$mean, decreasing = T),]

names1=paste0(comparisons$group1,"-", comparisons$group2)
p=comparisons$p.adj
names(p)=names1
letters=  multcompView::multcompLetters(p,  compare = "<", reversed =F)$Letters%>%as.data.frame()

colnames(letters)="L1"
letters$G2=rownames(letters)
Hypocreaceae_Iron3=dcast(Hypocreaceae_Iron, family+Year+BMc+Tillage+Int+Block+genus~., value.var = "RA", sum)%>%
  select(., RA=`.`, everything())%>%mutate(RA=100*RA/48690)%>%
  dcast(., family+Year+BMc+genus~., value.var = "RA", mean)%>%mutate(G2=paste0(BMc,Year))%>%
  full_join(letters)%>%full_join(orderlist)

Hypocreaceae_Iron4=dcast(Hypocreaceae_Iron, family+Year+BMc+Tillage+Int+Block~., value.var = "RA", sum)%>%
  select(., RA=`.`, everything())%>%mutate(RA=100*RA/48690)%>%
  dcast(., family+Year+BMc~., value.var = "RA", sd)%>% select(., sd=`.`, everything())

Hypocreaceae_Iron5=dcast(Hypocreaceae_Iron, family+Year+BMc+Tillage+Int+Block~., value.var = "RA", sum)%>%
  select(., RA=`.`, everything())%>%mutate(RA=100*RA/48690)%>%
  dcast(., family+Year+BMc~., value.var = "RA", mean)%>%
  select(., Mean=`.`, everything())%>%
  full_join(Comanadaceae_Iron4)

maxes= dcast(Hypocreaceae_Iron, family+Year+BMc+Tillage+Int+Block+genus~., value.var = "RA", sum)%>%
  select(., RA=`.`, everything())%>%
  dcast(., family+Year+BMc+Tillage+Int+genus+Block~., value.var = "RA", sum)%>%
  dcast(., family+Year+BMc+Tillage+Int+Block~., value.var = ".", sum)%>%
  select(., RA=`.`, everything())%>%mutate(RA=100*RA/48690)%>%
  dcast(., family+Year+BMc+Tillage+Int~., value.var = "RA", max)
Hypocreaceae_Iron4=full_join(Hypocreaceae_Iron4,Hypocreaceae_Iron3)
Hypocreaceae_Iron4$min=Hypocreaceae_Iron4$.-Hypocreaceae_Iron4$sd
Hypocreaceae_Iron4$max=Hypocreaceae_Iron4$.+Hypocreaceae_Iron4$sd

Hypocreaceae_Iron6=full_join(Hypocreaceae_Iron3,Hypocreaceae_Iron5)%>%full_join(letters)
fung_families_plot1= ggplot(Hypocreaceae_Iron6, aes(y=factor( paste0(BMc,"_20", Year),
levels = c("Ctrl_2020","BMc_2020", "Ctrl_2021","BMc_2021")),x=.))+
  geom_bar(stat = "summary", aes( fill=genus), alpha=0.8)+geom_pointrange(aes(x=Mean, xmin=Mean-sd,
xmax=Mean+sd),size=1, linewidth=2.5, colour="red")+
  theme_bw()+
  scale_fill_viridis_d(name="", direction = -1)+geom_text(aes(x=Mean+1.3+sd, label=L1), stat="summary", size=10)+
  
  theme( axis.text.y = element_text( size=45, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.text = element_text( size=45, face="bold", colour = "black")) + 
  theme( axis.title.y = element_text( size=15, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 45, face = "bold", colour="black")) +
  theme(legend.title =  element_text( size=45, face="bold", colour = "black"))+
  xlab("RA(%)")+ylab("")
fung_families_plot1
plotsalplha= cowplot::plot_grid(no_of_ASVs,no_of_fungal_ASVs,ironplot2, 
                                ASV_table_iron_plot_f2,families_plot1,
                                fung_families_plot1,
                                align = "hv", ncol = 2, 
                                labels = c("A)","B)","C)",
                                           "D)","E)","F)"),
                                label_size = 70)

plotsalplha

tiff(file = "Figure7.tiff",height = 7500, width = 8500, res=200)
print(plotsalplha)
dev.off()



####################################################################
##Save Comamonadacea ASVs as fasta for screning them against GTDB###
####################################################################

comamonadaceae_bac_Fe2= filter( important_bac_Fe2, tax.Family=="Comamonadaceae")
comamonadaceae_bac_Fe_string= Biostrings::DNAStringSet(
comamonadaceae_bac_Fe2$Sequence)
names(comamonadaceae_bac_Fe_string)=comamonadaceae_bac_Fe2$ASV
Biostrings::writeXStringSet(comamonadaceae_bac_Fe_string, 
filepath =paste0( "/home/ioannis.kampouris",
"/DiControl_JKI_Phase3_repos/",
"LTE_2_years/16S/iron-associated_comamonadaceae.fasta"),
format = "fasta")

