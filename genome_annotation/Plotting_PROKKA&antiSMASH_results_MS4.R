######
#Code for analyzing the results from Comparative Genomics.
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(reshape2)


############Upload all data from the alignment
file1=read_delim("allfna.tsv", delim = "\t" )%>%filter(Alignment_size_q...5>400)

#Upload the counted siderephore from antiSMASH prediction and obtaining via grep. 
filenames <- list.files(path = "antismash_counts",
                        pattern="txt")

#Create list of data frame names without the ".tsv" part 
fnames <-str_remove(filenames, pattern = "\\_siderophore.txt")

#Load all files and unify format.
counting_siderophores=list()
for(i in fnames){
  filepath <- file.path(paste0("antismash_counts/",i,"_siderophore.txt"))
  file2=read_delim(filepath,delim="\t",col_names = FALSE)
  file2$Sample=paste0(i, ".txt")
  counting_siderophores=rbind(counting_siderophores,file2)
  }

#Load all files and unify format
filenames <- list.files(path = "prokka_selection",
                        pattern="tsv")

#Create list of data frame names without the ".tsv" part 
fnames <-str_remove(filenames, pattern = "\\.tsv")
counting_prokka=list()
for(i in fnames){
  filepath <- file.path(paste0("prokka_selection/",i,".tsv"))
  file2=read_delim(filepath,delim="\t",col_names = FALSE)
  file2$Sample=paste0(i, ".txt")
  counting_prokka=rbind(counting_prokka,file2)
}

#Find the ACC-Deaminase genes by PROKKA.
counting_prokka2=as.data.frame(counting_prokka)%>%filter(str_detect(`X7`, "1-aminocyclopropane-1-carboxylate deaminase"))%>%filter(Sample%in%c(unique(file1$File)))
ACC_no= as.data.frame( table(counting_prokka2$Sample))  
colnames(ACC_no)=c("File", "ACC-Deaminase")
ACC_no=filter(ACC_no, File%in%(file1$File))
#Count siderophores.
counting_siderophores$Status=gsub("*.*receptor*.*", "Siderophore_Receptor",counting_siderophores$X1 )
counting_siderophores$Status[counting_siderophores$Status!="Siderophore_Receptor"]="Siderophore_Biosynthesis"
counting_siderophores$File=counting_siderophores$Sample
counting_siderophores=as.data.frame(counting_siderophores)%>%full_join(ACC_no)
#Iron dependent genes.
sidirition= full_join(file1, counting_siderophores)%>%dcast(Genome+ASV+Similarity+`ACC-Deaminase`~Status, value.var = "Status", length )
sidirition2=sidirition[,-(ncol(sidirition))]
sidirition3=sidirition2%>%gather(-c(Genome,ASV,Similarity), key="Status",value="No_of_CDS")
sidirition3$No_of_CDS[is.na(sidirition3$No_of_CDS)]=0
sidirition3$No_of_CDS=as.numeric(sidirition3$No_of_CDS)
sidirition3$No_of_CDS[sidirition3$No_of_CDS>0]="Presence"
sidirition3$No_of_CDS[sidirition3$No_of_CDS!="Presence"]="Absence" 
sidirition3$Genome= gsub("*\\.1","", sidirition3$Genome)
sidirition3$Genome= gsub(", whole gen*.*","", sidirition3$Genome)
#Last figure with simple comparative genomics.
ggplot(sidirition3,
aes(x=No_of_CDS , y=Genome))+geom_point(shape=21, fill="skyblue4", size=9)+facet_grid(~ASV~Status,space='free', scales = "free")+theme_bw()+
theme(strip.text.y.right =  element_text(angle = 0))+
theme(strip.text=  element_text(size = 15, colour = "white"))+theme(strip.background = element_rect(fill="skyblue4"))  + 
  xlab("Predicted Annotated Function")+
ylab("")+theme(axis.text.x =   element_text(angle = 0, size=20))+
  theme(axis.title.x =   element_text(angle = 0, size=20))+theme(axis.text.y =   element_text(angle = 0, size=15))



