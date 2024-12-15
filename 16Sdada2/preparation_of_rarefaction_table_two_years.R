library(ggplot2) #ggplot2: Fancy graphs
library(ggpubr) #ggplot expansion for more fancy graphs
library(readxl) # upload files  
library(readr) # read the files
library(dplyr)# Data handling
library(vegan)# Multivariate stat tool
library(tidyr) # Data handling
library(RColorBrewer) # Colours for fancy graphs
library(tibble)# Data handling: rownames/columnanmes
library(stringr)

library(Biostrings)
set.seed(17011990)
LTE_ASVs=read.csv(file = 
"/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/dada2_out/plastid_clean_low_read_clean_contaminant_free_ASV_Table_LTE_2020_2021.csv")
                
#Remove mitochondria 
rownames(LTE_ASVs)<-LTE_ASVs$ASV
Taxonomy=select(LTE_ASVs, ASV, Sequence, tax.Phylum, tax.Class, tax.Order, tax.Family, tax.Genus)

###Get the metadata out
LTE_metadata=data.frame(Sample=rownames(t(LTE_ASVs)), ID=rownames(t(LTE_ASVs)))%>%
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
  mutate(BMs= str_replace(Sample, ".*.t.I.*","BMs"))%>%
  mutate(BMs= str_replace(BMs, ".*.x.I.*","BMs"))%>%
  mutate(BMs= str_replace(BMs, ".*.n.I.*","BMs"))%>%
  mutate(Int= str_replace(Int, ".*Ex.*.","Extensive"))%>%
  mutate(Block=str_sub(Sample, -1, -1)  )%>% mutate(Root_Window=str_replace(Sample, ".*.RW.*.","RW"))
LTE_metadata$BMs[LTE_metadata$BMs!="BMs"]<-"Control"
LTE_metadata=LTE_metadata%>%na.omit()
library(parallel)
#Rarefaction curve
#rarefaction_curve=  rarecurve(t(LTE_ASVs_cleaned[, 10:ncol(LTE_ASVs_cleaned)]))
quickRareCurve <- function (x, step = 1, sample, xlab = "Sample Size",
                           ylab = "Species", label = TRUE, col, lty, max.cores = T, nCores = 1, ...)
{
  require(parallel)
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE))
    stop("function accepts only integers (counts)")
  if (missing(col))
    col <- par("col")
  if (missing(lty))
    lty <- par("lty")
  tot <- rowSums(x) # calculates library sizes
  S <- specnumber(x) # calculates n species for each sample
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  } # removes any empty rows
  nr <- nrow(x) # number of samples
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  # parallel mclapply
  # set number of cores
  mc <- getOption("mc.cores", ifelse(max.cores, detectCores(), nCores))
  message(paste("Using ", mc/3, " cores"))
  out <- mclapply(seq_len(nr), mc.cores = mc/3, function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i])
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab,
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample)
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"),
                                           y = z, xout = sample, rule = 1)$y)
    abline(h = rare, lwd = 0.5)
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) {
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}
metadata=data.frame(Sample=colnames(LTE_ASVs[, 11:ncol(LTE_ASVs)]), 
                    ID=colnames(LTE_ASVs[, 11:ncol(LTE_ASVs)]))%>%
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
metadata_RH=filter(metadata, Soil_Type=="RH"&Time_point!="T0"&Root_Window=="Normal-Sample"&ID!="RH.MP.Int.C1..2020.trimmed")


rarefaction_curve=  quickRareCurve(t(LTE_ASVs[, metadata_RH$ID]))

names(rarefaction_curve) <- paste("Sample", 1:ncol(LTE_ASVs[, metadata_RH$ID]), sep = "")
protox <- mapply(FUN = function(x, y) {
  mydf <- as.data.frame(x)
  colnames(mydf) <- "value"
  mydf$species <- y
  mydf$subsample <- attr(x, "Subsample")
  mydf
}, x = rarefaction_curve, y = as.list(names(rarefaction_curve)),
SIMPLIFY = FALSE)
xy <- do.call(rbind, protox)
rownames(xy) <- NULL  # pretty

rarecurve_plot= ggplot(xy, aes(x = subsample, y = value,colour=species)) +
  theme_bw() +
  scale_color_discrete(guide = "none") +  # turn legend on or off
  geom_line() + ylab("ASV Richness")+xlab("Sequencing Depth")
rarecurve_plot
png(file="rarecurve_plot16S_two.Years")
print(rarecurve_plot)
dev.off()
ASV_input=LTE_ASVs[, 11:ncol(LTE_ASVs)]
Sequencing_Depth=as.data.frame(colSums(ASV_input))

ggplot(Sequencing_Depth, 
       aes(x=Sequencing_Depth$`colSums(ASV_input)`)) +
  geom_histogram() +
  theme_bw()
#Change the coverage based on your script
Sequencing_Depth_25000=Sequencing_Depth%>%filter(.,`colSums(ASV_input)`>=20000)

ASVs=as.data.frame(t(ASV_input))
Have_you_Runned_previously_RF<-"NO"
if (Have_you_Runned_previously_RF=="NO"){
  list_of_1000_rarefaction=mclapply(1:1000, function(i){ rrarefy(ASVs[rownames(Sequencing_Depth_25000),], 20000)},  mc.cores = 10)
}

mean.dat <- as.data.frame(apply(simplify2array(lapply(list_of_1000_rarefaction, as.matrix)),1:2,mean))

write.table(mean.dat, file="/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/LTE_2_years/16S/Decontam_ASV_Rarefied_Table_1000_times.txt", sep = "\t")
distribution_ASV=list()
#Diagnostic for mean median and mode.

for (i in 1:1000){
  distribution_ASV[i]=list_of_1000_rarefaction[[i]][1,1]
  print(i)}
distribution_ASV= as.data.frame( unlist(distribution_ASV))
ggplot(distribution_ASV, aes(x = distribution_ASV$`unlist(distribution_ASV)`)) + geom_histogram() +
  geom_vline(xintercept =round(mean(distribution_ASV$`unlist(distribution_ASV)`)))+
  geom_vline(xintercept = round(median(distribution_ASV$`unlist(distribution_ASV)`)), linetype="dashed")

#The distribution is close to normal and the median, the mode and the mean are almost the same.
#This indicate that after rarefaction, this is the most probable value.

distribution_ASV_low=list()
#Diagnostic for mean median and mode.

for (i in 1:1000){
  distribution_ASV_low[i]=list_of_1000_rarefaction[[i]][1,"ASV9"]
  print(i)}
distribution_ASV_low= as.data.frame( unlist(distribution_ASV_low))
ggplot(distribution_ASV_low, aes(x = distribution_ASV_low$`unlist(distribution_ASV_low)`)) + density() +
  geom_vline(xintercept = round( mean(distribution_ASV_low$`unlist(distribution_ASV_low)`)))+
  geom_vline(xintercept = round(median(distribution_ASV_low$`unlist(distribution_ASV_low)`)), linetype="dashed") 


write.csv(mean.dat, file ="/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/LTE_2_years/16S/decontam_Rarefied_1000_times_ASV_table_2020_2021.csv" )
write.csv(Taxonomy, file ="/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/LTE_2_years/16S/ASV_taxonomy.csv" )
png(filename = "/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/LTE_2_years/16S/Rarefaction_curve.png", width=1000, height = 1000)
print(rarecurve_plot)
dev.off()
