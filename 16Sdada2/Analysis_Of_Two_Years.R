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
set.seed(01122024)


track1=read.csv(file = "/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/dada2_out/track_until_chimera.csv") 
track2=read.csv(file = "/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/dada2_out/track_chimera_removal.csv") 
track3=read.csv(file = "/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/dada2_out/track_cleaned_mit_chloroplast.csv") 
reads_follow_up= cbind(track1, track2, Cleaned=track3$.)
ggplot(reads_follow_up%>%select(-Ratio)%>%gather(key="step", value = "reads"),aes(x=reads))+geom_histogram()+facet_wrap(~step)

100*reads_follow_up$Cleaned/reads_follow_up$input
100*reads_follow_up$Chimeric/reads_follow_up$input


LTE_ASVs=read.csv(file = "/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/dada2_out/ASV_Table_LTE_2020_2021.csv") %>%
  select(ASV, everything())
#Remove mitochnondria 
rownames(LTE_ASVs)<-LTE_ASVs$ASV
Taxonomy=select(LTE_ASVs, ASV, Phylum, Class, Order, Family, Genus)

LTE_ASVs=LTE_ASVs
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
  mutate(Block=str_sub(Sample, -2, -2)  )%>% mutate(Root_Window=str_replace(Sample, ".*.RW.*.","RW"))
LTE_metadata$BMs[LTE_metadata$BMs!="BMs"]<-"Control"
LTE_metadata$Root_Window[LTE_metadata$Root_Window!="RW"]<-"Normal-Sample"
  
LTE_metadata$BMs[LTE_metadata$BMs!="BMs"]<-"Control"
LTE_metadata=LTE_metadata%>%na.omit()
library(parallel)

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
  message(paste("Using ", mc/10, " cores"))
  out <- mclapply(seq_len(nr), mc.cores = mc/10, function(i) {
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


 rarefaction_curve=  quickRareCurve(t(LTE_ASVs[, 10:ncol(LTE_ASVs)]))  
 names(rarefaction_curve) <- paste("Sample", 1:224, sep = "")
 # Coerce data into "long" form.
 protox <- mapply(FUN = function(x, y) {
   mydf <- as.data.frame(x)
   colnames(mydf) <- "value"
   mydf$species <- y
   mydf$subsample <- attr(x, "Subsample")
   mydf
 }, x = rarefaction_curve, y = as.list(names(rarefaction_curve)), SIMPLIFY = FALSE)
 
 xy <- do.call(rbind, protox)
 rownames(xy) <- NULL  # pretty
 
 
 
# Plot the rarefaction curve.
 rarecurve_plot= ggplot(xy, aes(x = subsample, y = value, color = species)) +
 theme_bw() +
 scale_color_discrete(guide = "none") +  # turn legend on or off
 geom_line() + ylab("ASV Richness")+xlab("Sequencing Depth")
 rarecurve_plot
 Sequencing_Depth=as.data.frame(colSums(LTE_ASVs[,10:ncol(LTE_ASVs)]))
ggplot(Sequencing_Depth,aes(x=`colSums(LTE_ASVs[, 10:ncol(LTE_ASVs)])`))+geom_histogram()
  minimum_seq=min(Sequencing_Depth)
  
pdf(  "rarecurve.pdf", width = 1000, height = 1000)
print(rarecurve_plot)
dev.off()

ASVs= t(LTE_ASVs[, 10:ncol(LTE_ASVs)])

Have_you_Runned_previously_RF<-"NO"
if (Have_you_Runned_previously_RF=="NO"){
  list_of_1000_rarefaction=mclapply(1:1000, function(i){ rrarefy(ASVs, minimum_seq)},  mc.cores = 3)
df= as.data.frame( list_of_1000_rarefaction[[1]])
for (i in 2:1000){
  df=df+list_of_1000_rarefaction[[i]]
  print(i)
}
Rarefied_Table=df/1000
}

write.table(Rarefied_Table, file="/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/LTE_2_years/16S/ASV_Rarefied_Table_1000_times.txt", sep = "\t")
write.csv(data_frame(Sample= LTE_metadata$ID, Sequences= Sequencing_Depth[LTE_metadata$ID,]), 
file = "/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/LTE_2_years/16S/SequencingDepth_ASV_table_two_years.csv")


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
  distribution_ASV_low[i]=list_of_1000_rarefaction[[i]][1,"ASV27786"]
  print(i)}
distribution_ASV_low= as.data.frame( unlist(distribution_ASV_low))
ggplot(distribution_ASV_low, aes(x = distribution_ASV_low$`unlist(distribution_ASV_low)`)) + geom_density() +
  geom_vline(xintercept = round( mean(distribution_ASV_low$`unlist(distribution_ASV_low)`)))+
  geom_vline(xintercept = round(median(distribution_ASV_low$`unlist(distribution_ASV_low)`)), linetype="dashed") 


########################
#MP 2021 Root window
LTE_metadata_RH_RW=LTE_metadata%>%filter(Year==21&Soil_Type=="RH"&Root_Window=="RW")
min(Sequencing_Depth[LTE_metadata_RH_RW$ID,])
Have_you_Runned_previously_RF<-"NO"

if (Have_you_Runned_previously_RF=="NO"){
  list_of_1000_rarefaction=mclapply(1:1000, function(i){ rrarefy(ASVs[LTE_metadata_RH_RW$ID,], min(Sequencing_Depth[LTE_metadata_RH_RW$ID,]))}, mc.cores = 3)
  df= as.data.frame( list_of_1000_rarefaction[[1]])
  for (i in 2:1000){
    df=df+list_of_1000_rarefaction[[i]]
    print(i)
  }
  MP_RW_Rarefied_Table=df/1000
}
write.csv(MP_RW_Rarefied_Table, file ="/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/MP_2021_Root_Windows/16S/Rarefied_1000_times_RW_ASV_table_2021.csv" )
write.csv(data_frame(Sample= LTE_metadata_RH_RW$ID, Sequences= Sequencing_Depth[LTE_metadata_RH_RW$ID,]), file = "/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/MP_2021_Root_Windows/16S/SequencingDepth_RW_ASV_table_2021.csv")

# 2021 General 
LTE_metadata_2021=LTE_metadata%>%filter(Year==21&Root_Window=="Normal-Sample")
min(Sequencing_Depth[LTE_metadata_2021$ID,])
Have_you_Runned_previously_RF<-"NO"

if (Have_you_Runned_previously_RF=="NO"){
  list_of_1000_rarefaction=mclapply(1:1000, function(i){ rrarefy(ASVs[LTE_metadata_2021$ID,], min(Sequencing_Depth[LTE_metadata_2021$ID,]))}, mc.cores = 30)
  df= as.data.frame( list_of_1000_rarefaction[[1]])
  for (i in 2:1000){
    df=df+list_of_1000_rarefaction[[i]]
    print(i)
  }
  LTE_2021_Rarefied_Table=df/1000
}


write.csv(LTE_2021_Rarefied_Table, file ="/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/all_treatments_2021/16S/Rarefied_1000_times_ASV_table_2021.csv" )
write.csv(data_frame(Sample= LTE_metadata_2021$ID, Sequences= Sequencing_Depth[LTE_metadata_2021$ID,]), file = "/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/all_treatments_2021/16S/SequencingDepth_ASV_table_2021.csv")

write.csv(Taxonomy, file ="/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/all_treatments_2021/16S/ASV_taxonomy.csv" )
write.csv(Taxonomy, file ="/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/MP_2021_Root_Windows/16S/ASV_taxonomy.csv" )
write.csv(Taxonomy, file ="/home/ioannis.kampouris/DiControl_JKI_Phase3_repos/LTE_2_years/16S/ASV_taxonomy.csv" )


#Tree creation via msa

