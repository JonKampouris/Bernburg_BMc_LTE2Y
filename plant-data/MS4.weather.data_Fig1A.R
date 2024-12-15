#The following code is about generating the Figure 1A.
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
set.seed(122022)

############################################################################### 
#Load the table with the Weather data
##############################################################################
LTE_a=  as.data.frame( 
  read_excel("MS4_Weather.xlsx"))

###############################################################################
# Figure 1A
###############################################################################

ghoutline <- data.frame(x = c(4,4,7.9,7.9,4),
                        y = c(1,95, 95, 1,1))

inoculations <- data.frame(x = c(5,6),
                           y = c(70,85))


Fig1A= ggplot(LTE_a , 
              aes(y=(mm), x=Month))+
  geom_bar(aes( fill=factor(Year)),
           stat="summary", position = position_dodge(width = 1),
           colour="black")+theme_bw()+
  
  geom_line(stat="summary",
            aes(y=`1981-2010`, x=Month, colour="Average (1981-2010)" ), linewidth=2)+
  scale_colour_manual(values = "black", name="")+
  geom_point(stat="summary", fill="white", colour="black", shape=21, size=5,
             aes(y=`1981-2010`, x=Month ))+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12), labels = c(
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec"
    
  ))+scale_fill_brewer(palette = "Paired", name="")+
  
  theme( axis.text.y = element_text( size=35, face = "bold", colour="black")) +
  theme(axis.text.x = element_text(size = 35, face = "bold", angle = 0, colour="black")) +
  theme(legend.text = element_text( size=35, face="bold", colour = "black")) + 
  theme(legend.position = "right",
        axis.title.y = element_text( size=35, face = "bold", colour="black")) + 
  theme(legend.title =  element_text( size=55, face="bold", colour = "black")) +
  theme(legend.key.size =    unit(1.5,"cm"))+
  theme(axis.title.x = element_text(size = 35, face = "bold", colour="black"))+
  ylab("Precipitation (mm)")+
  xlab("Month") +geom_path(data=ghoutline , aes(x=x,y=y), colour="gold", size=1.5)+
  geom_point(data=inoculations, aes(x=x, y=y), shape=23, size=17, fill="gold")

Fig1A

png(filename = "weatherplot.png", res=300, height = 4000, width = 9000)
Fig1A
dev.off()
