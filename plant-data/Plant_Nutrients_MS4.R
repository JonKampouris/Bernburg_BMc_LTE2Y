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
set.seed(172590)
##################################################################################
#Be warned, sometimes the path might have changed by the way I uploaded on Github#
##################################################################################
#Load the table with the Nutrient data
##############################################################################
LTE_a=  as.data.frame( 
  read_excel("LTE_2_years/LTE2020_M_central_2021_growth.xlsx", sheet = "Sheet1"))
LTE2=gather(LTE_a, -(colnames(LTE_a[,c(1:11, 22)])), key="Nutrient", 
            value = "Value")
LTE2$Group=paste0( LTE2$Nutrient)
LTE3=LTE2
###############################################################################
#Make some preselections for later
###############################################################################

LTE_a2=LTE_a
LTE_a2[,12:22]=scale( log10(LTE_a2[,12:22]))
LTE_a2$SDM=LTE_a2$Shoot_DW

###############################################################################
#Figures 1B, 1C, 1D: Shoot dry mass, Nutrient Shoot Concentrations & St.eq. Model
###############################################################################
file1=LTE_a # Label as file1 to make smoother modifications
file1$G=paste0(file1$BMc,"_",file1$Year) # make groups
#####################################################################################
#Perform Dunn's tests (pairwise wilcoxon tests with  Benjamini-Hochberg correction)#
####################################################################################
dunns= dunn.test::dunn.test(file1$Shoot_DW, file1$G, method = "BH")
comparisons=data.frame(comparisons=dunns$comparisons, p=dunns$P.adjusted)%>%
  separate(., "comparisons", into = c("row","col"), sep = " - ")
orderlist= aggregate( file1, Shoot_DW~G, max)
comparisons=full_join(comparisons, orderlist%>%mutate(row=G), by="row")%>%
  na.omit()
comparisons=comparisons[order(comparisons$Shoot_DW, decreasing = T),]
names1=paste0(comparisons$row,"-", comparisons$col)
p=comparisons$p
names(p)=names1
letters=  multcompView::multcompLetters(p,  
compare = "<", reversed =F)$Letters%>%as.data.frame()
colnames(letters)="L1"
letters$G=rownames(letters)
file1=full_join(file1,letters)
new_with_letters=(file1)
Fig1D= ggplot(new_with_letters%>%mutate(Int=stringr::str_replace(Int,"ensive","")), 
aes(y=(Shoot_DW), x=as.factor( paste0("20", Year)), 
fill=factor(BMc, levels = c("Ctrl","BMc") )))+
geom_bar(stat="summary", position = position_dodge(width = 1))+theme_bw()+
geom_errorbar(stat="summary", position = position_dodge2(width = 1))+
theme( strip.background = 
element_rect(fill="skyblue4"),  strip.text.x  =  element_text(angle=0, 
size = 22, colour="white", face="bold"),
strip.text=  element_text(angle=0, size = 22, colour="white", face="bold") )+
theme( axis.text.y = element_text( size=25, face = "bold", colour="black")) +
theme(axis.text.x = element_text(size = 35, face = "bold", angle = 0, colour="black")) +
theme(legend.text = element_text( size=35, face="bold", colour = "black")) + 
theme(legend.position = "right",
axis.title.y = element_text( size=35, face = "bold", colour="black")) + 
theme(legend.title =  element_text( size=55, face="bold", colour = "black")) +
theme(legend.key.size =    unit(1.5,"cm"))+
scale_fill_manual(name="",values=c( "#146eb4","#ff9900"))+
ylab( expression(bold(paste("SDM [g plant"^{-1},']'))))+xlab("")+
theme(axis.title.x = element_text(size = 35, face = "bold", colour="black"))+
geom_text(stat = "summary", 
aes(y=Shoot_DW+50, label=L1), size=9, position = position_dodge2(width = 0.9)) +
theme(legend.position = "top")

#Print Figure 1B
Fig1B
png(filename = "SDMplot.png", res=300, height = 3000, width = 2000)
Fig1B
dev.off()
###################################################################
#Linear Model with the whole Dataset
#####################################################################
f1=(lm( log10(Shoot_DW)~BMc*(Year)*Tillage*Int, data=LTE_a))
car::leveneTest( log10(Shoot_DW)~BMc*as.factor(Year)*Tillage*Int, data=LTE_a)
plot(density(resid(f1)))
df1= data.frame( anova(f1))
hist(resid(f1))
#Histogram indicates a slight deviance from normality
#Calculate the explained Variance
df1$R2=100*df1$Sum.Sq/sum(df1$Sum.Sq)
write.csv(file = "MS4_ANOVA_LM.csv", df1)
###############################################

#Conduct SEMs for concentration and content. Figure 1D
LTE_a2=LTE_a
for(i in c(
  "Concentration","Content"
)){
#Transform concentrations to content by multiplying with the dry mass
if(i=="Content" ){
LTE_a2[,12:21]*LTE_a2$Shoot_DW}
  LTE_a2[,12:22]=scale( log10(LTE_a2[,12:22]))
LTE_a2$SDM=LTE_a2$Shoot_DW
LTE_a3=LTE_a2
colnames1=colnames(LTE_a2)
colnames1=gsub("\\[.*.","",colnames1)
colnames(LTE_a3)=colnames1
LTE_a3$Year=paste0("Sampling_", LTE_a3$Year)

library(lavaan)
seme_formula <- '
SDM ~Fe+N+Mg+Ca+Mn+S+P+Zn+K+Cu
Fe~N+Mg+Ca+Mn+S+P+Zn+K+Cu
N~Mg+Ca+Mn+S+P+Zn+K+Cu
Mg~Ca+Mn+S+P+Zn+K+Cu
Ca~Mn+S+P+Zn+K+Cu
Mn~S+P+Zn+K+Cu
S~P+Zn+K+Cu
P~Zn+K+Cu
Zn~K+Cu
K~Cu'

seme= sem(seme_formula, data = LTE_a3, estimator="MLR")
sem_model=summary( seme)
sem_model$pe=sem_model$pe%>%na.omit()
sem_eq=data.frame(ARG=c(sem_model$pe$lhs), PR=c(sem_model$pe$rhs), Estimate=c(sem_model$pe$z),
                  pvalue=c(sem_model$pe$pvalue))%>%na.omit()#%>%filter(pvalue<0.05&PR!=ARG)

sem_datas=data.frame(ARG=c(sem_model$pe$lhs), PR=c(sem_model$pe$rhs), Estimate=c(sem_model$pe$z),
                  pvalue=c(sem_model$pe$pvalue))%>%na.omit()%>%filter(pvalue<0.05&PR!=ARG)

write.csv(file = paste0(i,"stuct_eq_model.csv"), as.data.frame(sem_datas))
covariates=as.data.frame( 
lavInspect(seme, what = "cor.all"))%>%
  rownames_to_column(var="PR")%>%gather(-PR, value = "cor", key="ARG")%>%
  filter(PR%in%c(sem_eq$PR)&ARG%in%c(sem_eq$ARG))
sem_eq= full_join(sem_eq, covariates)%>%na.omit()
sem_eq$R2=(sem_eq$cor)^2

library(igraph)
library(ggnetwork)
sem_model_stats=  summary(seme, rsquare = TRUE, fit.measures = TRUE, standardized = TRUE) 
sem_eq$pvalue2=sem_eq$pvalue
sem_eq$pvalue2[sem_eq$pvalue2<0.05]="p<0.05"
sem_eq$pvalue2[sem_eq$pvalue2!="p<0.05"]="p>0.05"

sem_g=graph_from_data_frame(as.data.frame(sem_eq%>%select(PR, everything())),directed = F)
gg= ggnetwork(sem_g, arrow.gap=0.05, layout=layout_with_fr(sem_g, niter = 1000))
gg$pvalue2[is.na( gg$pvalue2)]=""
sem_model_graph= ggplot(gg,aes(x = x, y = y,   xend = xend, linetype=pvalue2, yend = yend))+
  geom_point(size=20, colour="black")+
  scale_color_distiller(palette =
"RdBu",direction = 1, "St. Model Coefficient", limits = c(-1,1)*max(abs( na.omit( gg$Estimate))))+theme_blank() + 
   scale_linetype_manual( values = c(2,1,4), name="p-value")+
  ggtitle(paste0("RMSEA: ", round(sem_model_stats$fit["rmsea"][[1]],3), 
"\n","Chi-Square: ",round(  sem_model_stats$fit["chisq"][[1]],3 )))+theme(plot.title = element_text(size=20))+
  theme(legend.position = "right", 
legend.key.size = unit(2,"cm"), legend.text = element_text(size=30), 
legend.title = element_text(size=30))+geom_edges(alpha=12, curvature = 0, 
                         linewidth=2,
  aes(  colour=Estimate))+geom_nodelabel( size=7,aes(label=name))+ xlim(-0.2,1)
#SEM model Figure 1D
sem_model_graph+geom_edgelabel_repel(size=4, aes(label=round(100*R2,3)),  max.iter = 100000000000)
}

# Pairwise testing 
tests_=list()
for(j in unique(LTE2$Group)){
file1=filter(LTE2, Group==paste0(j))
print(paste0(j))  
t2= shapiro.test( resid(lm((file1$Value+1 )~file1$BMc)))
if(t2$p.value<0.05){
t1=t.test((file1$Value+1 )~file1$BMc)}else{
t2= shapiro.test( resid(lm(log10(file1$Value+1 )~file1$BMc)))
if(t2$p.value<0.05){
  t1=t.test(log10(file1$Value+1 )~file1$BMc)}else{t1=t.test(base::rank(file1$Value )~file1$BMc)}  

m1=data.frame(Group=paste0(j),t=paste0(t1$p.value),s=paste0(t2$p.value)) 
tests_=rbind(tests_, m1)
}}

write.csv(file = "pairwise_wilcox_tests.csv",tests_)
#################
#Figure 1C##########
##################
library(rstatix)
#Nutrients multiple testing
new_with_letters=list()
for(i in unique(paste0(LTE3$Nutrient))){
file1=filter(LTE3, Nutrient==paste0(i)) 
file1$G=paste0(file1$BMc,"_",file1$Year)
dunns= dunn.test::dunn.test(file1$Value, file1$G, method = "BH")
comparisons=data.frame(comparisons=dunns$comparisons, p=dunns$P.adjusted)%>%
  separate(., "comparisons", into = c("row","col"), sep = " - ")
orderlist= aggregate( file1, Value~G, max)
comparisons=full_join(comparisons, orderlist%>%mutate(row=G), by="row")%>%na.omit()
comparisons=comparisons[order(comparisons$Value, decreasing = T),]
names1=paste0(comparisons$row,"-", comparisons$col)
p=comparisons$p
names(p)=names1
letters=  multcompView::multcompLetters(p,  compare = "<", reversed =F)$Letters%>%as.data.frame()

colnames(letters)="L1"
letters$G=rownames(letters)
file1=full_join(file1,letters)
new_with_letters=rbind(new_with_letters,file1)
}

#New data with new letters for multiple groups testing
new_with_letters=as.data.frame(new_with_letters)
new_with_letters$L1[is.na(new_with_letters$L1)]=""
new_with_letters$BMc=factor(new_with_letters$BMc, levels = c("Ctrl", "BMc"))
new_with_letters$Year2=factor(new_with_letters$Year, levels = c(2021,2020))

    
new_with_letters$Nutrient2=factor(
  new_with_letters$Nutrient,
  levels = c("N[g/kg]", 
             "P[g/kg]",
             "Ca[g/kg]", 
             "K[g/kg]",
             "S[g/kg]",
             "Mg[g/kg]",
             "Cu[mg/kg]",
             "Fe[mg/kg]",
            "Mn[mg/kg]",
            "Zn[mg/kg]" ), 

labels = c( 
            "bold(N~g~kg^{-1}~SDM)", 
           "bold(P~g~kg^{-1}~SDM)",
           "bold(Ca~g~kg^{-1}~SDM)", 
           "bold(K~g~kg^{-1}~SDM)",
           "bold(S~g~kg^{-1}~SDM)",
           "bold(Mg~g~kg^{-1}~SDM)",
           "bold(Cu~mg~kg^{-1}~SDM)",
           "bold(Fe~mg~kg^{-1}~SDM)",
           "bold(Mn~mg~kg^{-1}~SDM)",
           "bold(Zn~mg~kg^{-1}~SDM)" ))



nutrients=ggplot(new_with_letters,
       aes(fill= BMc, x=factor(paste0("20", Year),
levels = c("2020", "2021")),
y=as.numeric(( Value))))+
  facet_wrap(~Nutrient2,strip.position = "left",nrow=2, scales = "free_y",
             labeller = label_parsed)+
  theme_bw()+
  theme( strip.background = 
          element_rect(fill="white", colour="white"),
         strip.text.x  =  element_text(angle=0, size = 22, colour="black", face="bold"),
        strip.text= element_text(angle=0, size = 35, colour="black", face="bold") )+
  geom_bar(stat="summary", position = position_dodge2(width = 5),
              colour="black")+
  geom_errorbar(stat="summary", position = position_dodge2(width = 13))+
  theme( axis.text.y = element_text( size=25, face = "bold", colour="black"),
         strip.placement = "outside") +
  theme(axis.text.x = element_text(size = 35, face = "bold", angle = 0, hjust = 0, vjust = 0, colour="black")) +
  theme(legend.text = element_text( size=35, face="bold", colour = "black")) + 
  theme(legend.position = "top",
axis.title.y = element_text( size=35, face = "bold", colour="black")) + 
  theme(axis.title.x = element_text(size = 35, face = "bold", colour="black")) + 
  theme(legend.title =  element_text( size=55, face="bold", colour = "black")) +
  theme(legend.key.size =    unit(1.5,"cm"))+
    scale_fill_manual(name="",values=c( "#146eb4","#ff9900"))+geom_text(stat = "summary", 
aes(y=Value*1.05, label=L1), size=8, vjust=-0.25, position = position_dodge2(width = 1)) +
  ylab( expression(bold(paste(""))))  +xlab("")


png(filename = "MS4Nutri.tiff", res=300, height = 4000, width = 9000)
print(nutrients)
dev.off()



############################################################################### 
#                                 End                                    ####  
#############################################################################







