#Farnaz Fouladi
#04-10-2020
#Comparison of phenotypes, shannon diversity, Bray-curtis distances between groups using either
#a linear model or a mixed linear model.

rm(list =ls())

#libraries
library(ggplot2)
library(vegan)
library(ggsignif)
library(cowplot)
library(nlme)

output<-"/Users/farnazfouladi/Google Drive/AnorexiaPaper11-11-19/paper/output/"
input<-"/Users/farnazfouladi/Google Drive/AnorexiaPaper11-11-19/paper/input/"

myT<-read.table(paste0(input,"svMeta.txt"),sep = "\t",header=TRUE,check.names=FALSE,na.strings = "NA" ,comment.char="")
finishAbundanceIndex<-which(colnames(myT)=="Sample")-1

names<-c("non-AN","AN T1","AN T2")
col<-c("blue","darkorange2","forestgreen")

myTM<-myT[myT$Sample.type=="Mouse.feces",]
myTM$Group<-relevel(myTM$Group, ref="HC")
myTM$DailyFood<-myTM$Cum.food.consumption/28

##Comparison of phenotypes between groups
#Mixed linear model
anova(lme(Body.wt.pct.change~Group+Week, random = ~1 | Donor,data=myTM,na.action = na.omit))
summary(lme(Body.wt.pct.change~Group+Week, random = ~1 | Donor,data=myTM,na.action = na.omit))
p.adjust(c(0.1618,0.3635,0.0338),method = "BH") #0.2427 0.3635 0.1014
anova(lme(Fat.mass.pct.change~Group+Week, random = ~1 | Donor,data=myTM,na.action = na.omit))
anova(lme(Lean.mass.pct.change~Group+Week, random = ~1 | Donor,data=myTM,na.action = na.omit))
anova(lme(DailyFood~Group+Week, random = ~1 | Donor,data=myTM,na.action = na.omit))

#Linear model
anova(lm(Body.wt.pct.change~Group+Week,data=myTM,na.action = na.omit))
myTM$Group<-relevel(myTM$Group, ref="T1") #For pairwise comparison, we can change the reference
summary(lm(Body.wt.pct.change~Group+Week,data=myTM,na.action = na.omit))
p.adjust(c(0.128,0.339,0.0126),method = "BH") # 0.1920 0.3390 0.0378
anova(lm(Fat.mass.pct.change~Group+Week,data=myTM,na.action = na.omit))
anova(lm(Lean.mass.pct.change~Group+Week,data=myTM,na.action = na.omit))
anova(lm(DailyFood~Group+Week,data=myTM,na.action = na.omit))

##Cecum weight, small intestin and gonadal fat at week 4
#Mixed linear model
myTM_4<-myTM[myTM$Week==4,]
anova(lme(Cecum.wt~Group, random = ~1 | Donor, data=myTM_4,na.action = na.omit))
summary(lme(Cecum.wt~Group, random = ~1 | Donor, data=myTM_4,na.action = na.omit))
#pvalue<-c("HC-T1","HC-T2","T1-T2")
pval<-c(0.0570,0.0478,0.9219)
adjustpval<-p.adjust(pval,method = "BH") #0.0855 0.0855 0.9219
anova(lme(SI.wt~Group, random = ~1 | Donor, data=myTM_4,na.action = na.omit))
anova(lme(Gonadal.fat.wt~Group, random = ~1 | Donor, data=myTM_4,na.action = na.omit))

#Linear model
anova(lm(Cecum.wt~Group, data=myTM_4,na.action = na.omit))
pairwise.t.test(myTM_4$Cecum.wt,myTM_4$Group,p.adjust.method = "BH")#0.00024 0.00013 0.78229
anova(lm(SI.wt~Group, data=myTM_4,na.action = na.omit))
anova(lm(Gonadal.fat.wt~Group, data=myTM_4,na.action = na.omit))


##Principal component analysis
myT1<-myT[myT$Sample.type=="Human.donor"|myT$Sample.type=="Mouse.feces",]
myT1$Sample_type<-factor(myT1$Sample.type,levels = c("Human.donor","Mouse.feces"))
myMDS<-capscale(myT1[,1:finishAbundanceIndex]~1,distance="bray")
percentVariance<-myMDS$CA$eig/sum(eigenvals(myMDS))*100
df<-data.frame(MDS1=myMDS$CA$u[,1],MDS2=myMDS$CA$u[,2],Sample_type=myT1$Sample_type,Donor=myT1$Donor,Group=myT1$Group)

df$Donor_newName<-sapply(as.character(df$Donor),function(x){
  if (substr(x,3,4)=="34" & substr(x,5,6)=="HC" ) return(paste0("non-AN","_","1"))
  else if (substr(x,3,4)=="40" & substr(x,5,6)=="HC" ) return(paste0("non-AN","_","2"))
  else if (substr(x,3,4)=="81" & substr(x,5,6)=="HC" ) return(paste0("non-AN","_","3"))
  else if (substr(x,3,4)=="70" & substr(x,6,7)=="HC" ) return(paste0("non-AN","_","4"))
  else if (substr(x,3,4)=="34" & substr(x,5,6)=="T1" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","5"))
  else if (substr(x,3,4)=="40" & substr(x,5,6)=="T1" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","6"))
  else if (substr(x,3,4)=="81" & substr(x,5,6)=="T1" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","7"))
  else if (substr(x,3,4)=="70" & substr(x,6,7)=="T1" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","8"))
  else if (substr(x,3,4)=="34" & substr(x,5,6)=="T2" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","5"))
  else if (substr(x,3,4)=="40" & substr(x,5,6)=="T2" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","6"))
  else if (substr(x,3,4)=="81" & substr(x,5,6)=="T2" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","7"))
  else if (substr(x,3,4)=="70" & substr(x,6,7)=="T2" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","8"))
  
})


plot1<-ggplot(data=df,aes(x=MDS1,y=MDS2))+geom_point(aes(col=factor(Group),shape=Sample_type))+
  scale_colour_manual(values=col[1:length(levels(factor(df$Group)))],name="Groups",labels=names)+
  scale_shape_manual(values=c(1,16),labels=c("Human fecal samples","Mouse fecal pellets"),name="Sample type")+
  labs(x=paste0("MDS1 (",format(percentVariance[1],digits = 4),"%)"),y=paste0("MDS2 (",format(percentVariance[2],digits = 4),"%)"))
  

myCol=c("red","blue","orchid","purple","black","lightblue","hotpink","cyan","pink","tan","darkgrey","gold")
plot2<-ggplot(data=df,aes(x=MDS1,y=MDS2))+geom_point(aes(col=factor(Donor_newName),shape=Sample_type))+
  scale_colour_manual(values=myCol[1:length(levels(factor(df$Donor_newName)))],name="Donors")+
  scale_shape_manual(values=c(1,16),labels=c("Human fecal samples","Mouse fecal pellets"),name="Sample type")+
  labs(x=paste0("MDS1 (",format(percentVariance[1],digits = 4),"%)"),y=paste0("MDS2 (",format(percentVariance[2],digits = 4),"%)"))+
  theme(legend.key.size = unit(1,"line"))
  
##Statitics on MDS1 and MD2 for mouse data
myT2<-myT[myT$Sample.type=="Mouse.feces",]
myMDS<-capscale(myT2[,1:finishAbundanceIndex]~1,distance="bray")
percentVariance<-myMDS$CA$eig/sum(eigenvals(myMDS))*100

#MDS1
#Mixed linear model
df<-data.frame(MDS1=myMDS$CA$u[,1],MDS2=myMDS$CA$u[,2],Group=myT2$Group,Time=as.numeric(myT2$Week),Donor=myT2$Donor)
df$Group<-relevel(df$Group, ref="HC") #When using summary function the reference can be changed for pairwise comparison
summary(lme(MDS1~Time+Group, random = ~1 | Donor, data=df))
anova(lme(MDS1~Time+Group, random = ~1 | Donor, data=df))
p.adjust(c(0.5708,0.5325,0.9529),method = "BH") # 0.8562 0.8562 0.9529  #p time =8205

#Linear model
df$Group<-relevel(df$Group, ref="HC")
anova(lm(MDS1~Time+Group, data=df))
summary(lm(MDS1~Time+Group, data=df))
p.adjust(c(4.02e-07,1.31e-09,0.3259),method = "BH") #6.030e-07 3.930e-09 3.259e-01

#MDS2
#Mixed linear model
summary(lme(MDS2~Time+Group, random = ~1 | Donor, data=df))
anova(lme(MDS2~Time+Group, random = ~1 | Donor, data=df))
p.adjust(c(0.0670,0.0548,0.9042),method = "BH") #[1] 0.1005 0.1005 0.9042 #p time =0.3920 

#Linear model
anova(lm(MDS2~Time+Group, data=df))
summary(lm(MDS2~Time+Group, data=df))
p.adjust(c(2e-16,2e-16,0.236),method = "BH") # 3.00e-16 3.00e-16 2.36e-01

setwd(paste0(output,"Figures"))
theme_set(theme_classic(base_size = 12))
pdf("PcoPlots.pdf",width = 5.5,height = 3)
plot_grid(plot1,plot2,ncol=2,nrow=1)
dev.off()

##Nray-curtis distances
myT2_week4<-myT2[myT2$Week==4,]
myT3<-myT2_week4[,1:finishAbundanceIndex]
rownames(myT3)<-myT2_week4$Sample
distances<-vegdist(myT3,distance="bray")
distanes1<-as.matrix(distances)

HC<-distanes1[myT2_week4$Group=="HC",myT2_week4$Group=="HC"]
T1<-distanes1[myT2_week4$Group=="T1",myT2_week4$Group=="T1"]
T2<-distanes1[myT2_week4$Group=="T2",myT2_week4$Group=="T2"]

withinHCdistance<-colSums(HC)/(nrow(HC)-1)
withinT1distance<-colSums(T1)/(nrow(T1)-1)
withinT2distance<-colSums(T2)/(nrow(T2)-1)

HCT1<-distanes1[myT2_week4$Group=="T1",myT2_week4$Group=="HC"]
HCT2<-distanes1[myT2_week4$Group=="T2",myT2_week4$Group=="HC"]
T1T2<-distanes1[myT2_week4$Group=="T2",myT2_week4$Group=="T1"]

betweenHCT1distance<-colMeans(HCT1)
betweenT1HCdistance<-rowMeans(HCT1)

betweenHCT2distance<-colMeans(HCT2)
betweenT2HCdistance<-rowMeans(HCT2)

betweenT1T2distance<-colMeans(T1T2)
betweenT2T1distance<-rowMeans(T1T2)

#Get the id 
myT2_week4_HC<-myT2_week4[myT2_week4$Group=="HC",] #sum(colnames(HC)==as.character(myT2_week4_HC$Sample)) 50
myT2_week4_T1<-myT2_week4[myT2_week4$Group=="T1",] #50
myT2_week4_T2<-myT2_week4[myT2_week4$Group=="T2",] #53


d=c(betweenHCT1distance,betweenT1HCdistance,
    betweenHCT2distance,betweenT2HCdistance,
    betweenT1T2distance,betweenT2T1distance,
    withinHCdistance,withinT1distance,withinT2distance)
Group=c(rep("non-AN vs. AN T1",length(betweenHCT1distance)+length(betweenT1HCdistance)),
        rep("non-AN vs. AN T2",length(betweenHCT2distance)+length(betweenT2HCdistance)),
        rep("AN T1 vs. AN T2",length(betweenT1T2distance)+length(betweenT2T1distance)),
        rep("non-AN vs. non-AN",length(withinHCdistance)),
        rep("AN T1 vs. AN T1",length(withinT1distance)),
        rep("AN T2 vs. AN T2",length(withinT2distance)))
distance<-c(rep("Between",306),rep("Within",153))
ID<-c(myT2_week4_HC$Donor,myT2_week4_T1$Donor,
      myT2_week4_HC$Donor,myT2_week4_T2$Donor,
      myT2_week4_T1$Donor,myT2_week4_T2$Donor,
      myT2_week4_HC$Donor,myT2_week4_T1$Donor,myT2_week4_T2$Donor)


df<-data.frame(d,Group,distance,ID) 
df$Group<-factor(df$Group,levels = c("non-AN vs. AN T1","non-AN vs. AN T2","AN T1 vs. AN T2",
                                     "non-AN vs. non-AN","AN T1 vs. AN T1","AN T2 vs. AN T2"))
theme_set(theme_classic(base_size = 16))
pdf("BrayCurtisDistances.pdf",width=7,height=4)
plot<-ggplot(data=df,aes(x=factor(distance),y=d))+geom_boxplot(aes(col=factor(Group)),outlier.shape = NA)+
  geom_jitter(shape=16,size =0.8,position=position_jitterdodge(),aes(col=factor(Group)))+
  labs(y="Average Bray-Curtis Distances",x="",col="")+scale_x_discrete(labels=c("Between Groups","Within Groups"))+
  geom_signif(y_position=c(0.8,0.815,0.72,0.735,0.75), xmin=c(0.75,0.95,1.75,1.75,1.95), xmax=c(1.15,1.15,1.95,2.15,2.15),annotation=c("***","***","***","***","***"), tip_length=0,vjust = 0.5,textsize =4)
print(plot)
dev.off()
write.table(df,"Bray-CurtisDistances.txt",sep="\t",quote =FALSE)

#Between
bt<-df[df$distance=="Between",]
#Linear model
anova(lm(d~Group,data=bt))
pairwise.t.test(bt$d, bt$Group,p.adjust.method = "BH")
#Mixed linear model
anova(lme(d~Group,random = ~1 | ID,data=bt))
bt$Group<-relevel(bt$Group,ref="non-AN vs. AN T2")
summary(lme(d~Group,random = ~1 | ID,data=bt))
p.adjust(c(0.3184,0.0000,0.0000),method = "BH")

mean(bt[bt$Group=="non-AN vs. AN T1",]$d)
sd(bt[bt$Group=="non-AN vs. AN T1",]$d)

mean(bt[bt$Group=="non-AN vs. AN T2",]$d)
sd(bt[bt$Group=="non-AN vs. AN T2",]$d)

mean(bt[bt$Group=="AN T1 vs. AN T2",]$d)
sd(bt[bt$Group=="AN T1 vs. AN T2",]$d)

#Within
wt<-df[df$distance=="Within",]
#Linear model
anova(lm(d~Group,data=wt))
pairwise.t.test(wt$d, wt$Group,p.adjust.method = "BH")
#Mixed linear model
anova(lme(d~Group,random = ~1 | ID,data=wt))
wt$Group<-relevel(wt$Group,ref="AN T1 vs. AN T1")
summary(lme(d~Group,random = ~1 | ID,data=wt))
p.adjust(c(0.0118,0.0004,0.0448),method = "BH") #[1] 0.0177 0.0012 0.0448

##Shannon Diversity
myT_nonNormal<-read.table(paste0(input,"non-normalizedsvMeta.txt"),header=TRUE,sep="\t")
finishAbundanceIndex_nn<-which(colnames(myT_nonNormal)=="Sample")-1
myT_nonNormal$div<-diversity(myT_nonNormal[,1:finishAbundanceIndex_nn],index = "shannon")
myT_nonNormal1<-myT_nonNormal[myT_nonNormal$Sample.type=="Mouse.feces",]
myT_nonNormal1$Week<-as.numeric(myT_nonNormal1$Week)
myT_nonNormal1$Group<-relevel(myT_nonNormal1$Group,ref="HC")

#Linear mixed model
anova(lme(div~Group+Week, random = ~1 | Donor, data=myT_nonNormal1,na.action = na.omit))
summary(lme(div~Group+Week, random = ~1 | Donor, data=myT_nonNormal1,na.action = na.omit))
p.adjust(c(0.0581,0.0411,0.8379),method = "BH") #0.08715 0.08715 0.83790

#Linear model
myT_nonNormal1$Group<-relevel(myT_nonNormal1$Group,ref="T1")
anova(lm(div~Group+Week, data=myT_nonNormal1,na.action = na.omit))
summary(lm(div~Group+Week, data=myT_nonNormal1,na.action = na.omit))
p.adjust(c(2e-16,2e-16,0.368   ),method = "BH") #3.00e-16 3.00e-16 3.68e-01

#paiwise comparison at each time point following colonization
summary(lm(div~Group, data=myT_nonNormal1[myT_nonNormal1$Week==4,],na.action = na.omit))
#Week1
p.adjust(c(0.000234,0.000161,0.948257),method = "BH") #0.000351 0.000351 0.948257
#Week2
p.adjust(c(0.000228,6.72e-06,0.380948),method = "BH") #0.00034200 0.00002016 0.38094800
#Week3
p.adjust(c(5.36e-06,2.01e-07,0.508 ),method = "BH") #8.04e-06 6.03e-07 5.08e-01
#Week4
p.adjust(c(6.44e-08,1.39e-08,0.817),method = "BH") #9.66e-08 4.17e-08 8.17e-01

