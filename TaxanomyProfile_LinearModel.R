#Farnaz Fouladi
#04-10-2020
#Comparison of taxanomic compositions between non-AN, AN T1, and AN T2 using linear model

rm(list=ls())

output<-"/Users/farnazfouladi/Google Drive/AnorexiaPaper11-11-19/paper/output/"
input<-"/Users/farnazfouladi/Google Drive/AnorexiaPaper11-11-19/paper/input/"
taxaNames<-c("Phylum","Class","Order","Family","Genus","SV")

for (t in taxaNames){
  
  setwd(paste0(output,t))
  anovaResult<-vector()
  d<-vector()
  pval_groupT1HC<-vector()
  pval_groupT2HC<-vector()
  pval_groupT1T2<-vector()
  time<-vector()
  bugName<-vector()
  index<-1
  
  myT<-read.table(paste0(input,t,"_norm_table.txt"), header=TRUE, sep="\t")
  finishAbundanceIndex<-which(colnames(myT)=="Sample")-1
  myT1<-myT[myT$Sample.type=="Mouse.feces",]
  
  for (i in c(1:finishAbundanceIndex)){
    
    bug<-myT1[,i]
    
    if (mean(bug>0)>0.1){
      
      
      myData<-data.frame(bug,Group= myT1$Group,Week=myT1$Week,Donor=myT1$Donor)
      myData$Group<-relevel(myData$Group,ref = "HC")
      
      fit<-anova(lm(bug ~ Group + Week , data=myData))
      fit1<-summary(lm(bug ~ Group + Week , data=myData))
      
      
      anovaResult[index]<-fit[1,5]
      time[index]<-fit[2,5]
      
      pval_groupT1HC[index]<-fit1$coefficients[2,4]
      pval_groupT2HC[index]<-fit1$coefficients[3,4]
      
      myData$Group<-relevel(myData$Group,ref = "T1")
      fit1<-summary(lm(bug ~ Group + Week , data=myData))
      
      pval_groupT1T2[index]<-fit1$coefficients[3,4]
      
      d[index]<-fit1$r.squared
      bugName[index]<-colnames(myT1)[i]
      index<-index+1
    }
  }
  
  df<-data.frame(bugName,anovaResult,pval_groupT1HC,pval_groupT2HC,pval_groupT1T2,d)
  df$AdjustedAnovaResult<-p.adjust(df$anovaResult,method = "BH")
  df$Adjustedpval_groupT1HC<-p.adjust(df$pval_groupT1HC, method = "BH")
  df$Adjustedpval_groupT2HC<-p.adjust(df$pval_groupT2HC, method = "BH")
  df$Adjustedpval_groupT1T2<-p.adjust(df$pval_groupT1T2, method = "BH")
  write.table(df,paste0(t,"_Group_Comparison.txt"),sep="\t",row.names = FALSE)
}

#Family
t="Family"
myT<-read.table(paste0(input,t,"_norm_table.txt"), header=TRUE, sep="\t")
myT1<-myT[myT$Sample.type=="Mouse.feces",]
names<-c("non-AN","AN T1","AN T2")

result<-read.table(paste0(output,"Family/Family_Group_Comparison.txt"),sep="\t",header=TRUE)
result<-result[(order(result$Adjustedpval_groupT1HC)),]

library(ggplot2)
library(ggsignif)
library(cowplot)
theme_set(theme_classic(base_size = 12))
getplot<-function(data,taxa,d){
  
  ggplot(data=data,aes(x=Group,y=data[,taxa]))+geom_boxplot(outlier.shape = NA)+
    geom_jitter(shape=16, position=position_jitter(0.1))+labs(y=expression(log[10]~"normalized count"),x="",title=taxa,subtitle = bquote(R^2== ~.(d)) )+
    scale_x_discrete(labels=names)
}

plot1<-getplot(myT1,"Peptococcaceae",0.33)+
  geom_signif(y_position=c(3.5,3.5), xmin=c(1,2.1), xmax=c(1.9,3),annotation=c("***","***"), tip_length=0,vjust = 0.5,textsize =4)
plot2<-getplot(myT1,"Enterobacteriaceae",0.25)+
  geom_signif(y_position=c(3.9,4.2), xmin=c(1,1), xmax=c(2,3),annotation=c("***","***"), tip_length=0,vjust = 0.5,textsize =4)
plot3<-getplot(myT1,"Defluviitaleaceae",0.21)+
  geom_signif(y_position=c(1.7,1.9,2.1), xmin=c(1,1,2), xmax=c(2,3,3),annotation=c("***","***","***"), tip_length=0,vjust = 0.5,textsize =4)
plot4<-getplot(myT1,"Verrucomicrobiaceae",0.11)+
  geom_signif(y_position=c(4.5,4.5), xmin=c(1,2.1), xmax=c(1.9,3),annotation=c("***","***"), tip_length=0,vjust = 0.5,textsize =4)
plot5<-getplot(myT1,"Rikenellaceae",0.10)+
  geom_signif(y_position=c(3.9,4.2), xmin=c(1,1), xmax=c(2,3),annotation=c("***","***"), tip_length=0,vjust = 0.5,textsize =4)
plot6<-getplot(myT1,"Alcaligenaceae",0.08)+
  geom_signif(y_position=c(3.5,3.7,3.9), xmin=c(1,1,2), xmax=c(2,3,3),annotation=c("***","*","***"), tip_length=0,vjust = 0.5,textsize =4)
plot7<-getplot(myT1,"Christensenellaceae",0.15)+
  geom_signif(y_position=c(3,3.3), xmin=c(1,1), xmax=c(2,3),annotation=c("***","***"), tip_length=0,vjust = 0.5,textsize =4)
plot8<-getplot(myT1,"Veillonellaceae",0.10)+
  geom_signif(y_position=c(1.8,1.8), xmin=c(1,2.1), xmax=c(1.9,3),annotation=c("***","***"), tip_length=0,vjust = 0.5,textsize =4)
plot9<-getplot(myT1,"Ruminococcaceae",0.06)+
  geom_signif(y_position=c(4.3,4.6), xmin=c(1,1), xmax=c(2,3),annotation=c("***","***"), tip_length=0,vjust = 0.5,textsize =4)

setwd(paste0(output,"Figures"))
png("FamilyComparison.png",units="in", width=12, height=10,res=300)
plot_grid(plot1,plot2,plot3,
          plot4,plot5,plot6,plot7,
          plot8,plot9,
          ncol=3,nrow=3,scale = 0.9)
dev.off()

