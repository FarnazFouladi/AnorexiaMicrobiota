#Farnaz Fouladi
#09-15-2020
#Comparison of taxonomic compositions between non-AN, AN T1, and AN T2 using 
#mixed linear model. taxa ~ group + week

rm(list=ls())

library(nlme)

output<-"./output/"
input<-"./input/"
taxaNames<-c("Phylum","Class","Order","Family","Genus","SV")

for (t in taxaNames){
  
  pval_groupT1HC<-vector()
  pval_groupT2HC<-vector()
  pval_groupT1T2<-vector()
  time<-vector()
  bugName<-vector()
  index<-1
  
  myT<-read.table(paste0(input,t,"_norm_table.txt"), header=TRUE, sep="\t")
  finishAbundanceIndex<-which(colnames(myT)=="Sample")-1
  myT1<-myT[myT$Sample.type=="Mouse.feces",]

  
  for (i in 1:finishAbundanceIndex){
    
    bug<-myT1[,i]
    
    if (mean(bug>0)>0.1){
      
      myData<-data.frame(bug,Group= myT1$Group,Donor=myT1$Donor,Week=myT1$Week)
      myData$Group<-relevel(factor(myData$Group), ref="HC")
      fit1<-summary(lme(bug~Group+Week,random = ~1 | Donor,data = myData ))
      
      pval_groupT1HC[index]<-fit1$tTable[2,5]
      pval_groupT2HC[index]<-fit1$tTable[3,5]
      bugName[index]<-colnames(myT1)[i]
      time[index]<-fit1$tTable[4,5]
      
      myData$Group<-relevel(myData$Group, ref="T1")
      fit1<-summary(lme(bug~Group+Week,random = ~1 | Donor,data = myData ))
      
      pval_groupT1T2[index]<-fit1$tTable[3,5]
      index<-index+1
      
    }
  }
  
  
  df<-data.frame(bugName,pval_groupT1HC,pval_groupT2HC,pval_groupT1T2,time)
  df$Adjustedpval_groupT1HC<-p.adjust(df$pval_groupT1HC, method = "BH")
  df$Adjustedpval_groupT2HC<-p.adjust(df$pval_groupT2HC, method = "BH")
  df$Adjustedpval_groupT1T2<-p.adjust(df$pval_groupT1T2, method = "BH")
  df$Adjustedtime<-p.adjust(df$time, method = "BH")
  write.table(df,paste0(output,t,"/",t,paste0("_Group_Comparison_MixedLinearModel.txt")),sep="\t",row.names = FALSE)
}
