#Farnaz Fouladi
#04-10-2020
#Association between phenotypes and taxa using linear model
#Model: taxa ~ Group*Phenotype+Week

rm(list=ls())

output<-"/Users/farnazfouladi/Google Drive/AnorexiaPaper11-11-19/paper/output/"
input<-"/Users/farnazfouladi/Google Drive/AnorexiaPaper11-11-19/paper/input/"
taxaNames<-c("Phylum","Class","Order","Family","Genus")

for (t in taxaNames){
  
  setwd(paste0(output,t))
  pval_phenotype<-vector()
  pval_group<-vector()
  pval_time<-vector()
  pval_interaction<-vector()
  d<-vector()
  pval_groupT1_HC<-vector()
  pval_groupT2_HC<-vector()
  pval_groupT2_T1<-vector()
  
  phenotype<-vector()
  bugName<-vector()
  index<-1
  
  myT<-read.table(paste0(input,t,"_norm_table.txt"), header=TRUE, sep="\t")
  finishAbundanceIndex<-which(colnames(myT)=="Sample")-1
  
  myTM<-myT[myT$Sample.type=="Mouse.feces",]
  
  for (i in c(1:finishAbundanceIndex)){
    
    bug<-myTM[,i]
    
    if (mean(bug>0)>0.1){
      
      for (p in c("Body.wt.pct.change","Fat.mass.pct.change","Lean.mass.pct.change","Cecum.wt","SI.wt","Gonadal.fat.wt")){
        
        
        myData<-data.frame(bug,Group= myTM$Group, Donor=myTM$Donor, Phenotype=myTM[,p],Week=myTM$Week)
        
        if(p!="Cecum.wt" & p!="SI.wt" & p!="Gonadal.fat.wt"){
        
        fit1<-anova(lm(bug ~ Group*Phenotype+Week, data = myData, na.action = na.omit))
        fit2<-summary(lm(bug ~ Group*Phenotype+Week, data = myData, na.action = na.omit))
        
        
        pval_group[index]<-fit1[1,5]
        pval_phenotype[index]<-fit1[2,5]
        pval_time[index]<-fit1[3,5]
        pval_interaction[index]<-fit1[4,5]
        
        d[index]<-summary(lm(bug ~ Group*Phenotype+Week, data = myData, na.action = na.omit))$r.squared
        pval_groupT1_HC[index]<-fit2$coefficients[2,4]
        pval_groupT2_HC[index]<-fit2$coefficients[3,4]
        
        myData$Group<-relevel(myData$Group,ref = "T1")
        fit2<-summary(lm(bug ~ Group*Phenotype+Week, data = myData, na.action = na.omit))
        pval_groupT2_T1[index]<-fit2$coefficients[3,4]
        
        phenotype[index]<-p
        bugName[index]<-colnames(myTM)[i]
        index<-index+1
        }else {
          
          myData<-myData[!is.na(myData$Phenotype),]
          
          fit1<-anova(lm(bug ~ Group*Phenotype, data = myData, na.action = na.omit))
          fit2<-summary(lm(bug ~ Group*Phenotype, data = myData, na.action = na.omit))
          
          
          pval_group[index]<-fit1[1,5]
          pval_phenotype[index]<-fit1[2,5]
          pval_time[index]<-NA
          pval_interaction[index]<-fit1[3,5]
          
          d[index]<-summary(lm(bug ~ Group*Phenotype, data = myData, na.action = na.omit))$r.squared
          pval_groupT1_HC[index]<-fit2$coefficients[2,4]
          pval_groupT2_HC[index]<-fit2$coefficients[3,4]
          
          myData$Group<-relevel(myData$Group,ref = "T1")
          fit2<-summary(lm(bug ~ Group*Phenotype, data = myData, na.action = na.omit))
          pval_groupT2_T1[index]<-fit2$coefficients[3,4]
          
          phenotype[index]<-p
          bugName[index]<-colnames(myTM)[i]
          index<-index+1
          
          
        }
      }
    }
  }
  
  
  df<-data.frame(bugName,phenotype,
                 pval_group,pval_phenotype,pval_time,
                 pval_interaction,pval_groupT1_HC,pval_groupT2_HC,pval_groupT2_T1,d)
  df$Adjustedpval_group<-p.adjust(df$pval_group, method = "BH")
  df$Adjustedpval_phenotype<-p.adjust(df$pval_phenotype, method = "BH")
  df$Adjustedpval_interaction<-p.adjust(df$pval_interaction, method = "BH")
  df$Adjustedpval_time<-p.adjust(df$pval_time, method = "BH")
  df$Adjustedpval_groupT1_HC<-p.adjust(df$pval_groupT1_HC, method = "BH")
  df$Adjustedpval_groupT2_HC<-p.adjust(df$pval_groupT2_HC, method = "BH")
  df$Adjustedpval_groupT2_T1<-p.adjust(df$pval_groupT2_T1, method = "BH")
  write.table(df,paste0(t,"_GroupPhenotypeTimeModel.txt"),sep="\t",row.names = FALSE)
}

