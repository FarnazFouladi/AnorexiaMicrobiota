#Farnaz Fouladi
#04-10-2020
#Association between phenotypes and taxa using mixed linear model
#Model: taxa ~ Group*Phenotype+Week,  random = ~1 | Donor

rm(list=ls())
output<-"/Users/farnazfouladi/Google Drive/AnorexiaPaper11-11-19/paper/output/"
input<-"/Users/farnazfouladi/Google Drive/AnorexiaPaper11-11-19/paper/input/"
taxaNames<-c("Phylum","Class","Order","Family","Genus")
library(nlme)

for (t in taxaNames){
  
  setwd(paste0(output,t))
  pval_mixedModel_phenotype<-vector()
  pval_mixedModel_group<-vector()
  pval_mixedModel_time<-vector()
  pval_mixedModel_interaction<-vector()
  pval_mixedModel_groupT1_HC<-vector()
  pval_mixedModel_groupT2_HC<-vector()
 
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
              
              fit1<-anova(lme(bug ~ Group*Phenotype+Week, random = ~1 | Donor, data = myData, na.action = na.omit))
              fit2<-summary(lme(bug ~ Group*Phenotype+Week, random = ~1 | Donor, data = myData, na.action = na.omit))
              
              pval_mixedModel_group[index]<-fit1[2,4]
              pval_mixedModel_phenotype[index]<-fit1[3,4]
              pval_mixedModel_time[index]<-fit1[4,4]
              pval_mixedModel_interaction[index]<-fit1[5,4]
              
              pval_mixedModel_groupT1_HC[index]<-fit2$tTable[2,5]
              pval_mixedModel_groupT2_HC[index]<-fit2$tTable[3,5]
              
              phenotype[index]<-p
              bugName[index]<-colnames(myTM)[i]
              index<-index+1
              
            } else{
              
              myData<-myData[!is.na(myData$Phenotype),]
              
              fit1<-anova(lme(bug ~ Group*Phenotype, random = ~1 | Donor, data = myData, na.action = na.omit))
              fit2<-summary(lme(bug ~ Group*Phenotype, random = ~1 | Donor, data = myData, na.action = na.omit))
              
              pval_mixedModel_group[index]<-fit1[2,4]
              pval_mixedModel_phenotype[index]<-fit1[3,4]
              pval_mixedModel_time[index]<-NA
              pval_mixedModel_interaction[index]<-fit1[4,4]
              
              pval_mixedModel_groupT1_HC[index]<-fit2$tTable[2,5]
              pval_mixedModel_groupT2_HC[index]<-fit2$tTable[3,5]
              
              phenotype[index]<-p
              bugName[index]<-colnames(myTM)[i]
              index<-index+1
              
            }
        }
      }
    }
  
  
  df<-data.frame(bugName,phenotype,
                 pval_mixedModel_group,pval_mixedModel_phenotype,pval_mixedModel_time,
                 pval_mixedModel_interaction,pval_mixedModel_groupT1_HC,pval_mixedModel_groupT2_HC)
  df$Adjustedpval_mixedModel_group<-p.adjust(df$pval_mixedModel_group, method = "BH")
  df$Adjustedpval_mixedModel_phenotype<-p.adjust(df$pval_mixedModel_phenotype, method = "BH")
  df$Adjustedpval_mixedModel_interaction<-p.adjust(df$pval_mixedModel_interaction, method = "BH")
  df$Adjustedpval_mixedModel_time<-p.adjust(df$pval_mixedModel_time, method = "BH")
  df$Adjustedpval_mixedModel_groupT1_HC<-p.adjust(df$pval_mixedModel_groupT1_HC, method = "BH")
  df$Adjustedpval_mixedModel_groupT2_HC<-p.adjust(df$pval_mixedModel_groupT2_HC, method = "BH")
  write.table(df,paste0(t,"_GroupPhenotypeTimeModel_MixedLinearModel.txt"),sep="\t",row.names = FALSE)
}


#Family
t="Family"
myT<-read.table(paste0(input,t,"_norm_table.txt"), header=TRUE, sep="\t")
myT1<-myT[myT$Sample.type=="Mouse.feces",]
myT2<-myT1[myT1$Week==4,]

library(ggplot2)
library(ggsignif)
library(cowplot)
library(grid)
theme_set(theme_classic(base_size = 11))
col<-c("blue","darkorange2","forestgreen")
getplot<-function(data,taxa,phenotype,nameOfPhenotype){
  
  ggplot(data=data,aes(x=data[,phenotype],y=data[,taxa]))+geom_point(aes(color=Group))+
    scale_color_manual(values = col)+theme(legend.position = "none")+
    labs(y=expression(log[10]~"normalized count"),x=nameOfPhenotype,title=taxa)

}

plot1<-getplot(myT2,"Rikenellaceae","Cecum.wt","Cecum weight")+labs(subtitle="Group p=0.75\nPhenotype p<0.001\nInteraction p<0.001")
plot2<-getplot(myT2,"Family_XIII","Cecum.wt","Cecum weight")+labs(subtitle="Group p=0.99\nPhenotype p<0.001\nInteraction p=0.17")
plot3<-getplot(myT2,"Ruminococcaceae","Cecum.wt","Cecum weight")+labs(subtitle="Group p=0.75\nPhenotype p<0.001\nInteraction p<0.001")
plot4<-getplot(myT2,"Erysipelotrichaceae","Cecum.wt","Cecum weight")+labs(subtitle="Group p=0.99\nPhenotype p=0.007\nInteraction p=0.22")
plot5<-getplot(myT2,"Coriobacteriaceae","Cecum.wt","Cecum weight")+labs(subtitle="Group p=0.96\nPhenotype p=0.012\nInteraction p=0.08")
plot6<-getplot(myT2,"Acidaminococcaceae","Cecum.wt","Cecum weight")+labs(subtitle="Group p=0.99\nPhenotype p=0.043\nInteraction p=0.87")

plot7<-getplot(myT1,"Bifidobacteriaceae","Body.wt.pct.change","Change in body weight")+labs(subtitle="Group p=0.99\nPhenotype p=0.002\nInteraction p=0.13\nTime=0.004")
plot8<-getplot(myT1,"Enterococcaceae","Body.wt.pct.change","Change in body weight")+labs(subtitle="Group p=0.98\nPhenotype p=0.021\nInteraction p=0.37\nTime=0.16")
plot9<-getplot(myT1,"Christensenellaceae","Fat.mass.pct.change","Change in Fat mass")+labs(subtitle="Group p=0.97\nPhenotype p=0.002\nInteraction p=0.87\nTime=0.006")
plot10<-getplot(myT1,"Bifidobacteriaceae","Fat.mass.pct.change","Change in Fat mass")+labs(subtitle="Group p=0.98\nPhenotype p=0.024\nInteraction p=0.89\nTime<0.001")
plot11<-getplot(myT1,"Peptostreptococcaceae","Fat.mass.pct.change","Change in Fat mass")+labs(subtitle="Group p=0.99\nPhenotype p=0.024\nInteraction p=0.76\nTime=0.09")
plot12<-getplot(myT1,"Rikenellaceae","Fat.mass.pct.change","Change in Fat mass")+labs(subtitle="Group p=0.96\nPhenotype p=0.043\nInteraction p=0.36\nTime=0.84")

setwd(paste0(output,"Figures"))
png("PhenotypComparison.png",units="in", width=10, height=10,res=300)
plot_grid(plot1,plot2,plot3,
          plot4,plot5,plot6,
          plot7,plot8,plot9,plot10,plot11,plot12,ncol=4,nrow=3,scale = 0.8)
dev.off()

#Extract legend
getplotWithLegend<-function(data,taxa,phenotype,nameOfPhenotype){
  
  ggplot(data=data,aes(x=data[,phenotype],y=data[,taxa]))+geom_point(aes(color=Group))+
    scale_color_manual(values = col,labels=names)+
    labs(y=expression(log[10]~"normalized count"),x=nameOfPhenotype,title=taxa)
  
}
names<-c("non-AN","AN T1","AN T2")
myPlot<-getplotWithLegend(myT1,"Rikenellaceae","Cecum.wt","Cecum weight")+labs(subtitle="Group p=0.55\nPhenotype p<0.001\nInteraction p<0.001",
                                                                              col="Groups")
png("legend.png", units="in", width=8, height=8,res=300)
legendp <- cowplot::get_legend(myPlot)
grid.newpage()
grid.draw(legendp)
dev.off()

