#Farnaz Fouladi
#04-10-2020 
#This R code generates taxanomy bar plots

rm(list=ls())
output<-"/Users/farnazfouladi/Google Drive/AnorexiaPaper11-11-19/paper/output/"
input<-"/Users/farnazfouladi/Google Drive/AnorexiaPaper11-11-19/paper/input/"
setwd(paste0(output,"Figures"))
taxaNames<-c("Phylum","Class","Order","Family")
names<-c("non-AN","AN T1","AN T2")

myList<-list()
index<-1

library(reshape2)
library(ggplot2)
library(cowplot)

for (t in taxaNames){
  
  myT<-read.table(paste0(input,t,"_NonNorm_table.txt"), header=TRUE, sep="\t")
  myT_mouse<-myT[myT$Sample.type=="Mouse.feces",]
  myT_mouseWeek4<-myT_mouse[myT_mouse$Week==4,]
  
  finishAbundanceIndex<-which(colnames(myT_mouseWeek4)=="Sample")-1
  meta<-myT_mouseWeek4[,(finishAbundanceIndex+1):ncol(myT_mouseWeek4)]
  myT1<-myT_mouseWeek4[,1:finishAbundanceIndex]
  
  myT2<-sweep(myT1,1,rowSums(myT1),"/")
  myT3<-rowsum(myT2,group = meta$Group)
  myT4<-rbind(myT3[1,]/sum(meta$Group=="HC"),myT3[2,]/sum(meta$Group=="T1"),myT3[3,]/sum(meta$Group=="T2"))
  myT5_high<-myT4[,apply(myT4,2,max)>=0.01] #Selecting bugs with at least average relative abundance 0.01
  myT5_low<-myT4[,apply(myT4,2,max)<0.01]
  myT5_high$Others<-rowSums(myT5_low)
  myT5_high$Sample<-rownames(myT5_high)
  
  if(t=="Phylum"){
    myT_long <- melt(myT5_high, id.vars = "Sample", variable.name = "Phylum")
    plot<-ggplot(myT_long, aes(x = Sample, y = value, fill= Phylum)) + 
      geom_bar(stat = "identity" ,width = 0.5)+
      theme(legend.text=element_text(size=8),legend.key.size = unit(0.8,"line"))+labs(x="",y="Relative abundance")+
      scale_fill_manual(values = c("steelblue2","darkgoldenrod1","deeppink","mediumpurple3","seagreen3","coral"))+
      scale_x_discrete(labels=names)+scale_y_continuous(expand = c(0, 0)) 
    myList[[1]]<-plot
  }
  if(t=="Class"){
    myT_long <- melt(myT5_high, id.vars = "Sample", variable.name = "Class")
    plot<-ggplot(myT_long, aes(x = Sample, y = value, fill= Class)) + 
      geom_bar(stat = "identity" ,width = 0.5)+
      theme(legend.text=element_text(size=8),legend.key.size = unit(0.8,"line"))+labs(x="",y="Relative abundance")+
      scale_fill_manual(values = c("steelblue2","darkgoldenrod1","deeppink","mediumpurple3","seagreen3",
                                   "coral","chocolate4","darkgreen"))+
      scale_x_discrete(labels=names)+scale_y_continuous(expand = c(0, 0)) 
    myList[[2]]<-plot
  }
  
  if(t=="Order"){
    
    myT_long <- melt(myT5_high, id.vars = "Sample", variable.name = "Order")
    plot<-ggplot(myT_long, aes(x = Sample, y = value, fill= Order)) + 
      geom_bar(stat = "identity" ,width = 0.5)+
      theme(legend.text=element_text(size=8),legend.key.size = unit(0.8,"line"))+labs(x="",y="Relative abundance")+
      scale_fill_manual(values = c("steelblue2","darkgoldenrod1","deeppink","mediumpurple3","seagreen3",
                                   "coral","chocolate4","darkgreen"))+
      scale_x_discrete(labels=names)+scale_y_continuous(expand = c(0, 0)) 
    myList[[3]]<-plot
  }
  
  if (t=="Family"){
    
    myT_long <- melt(myT5_high, id.vars = "Sample", variable.name = "Family")
    plot<-ggplot(myT_long, aes(x = Sample, y = value, fill= Family)) + 
      geom_bar(stat = "identity" ,width = 0.5)+
      theme(legend.text=element_text(size=8),legend.key.size = unit(0.8,"line"))+labs(x="",y="Relative abundance")+
      scale_fill_manual(values = c("steelblue2","darkgoldenrod1","deeppink","mediumpurple3","seagreen3","coral",
                                   "chocolate4","darkgreen","gray", "darkorange","darkkhaki","darkolivegreen1"))+
      scale_x_discrete(labels=names)+scale_y_continuous(expand = c(0, 0)) 
    myList[[4]]<-plot
    
  }
}

pdf("Taxanomy.pdf",height = 9)
plot_grid(myList[[1]],myList[[2]],nrow = 2,ncol=1)
plot_grid(myList[[3]],myList[[4]],nrow = 2,ncol=1)
dev.off()

png("FamilyBarPlot.png",units="in", width=6, height=6,res=300)
myList[[4]]
dev.off()

