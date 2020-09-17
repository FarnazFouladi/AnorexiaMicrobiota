#Farnaz Fouladi
#04-10-2020 
#This R code generates taxonomic count tables

rm(list=ls())

output<-"./output/"
input<-"./input/"
source("./Rcode/AnorexiaMicrobiota/functions.R")

myT_nonNormal<-read.table(paste0(input,"non-normalizedsvMeta.txt"),header=TRUE,sep="\t")
finishAbundanceIndex_nn<-which(colnames(myT_nonNormal)=="Sample")-1
count_table<-myT_nonNormal[,1:finishAbundanceIndex_nn]
meta<-myT_nonNormal[,(finishAbundanceIndex_nn+1):ncol(myT_nonNormal)]
taxa<-read.table(paste0(input,"taxForwardReads.txt"),header = TRUE,sep = "\t",na.strings = "NA")

taxaNames<-c("Phylum","Class","Order","Family","Genus")

for (t in taxaNames){
  
  t1<-getTaxaTable(count_table,taxa,t)
  t1_norm<-norm(t1)
  t1_normMeta<-cbind(t1_norm,meta)
  write.table(t1_normMeta,paste0(input,t,"_norm_table.txt"),sep = "\t",row.names = TRUE,quote = FALSE)
}

for (t in taxaNames){
  
  t1<-getTaxaTable(count_table,taxa,t)
  t1_Meta<-cbind(t1,meta)
  write.table(t1_Meta,paste0(input,t,"_NonNorm_table.txt"),sep = "\t",row.names = TRUE,quote = FALSE)
}

