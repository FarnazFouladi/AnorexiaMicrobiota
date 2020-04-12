#Farnaz Fouladi
#04-10-2020 
#This R code generates taxanomy count tables

rm(list=ls())

#Generating taxanomy tables
output<-"/Users/farnazfouladi/Google Drive/AnorexiaPaper11-11-19/paper/output/"
input<-"/Users/farnazfouladi/Google Drive/AnorexiaPaper11-11-19/paper/input/"
setwd(input)
myT_nonNormal<-read.table(paste0(input,"non-normalizedsvMeta.txt"),header=TRUE,sep="\t")
finishAbundanceIndex_nn<-which(colnames(myT_nonNormal)=="Sample")-1
count_table<-myT_nonNormal[,1:finishAbundanceIndex_nn]
meta<-myT_nonNormal[,(finishAbundanceIndex_nn+1):ncol(myT_nonNormal)]
taxa<-read.table(paste0(input,"taxForwardReads.txt"),header = TRUE,sep = "\t",na.strings = "NA")
source("/Users/farnazfouladi/CarrollMouseTransfer/Rcodes/MicrobiotaTransfer/functions.R")
taxaNames<-c("Phylum","Class","Order","Family","Genus")

for (t in taxaNames){
  
  t1<-getTaxaTable(count_table,taxa,t)
  t1_norm<-norm(t1)
  t1_normMeta<-cbind(t1_norm,meta)
  write.table(t1_normMeta,paste0(t,"_norm_table.txt"),sep = "\t",row.names = TRUE,quote = FALSE)
}

for (t in taxaNames){
  
  t1<-getTaxaTable(count_table,taxa,t)
  t1_Meta<-cbind(t1,meta)
  write.table(t1_Meta,paste0(t,"_NonNorm_table.txt"),sep = "\t",row.names = TRUE,quote = FALSE)
}

