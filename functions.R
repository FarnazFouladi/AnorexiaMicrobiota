#Farnaz Fouladi
#04-10-2020 
##Functions for generating taxonomic count tables in Taxanomy script

getTaxaTable<-function(svTable,taxaTable,taxa){
  colnames(svTable)<-taxaTable[,taxa]
  svTable<-t(svTable)
  tab<-rowsum(svTable,group=rownames(svTable))
  tab<-t(tab)
  colnames(tab)[ncol(tab)]<-"others"
  return(tab)
}

norm<-function(table){
  table<-table[rowSums(table)>1000,]
  average<-sum(rowSums(table))/nrow(table)
  table<-sweep(table,1,rowSums(table),"/")
  table<-log10(table*average + 1)
  return(table)
}
