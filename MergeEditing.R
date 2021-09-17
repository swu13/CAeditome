#!/usr/bin/R
#Input File1:RNA editing file, row: RNA editing site, col:ENSG,GeneName,CAeditome,Tumor samples..
#Input File2:ADAR File, row: Tumor samples, col: RNA editing samples, RNA editing sample column, Expression samples, ADAR1, ADAR2, ADAR3
library(reshape2)
args=commandArgs(T)
filename=args[1]
outputName=args[2]
data=read.csv(filename,header=F,sep="\t")
colnames(data)=c('sample','Editing','Ratio')
Editing=dcast(data,formula=Editing~sample,mean)
write.table(Editing, file=outputName,sep="\t",quote=FALSE,row.names=FALSE)