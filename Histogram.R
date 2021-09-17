#!/usr/bin/R
#Input: 
#parameter1:Input
#parameter2:Output
library(ggplot2)
args=commandArgs(T)
InputFile=args[1]
OutputFile=args[2]
FigureCol=as.numeric(args[3])
FigureRow=as.numeric(args[4])
data=read.table(InputFile,sep='\t',header=FALSE)
colnames(data)=c("Name","Group","Value")
png(OutputFile,res = 300,width = FigureCol*500, height = FigureRow*500)
print(ggplot(data,aes(Value),log= "y")+geom_histogram(binwidth = 0.5,fill="Navy")+scale_y_continuous(trans="log",breaks=c(7,148,2981,59874))+facet_wrap(~Group, ncol = FigureCol,as.table=TRUE))
dev.off()