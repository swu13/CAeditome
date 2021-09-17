#!/usr/bin/R
#Input: 
#parameter1:Input
#parameter2:Output
library(ggplot2)
args=commandArgs(T)
InputFile=args[1]
OutputFile=args[2]
Colors=c("#BFBFBF","#00BFFF","#BFBFBF","#9400D3","#BFBFBF","#FF0000","#BFBFBF","#0000FF","#BFBFBF","#6495ED","#BFBFBF","#006400","#BFBFBF","#8B0000","#BFBFBF","#696969","#BFBFBF","#FF7F50","#BFBFBF","#A52A2A","#BFBFBF","#008000","#BFBFBF","#800080","#BFBFBF","#48D1CC","#BFBFBF","#DC143C","#BFBFBF","#7CFC00","#BFBFBF","#FFD700","#BFBFBF","#DA70D6","#BFBFBF","#DA503C","#BFBFBF","#FFFF00","#BFBFBF","#FFA500","#BFBFBF","#6B8E23","#BFBFBF","#FFEBCD","#BFBFBF","#DDA0DD","#BFBFBF","#00008B","#BFBFBF","#F5DEB3","#BFBFBF","#808080","#BFBFBF","#F0E68C","#BFBFBF","#87CEEB","#BFBFBF","#B0C4DE","#BFBFBF","#40E0D0","#BFBFBF","#CD5C5C","#BFBFBF","#FFC0CB","#BFBFBF","#000000")
data=read.table(InputFile,sep='\t',header=TRUE)
colnames(data)=c("RowName","Group","Value","ColorName")
png(OutputFile,width=2000,height=500)
print(ggplot(data,aes(x=RowName,y=Value,fill=ColorName))+geom_bar(stat="identity",width=0.5,position='dodge')+scale_fill_manual(values=Colors)+theme(legend.position="none")+theme(panel.background = element_rect(fill = "white"),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
dev.off()