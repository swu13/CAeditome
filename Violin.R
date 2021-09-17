#!/usr/bin/R
#Input: 
#parameter1:InputFile
#parameter2: OutputFile
#librarys
library(ggpubr)
args=commandArgs(T)
#input files
InputFile=args[1]
OutputFile=args[2]
PngWidth=as.numeric(args[3])
PngHeight=as.numeric(args[4])
figure_data = read.table(InputFile,sep="\t",header=FALSE) 
colnames(figure_data)=c("Type","rowname","Value")
figure_data$Type <- factor(figure_data$Type, levels=unique(figure_data$Type))
png(file=OutputFile,width=PngWidth,height=PngHeight)
print(ggplot(figure_data, aes(x =Type, y =Value))+geom_violin(fill="navy")+border()+theme(legend.position="none",plot.title=element_text(size=0,hjust=0.5,face = "bold"),axis.text=element_text(face = "bold",hjust=0.5),axis.title.x=element_text(size=0),axis.title.y=element_text(face = "bold"))+theme_bw()+theme(panel.grid=element_blank())+stat_compare_means(aes(label = ..p.signif..),method = "t.test", ref.group = "0"))
dev.off()