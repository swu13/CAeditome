#!/usr/bin/R
#Input: 
#parameter1:Input
#parameter2:Output
library(ggplot2)
args=commandArgs(T)
InputFile=args[1]
OutputFile=args[2]
cut_off_pvalue=as.numeric(args[3])
cut_off_logFC=as.numeric(args[4])
FigureCol=as.numeric(args[5])
FigureRow=as.numeric(args[6])
ymax=as.numeric(args[7])
data=read.table(InputFile,sep="\t",header=FALSE)
colnames(data)=c("Group","DotName","Column","Row")
data$change = ifelse(data$Column < cut_off_pvalue & abs(data$Row) >= cut_off_logFC,ifelse(data$Row> cut_off_logFC ,'Up','Down'),'Stable')
png(OutputFile,width = FigureCol*500, height = FigureRow*500)
print(ggplot(data,aes(x = Row, y=-log10(Column),colour=change))+geom_point(alpha=0.4, size=3) +scale_color_manual(values=c("#008000","#000000","#800000"))+geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8)+geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +ylim(0,ymax)+ theme(legend.position="none",plot.title = element_text(size=30,family = "Times",face = "bold"),axis.text=element_text(family = "Times",face = "bold",hjust=0.5,size=30),axis.title = element_text(size=0),strip.text = element_text(family = "Times",face = "bold",size = 25))+facet_wrap(~Group, ncol = FigureCol,as.table=TRUE))
dev.off()