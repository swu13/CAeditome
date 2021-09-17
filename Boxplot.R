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
Flag=args[5]
ymin=as.numeric(args[6])
ymax=as.numeric(args[7])
figure_data = read.table(InputFile,sep="\t",header=FALSE) 
colnames(figure_data)=c("Type","rowname","Value")
figure_data$Type <- factor(figure_data$Type, levels=unique(figure_data$Type))
png(file=OutputFile,width=PngWidth,height=PngHeight)
print(ggboxplot(figure_data, x = "Type", y = "Value",color = "Type", palette = "rainbow12",add = "none")+border()+ylim(ymin,ymax)+theme(legend.position="none",plot.title=element_text(size=0,hjust=0.5,face = "bold"),axis.text=element_text(face = "bold",hjust=0.5),axis.title.x=element_text(size=0),axis.title.y=element_text(face = "bold"))+stat_compare_means(aes(label = ..p.signif..),method = "t.test", ref.group = "0")+rotate_x_text(angle = 90))
dev.off()
if(Flag=="TRUE"){
types=unique(figure_data$Type)
Output=rep(0,length(types))
for(i in 1:length(types)){
ty=types[i]
Output[i]=length(unique(figure_data$rowname[which(figure_data$Type==ty)]))
} 
print(types)
print(Output)
}