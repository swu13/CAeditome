#!/usr/bin/R
#Input File1:Expression file, row: Gene, col:samples with header
#Input File2:clinical stage file, row: Tumor samples, col: RNA editing samples, RNA editing sample column, Clinical samples, Clinical stage or pathological stages roman, stage number without header
#Input File3:stage index in the clinical file
#Input Fold4:Correlation folder
#Inout File5:correlation results file
library(ggpubr)
args=commandArgs(T)
RawDataFile=args[1]
StageFile=args[2]
StageIndex=as.numeric(args[3])
OutputFolder=args[4]
OutputFile=args[5]
#Read files
RawData=read.table(RawDataFile,sep="\t",header=TRUE)
Stagedata=read.table(StageFile,sep="\t",header=FALSE)
Stage=as.numeric(as.vector(Stagedata[,StageIndex]))
Outputdata=RawData[,1:3]
colnames(Outputdata)=c("Gene","P","R")
for(i in 1:dim(RawData)[1]){
Eachdata=as.double(as.vector(RawData[i,2:dim(RawData)[2]]))
InformativeIndex=intersect(which(is.na(Eachdata)=="FALSE"),which(is.na(Stage)=="FALSE"))
if(length(InformativeIndex)==0 || all(Eachdata[InformativeIndex] == Eachdata[InformativeIndex[1]])){
Outputdata[i,2]=NA
Outputdata[i,3]=NA
}else{
correlation_results=cor.test(Eachdata,Stage,method="spearman")
Outputdata[i,2]=correlation_results$p.value
Outputdata[i,3]=correlation_results$estimate
}
if(!is.na(Outputdata[i,2]) && Outputdata[i,2]<0.05 && abs(Outputdata[i,3]) >0.3){
Figure_data=data.frame(Stage=Stage,Expression=Eachdata)
png(file=paste0(OutputFolder,"/",RawData[i,1],".png"),width=300,height=300)
print(ggplot(data = Figure_data, mapping = aes(x = Figure_data[,1],y = Figure_data[,2]))+ border()+ geom_point(colour = "black", size = 2,shape=1) + geom_smooth(method = lm,colour='red',fill='#E7E1D7') + stat_cor(method = "spearman")+ xlab("Tumor Stage") + theme(axis.text=element_text(face = "bold",hjust=0.5),axis.title.x = element_text(face = "bold",vjust = 0.5,hjust = 0.5))+ ylab("Expression")+ theme(axis.title.y = element_text(face = "bold",vjust = 0.5,hjust = 0.5))+theme(panel.background = element_rect(fill = "white"),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
#+annotate("text",x=2,y=0.9,label=paste0("P=",round(Outputdata[i,2],5),"; R=",round(Outputdata[i,3],5)))
dev.off()
}
}
write.table(Outputdata,OutputFile,quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)