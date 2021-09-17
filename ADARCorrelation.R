#!/usr/bin/R
#Input: 
#parameter1:RNA editing file, row: RNA editing site, col:Samples without header
#parameter2:Expression File, row: Gene mathched with RNA editing file, col: samples matched with RNA editing file without header
#parameter3:boxplot folder
#parameter4:correlation plot folder
#parameter5:analysis results
#librarys
library(ggpubr)
args=commandArgs(T)
#input files
RawRNAeditingFile=args[1]
RawExpressionFile=args[2]
BoxplotFold=args[3]
CorrelationPlotFold=args[4]
OutputFile=args[5]
#read files
RawRNAeditingData=read.table(RawRNAeditingFile,sep="\t",header=FALSE)
RawExpressionData=read.table(RawExpressionFile,sep="\t",header=FALSE)
Outputdata=RawRNAeditingData[,c(1,1:7)]
Outputdata[,2]=RawExpressionData[,1]
colnames(Outputdata)=c("EditingInformation","Gene","TtestP","MeanExpressionEdited","MeanExpressionNonedited","logFC","CorreP","CorreR")
for(i in 1:dim(RawRNAeditingData)[1]){
EachEditing=as.double(as.vector(RawRNAeditingData[i,2:dim(RawRNAeditingData)[2]]))
EachExpression=as.double(as.vector(RawExpressionData[i,2:dim(RawExpressionData)[2]]))
#Ttest
EditedIndex=intersect(which(is.na(EachEditing)=="FALSE"),which(is.na(EachExpression)=="FALSE"))
NonEditedIndex=intersect(which(is.na(EachEditing)=="TRUE"),which(is.na(EachExpression)=="FALSE"))
if(length(EditedIndex)<=2 || length(NonEditedIndex)<=2 || all(EachExpression[EditedIndex]==EachExpression[EditedIndex][1]) || all(EachExpression[NonEditedIndex]==EachExpression[NonEditedIndex][1])){
Outputdata[i,3]=NaN
}else{
Outputdata[i,3]=t.test(EachExpression[EditedIndex],EachExpression[NonEditedIndex])$p.value
}
Outputdata[i,4]=mean(EachExpression[EditedIndex])
Outputdata[i,5]=mean(EachExpression[NonEditedIndex])
Outputdata[i,6]=log2(Outputdata[i,4]/Outputdata[i,5])
#Correlation
if(length(EditedIndex)<=2 || all(EachExpression[EditedIndex] == EachExpression[EditedIndex[1]])){
Outputdata[i,7:8]=c(NaN,NaN)
}else{
correlation_results=cor.test(EachEditing[EditedIndex],EachExpression[EditedIndex],method="pearson")
Outputdata[i,7]=correlation_results$p.value
Outputdata[i,8]=correlation_results$estimate
}
#Boxplot
if(!is.na(Outputdata[i,3]) && Outputdata[i,3]<0.05 && !is.na(Outputdata[i,6]) && abs(Outputdata[i,6])>2){
figure_data=data.frame(Group=c(rep("Edited",length(EditedIndex)),rep("NonEdited",length(NonEditedIndex))),Expression=c(EachExpression[EditedIndex],EachExpression[NonEditedIndex]))
png(file=paste0(BoxplotFold,"/",RawRNAeditingData[i,1],".",RawExpressionData[i,1],".png"),width=300,height=300)
print(ggboxplot(figure_data, x = "Group", y = "Expression",color = "Group", palette = c("#800000","#000080"),add = "jitter")+border()+ylab("Expression")+theme(legend.position="none",axis.text=element_text(face = "bold",hjust=0.5),axis.title.x=element_text(size=0),axis.title.y=element_text(face = "bold"))+stat_compare_means(aes(group=Group),label = "p.signif",label.y=max(figure_data$Expression)*0.9,label.x=1.5,face="bold",hide.ns=T, paired=F,method="t.test"))
dev.off()
}
#Correlation Plot
if(!is.na(Outputdata[i,7]) && Outputdata[i,7]<0.05 && !is.na(Outputdata[i,8]) && abs(Outputdata[i,8])>0.3){
Figure_data=data.frame(Editing=EachEditing[EditedIndex],Expression=EachExpression[EditedIndex])
png(file=paste0(CorrelationPlotFold,"/",RawRNAeditingData[i,1],".",RawExpressionData[i,1],".png"),width=300,height=300)
print(ggplot(data = Figure_data, mapping = aes(x = Figure_data[,2],y = Figure_data[,1])) +border()+ geom_point(colour = "black", size = 2,shape=1) + geom_smooth(method = lm,colour='red',fill='#E7E1D7') + stat_cor(method = "pearson")+ xlab("Expression") + theme(axis.text=element_text(face = "bold",hjust=0.5),axis.title.x = element_text(face = "bold",vjust = 0.5,hjust = 0.5))+ ylab("Editing Frequency")+ theme(axis.title.y = element_text(face = "bold",vjust = 0.5,hjust = 0.5))+theme(panel.background = element_rect(fill = "white"),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
dev.off()
}
}
write.table(Outputdata,OutputFile,quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)