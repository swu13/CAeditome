#!/usr/bin/R
#Input: 
#parameter1:RNA editing file, row: RNA editing site, col:Samples without header
#parameter2:Expression File, row: Gene mathched with RNA editing file, col: samples matched with RNA editing file without header
#parameter3:analysis results
#librarys
args=commandArgs(T)
#input files
RawRNAeditingFile=args[1]
RawExpressionFile=args[2]
OutputFile=args[3]
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
}
write.table(Outputdata,OutputFile,quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)