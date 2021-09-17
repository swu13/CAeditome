#!/usr/bin/R
#Input: 
#parameter1: Raw data file without header
#parameter2: roman column
#parameter3: OutputFile for transfered results
#Input parameter
args=commandArgs(T)
RawDataName=args[1]
colname=as.numeric(args[2])
OutputFile=args[3]
#read the files
RawData=read.table(RawDataName,sep="\t",header=FALSE)
for(i in 1:dim(RawData)[1]){
if(i==1){
Output=as.numeric(as.roman(as.character(RawData[i,colname])))
}else{
Output=c(Output,as.numeric(as.roman(as.character(RawData[i,colname]))))
}
}
RawData$roman2num=Output
write.table(RawData,OutputFile,quote = FALSE,sep="\t",row.names=FALSE,col.names=FALSE) 