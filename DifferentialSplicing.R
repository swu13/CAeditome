#!/usr/bin/R
#Input: 
#parameter1: splicing file,Column:Samples; Row: RNA editing without header
#parameter2: group file without header
#parameter3: OutputFile for Test results
#parameter4: BoxplotOutput with the sequence started with the keyword
#
#Input parameter
library(ggpubr)
args=commandArgs(T)
RawDataName=args[1]
GroupName=args[2]
OutputFile=args[3]
BoxplotFold=args[4]
#read the files
RawData=read.table(RawDataName,sep="\t",header=TRUE)
Group=read.table(GroupName,sep="\t",header=FALSE)
OutputData=RawData[,1:7]
colnames(OutputData)=c("Gene","P","MeanTumor","MeanNormal","logFC","TumorSamples","NormalSamples")
#Test
for(i in 1:dim(RawData)[1]){
EachData=as.numeric(RawData[i,2:dim(RawData)[2]])
TumorIndex=which(Group$V1=="Tumor")
NormalIndex=which(Group$V1=="Normal")
OutputData[i,6]=length(TumorIndex)
OutputData[i,7]=length(NormalIndex)
TumorNonNAindex=intersect(TumorIndex,which(is.na(EachData)=="FALSE"))
NormalNonNAindex=intersect(NormalIndex,which(is.na(EachData)=="FALSE"))
if(length(TumorNonNAindex)==0){
OutputData[i,3]=NA
}else{
OutputData[i,3]=mean(EachData[TumorNonNAindex])
}
if(length(NormalNonNAindex)==0){
OutputData[i,4]=NA
}else{
OutputData[i,4]=mean(EachData[NormalNonNAindex])
}
if(length(TumorNonNAindex)<=2 || length(NormalNonNAindex)<=2 || all(EachData[TumorNonNAindex]==EachData[TumorNonNAindex][1]) || all(EachData[NormalNonNAindex]==EachData[NormalNonNAindex][1])){
OutputData[i,2]=NA
}else{
OutputData[i,2]=t.test(EachData[TumorNonNAindex],EachData[NormalNonNAindex])$p.value
}
OutputData[i,5]=log2(OutputData[i,3]/OutputData[i,4])
if(!is.na(OutputData[i,2]) && OutputData[i,2]<0.05){
status=rep("Tumor",length(EachData))
status[NormalIndex]="Normal"
figure_data=data.frame(Group=status,Editing=EachData)
png(file=paste0(BoxplotFold,"/",RawData[i,1],".png"),width=300,height=300)
print(ggboxplot(figure_data, x = "Group", y = "Editing",color = "Group", palette = c("#800000", "#000080"),add = "jitter")+border()+ylab("PSI")+theme(legend.position="none",axis.text=element_text(face = "bold",hjust=0.5),axis.title.x=element_text(size=0),axis.title.y=element_text(face = "bold"))+stat_compare_means(aes(group=Group),label = "p.signif",label.y=0.9,label.x=1.5,face="bold",hide.ns=T, paired=F,method="t.test"))
dev.off()
}
}
write.table(OutputData,OutputFile,quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE) 