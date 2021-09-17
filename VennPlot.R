#!/usr/bin/R
#Input: 
#parameter1:Input
#parameter2:Output
#
library(VennDiagram)# install.packages("VennDiagram")
args=commandArgs(T)
InputFile=args[1]
OutputFile=args[2]
Flag=args[3]
data=read.table(InputFile,sep="\t",header=FALSE)
colnames(data)=c("Editing","Gene","group")
Groups=unique(data$group)
for(i in 1:length(Groups)){
group=Groups[i]
index=which(data$group==group)
if(i==1){
All1=list(group = unique(data$Editing[index]))
}else{
All1=c(All1,list(group=unique(data$Editing[index])))
}
}
venn.diagram(All1,filename=OutputFile,lwd=1,lty=2,height = 300, width = 300, imagetype="png", category.names=Groups,col ="transparent", fill=c('#FBE5D6','#DEEBF7','#E2F0D9'),cat.cex=0.15,cat.fontfamily="Times",cat.fontface=1,alpha = 0.5,cex=0.15,fontfamily="Times",fontface=1,reverse=TRUE)
#calculate genes number for each group
Output=rep(0,7)
if(Flag=="TRUE"){
Index=list(A=data$Editing[which(data$group==Groups[1])],B=data$Editing[which(data$group==Groups[2])],C=data$Editing[which(data$group==Groups[3])])
editing=intersect(intersect(Index$A,Index$B),Index$C)
Output[1]=length(unique(ifelse(data$Editing %in% editing,data$Gene,"")))-1
editing=setdiff(intersect(Index$A,Index$B),intersect(intersect(Index$A,Index$B),Index$C))
Output[2]=length(unique(ifelse(data$Editing %in% editing,data$Gene,"")))-1
editing=setdiff(intersect(Index$C,Index$A),intersect(intersect(Index$A,Index$B),Index$C))
Output[3]=length(unique(ifelse(data$Editing %in% editing,data$Gene,"")))-1
editing=setdiff(intersect(Index$C,Index$B),intersect(intersect(Index$A,Index$B),Index$C))
Output[4]=length(unique(ifelse(data$Editing %in% editing,data$Gene,"")))-1
editing=setdiff(setdiff(setdiff(Index$A,intersect(Index$A,Index$B)),intersect(Index$A,Index$C)),intersect(intersect(Index$A,Index$B),Index$C))
Output[5]=length(unique(ifelse(data$Editing %in% editing,data$Gene,"")))-1
editing=setdiff(setdiff(setdiff(Index$B,intersect(Index$A,Index$B)),intersect(Index$B,Index$C)),intersect(intersect(Index$A,Index$B),Index$C))
Output[6]=length(unique(ifelse(data$Editing %in% editing,data$Gene,"")))-1
editing=setdiff(setdiff(setdiff(Index$C,intersect(Index$C,Index$B)),intersect(Index$A,Index$C)),intersect(intersect(Index$A,Index$B),Index$C))
Output[7]=length(unique(ifelse(data$Editing %in% editing,data$Gene,"")))-1
print(Output)
print(Groups)
}