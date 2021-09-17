#!/usr/bin/R
####################################################
#draw the heatmap for any cancer-related files
#parameter1: Input File containing the heatmap_row and heatmap_col and values in each heatmap: Also containing the gene file
#parameter2: The col for the cancer type in the file
#parameter3: The col for the genes in the file
#parameter4: The col for the editing ID in the file
#parameter5:The col for the values in the file
#parameter6: the heatmap folder
####################################################
library(pheatmap) #install.packages("pheatmap")
args=commandArgs(T)
Filename=args[1] #/data3/swu/CAeditome/Differential/Tmp/Heatmap.data
CancerCol=as.numeric(args[2])
GeneCol=as.numeric(args[3])
EditingCol=as.numeric(args[4])
ValueCol=as.numeric(args[5])
DifferenatialNumber=as.numeric(args[6])
UpMax=as.numeric(args[7])
FigureFolder=args[8]
RawData=read.table(Filename,sep="\t",header=FALSE)
Cancers=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
Unique_ENSG=unique(RawData[,GeneCol])
for(i in 1:length(Unique_ENSG)){
index=which(RawData[,GeneCol]==Unique_ENSG[i])
Data_gene=RawData[index,]
Unique_Editing=unique(Data_gene[,EditingCol])
FigureData=c(rep(NA,length(Cancers)))
for(j in 1:length(Unique_Editing)){
index=which(Data_gene[,EditingCol]==Unique_Editing[j])
EachEditing=Data_gene[index,]
for(k in 1:length(Cancers)){
index=which(EachEditing[,CancerCol]==Cancers[k])
if(length(index)==0 && k==1){
FigureEach=NA
}else if(length(index)==0){
FigureEach=c(FigureEach,NA)
}else if(k==1){
FigureEach=as.numeric(EachEditing[index,ValueCol])
}else{
FigureEach=c(FigureEach,as.numeric(EachEditing[index,ValueCol]))
}
}
FigureEach[which(FigureEach>UpMax)]=UpMax
FigureData=rbind(FigureData,FigureEach)
}
FigureData=as.matrix(FigureData[-1,])
if(dim(FigureData)[1]==33){
FigureData=t(as.data.frame(FigureData))
}else{
FigureData=as.data.frame(FigureData)
}
colnames(FigureData)=Cancers
rownames(FigureData)=Unique_Editing
unique_value=unique(as.numeric(as.matrix(FigureData)))
unique_value=unique_value[which(is.na(unique_value)==FALSE)]
if(length(unique_value)>0){
figurename=paste0(FigureFolder,"/",Unique_ENSG[i],".png")
if(min(unique_value)<DifferenatialNumber){
bk1=seq(min(unique_value),max(DifferenatialNumber-0.001,min(unique_value)),by=0.001)
}else{
bk1=vector()
}
if(max(unique_value)>DifferenatialNumber){
bk2=seq(DifferenatialNumber,max(unique_value),by=0.001)
}else{
bk2=vector()
}
color=c(colorRampPalette(colors = c("darkgreen","white"))(length(bk1)),colorRampPalette(colors = c("white","firebrick3"))(length(bk2)))
if(length(color)>1){
pheatmap(FigureData,cluster_row = FALSE,cluster_cols=FALSE,color= color,cellwidth = 5,cellheight = 3,fontsize=3,file=figurename,na_col = "grey50",breaks=c(bk1,bk2))
}
}
}