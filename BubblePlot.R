#!/usr/bin/R
#Input: 
#parameter1:fold
#parameter2:mRNA Expression File, row: Gene mathched with RNA editing file, col: samples matched with RNA editing file without header
#parameter3:lncRNA Expression File, row: Gene mathched with RNA editing file, col: samples matched with RNA editing file without header
#parameter5:analysis results
#librarys
library(ggplot2)
library(ggrepel)
library(grid)
args=commandArgs(T)
#input files
InputFile=args[1]
DiscriValue=as.numeric(args[2])
ShowNum=as.numeric(args[3])
OutputFile=args[4]
FigureCol=as.numeric(args[5])
FigureRow=as.numeric(args[6])
cancers=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
data = read.table(InputFile,sep="\t",header=FALSE) 
colnames(data)=c("rowname","colname","value","size","cancer") 
png(OutputFile,width = FigureCol*500, height = FigureCol*300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(FigureRow,FigureCol)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}
num=0
data$color=ifelse(data$value>DiscriValue,"P","N")
for(cancer in cancers){
EachCancer=data[which(data$cancer==cancer),]
if(dim(EachCancer)[1]>0){
orders=order(EachCancer$size)
orders=orders[max((length(orders)-ShowNum+1),1):length(orders)]
FigureCancer=EachCancer[orders,]
FigureCancer$rowname=factor(FigureCancer$rowname,levels=FigureCancer$rowname)
if(length(unique(FigureCancer$color))==1){
if(data$color[1]=="N"){
color=c("#008000")
}else{
color=c("#800000")
}
}else{
color=c("#008000","#800000")
}
num=num+1
x_name=paste0("p",num)
assign(x_name,ggplot(FigureCancer,aes(colname,rowname,colour=color))+geom_point(aes(size=size))+scale_color_manual(values=color)+labs(size="size",title=cancer)+theme(legend.position = "none", plot.title = element_text(size=18,family = "Times",face = "bold"),axis.text=element_text(family = "Times",face = "bold",hjust=0.5,size=18),axis.title = element_text(size=0,face = "bold",vjust = 0.5,hjust = 0.5))+theme(panel.background = element_rect(fill = "white"),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
print(get(x_name), vp = vplayout(floor((num-1)/FigureCol)+1,((num-1)%%FigureCol+1)))
}
}
dev.off() 