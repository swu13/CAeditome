#!/usr/bin/R
#Input: 
#parameter1:gene list
#parameter2:png name
#librarys
library(enrichR)
library(ggplot2)
library(ggrepel)
library(grid)
args=commandArgs(T)
#input files
GeneListFile=args[1]
OutName=args[2]
OutFile=args[3]
EnrichIndex=as.numeric(args[4])
CutoffP=as.numeric(args[5])
CutoffQ=as.numeric(args[6])
FigureCol=as.numeric(args[7])
FigureRow=as.numeric(args[8])
data=read.table(GeneListFile,sep="\t",header=FALSE)
colnames(data)=c("Gene","Group")
Groups=unique(data$Group)
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018","KEGG_2021_Human")
colors=c("Maroon", "darkgreen", "navy","purple")
#enrichment results for each group
file_data=data.frame(Name="",p_value=0,q_value=0,score=0,type="",group="")
file_data=file_data[-1,]
for(group in Groups){
genes=unique(data$Gene[which(data$Group==group)])
enriched <- enrichr(genes, dbs[EnrichIndex])
enriched=as.data.frame(enriched)
colnames(enriched)=c("Term","Overlap","P.value","Adjusted.P.value","Old.P.value","Old.Adjusted.P.value","Odds.Ratio","Combined.Score","Genes")
if(dim(enriched)[1]>0){
file_data=rbind(file_data,data.frame(Name=paste0(enriched$Term,"(",enriched$Overlap,")"),p_value=enriched$P.value,q_value=enriched$Adjusted.P.value,score=enriched$Combined.Score,type=dbs[EnrichIndex],group=group))
}
}
#enrichment results for all 
genes=unique(data$Gene)
enriched <- enrichr(genes, dbs[EnrichIndex])
enriched=as.data.frame(enriched)
colnames(enriched)=c("Term","Overlap","P.value","Adjusted.P.value","Old.P.value","Old.Adjusted.P.value","Odds.Ratio","Combined.Score","Genes")
if(dim(enriched)[1]>0){
file_data=rbind(file_data,data.frame(Name=paste0(enriched$Term,"(",enriched$Overlap,")"),p_value=enriched$P.value,q_value=enriched$Adjusted.P.value,score=enriched$Combined.Score,type=dbs[EnrichIndex],group="All"))
}
write.table(file_data,OutFile,quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
#enrichment figure data
index=intersect(which(file_data$p_value<CutoffP),which(file_data$q_value<CutoffQ))
if(length(index)>0){
figure_data=file_data[index,]
figure_data=figure_data[order(figure_data$score),]
}
# enrichment figure for each group
png(paste0(OutName,".EachCancer.png"),width = FigureCol*600, height = FigureRow*300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(FigureRow,FigureCol)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}
num=0;
FigureGroups=sort(unique(figure_data$group))
FigureGroups=FigureGroups[which(FigureGroups!="All")]
FigureEach=figure_data[which(figure_data$group!="All"),]
for(i in 1:length(FigureGroups)){
num=num+1;
x_name=paste0("p",i)
FigureEachCancer=FigureEach[FigureEach$group==FigureGroups[i],]
FigureEachCancer$Name <- factor(FigureEachCancer$Name, levels=unique(FigureEachCancer$Name))
assign(x_name,ggplot(data = FigureEachCancer, mapping = aes(x = Name, y = score, fill = type)) + geom_bar(stat = 'identity',width = 0.5)+scale_fill_manual(values=colors[EnrichIndex],guide = FALSE)+coord_flip()+labs(title=FigureGroups[i],x="",y="")+theme(axis.text.x = element_blank(),axis.text.y=element_text(face = "bold",size=15),axis.ticks = element_blank()))
print(get(x_name), vp = vplayout(floor((num-1)/FigureCol)+1,((num-1)%%FigureCol)+1))
}
dev.off()
#enrichment figure for all
png(paste0(OutName,".All.png"),width = 400, height = 300)
FigureAll=figure_data[which(figure_data$group=="All"),]
FigureAll$Name <- factor(FigureAll$Name, levels=unique(FigureAll$Name))
print(ggplot(data = FigureAll, mapping = aes(x = Name, y = score, fill = type)) + geom_bar(stat = 'identity',width = 0.5)+scale_fill_manual(values=colors[EnrichIndex],guide = FALSE)+coord_flip()+labs(title="",x="",y="")+theme(axis.text.x = element_blank(),axis.text.y=element_text(face = "bold",size=15),axis.ticks = element_blank()))
dev.off()