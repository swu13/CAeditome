#!/usr/bin/R
library(igraph)
args=commandArgs(T)
NodeFile=args[1]
EdgeFile=args[2]
outpath=args[3]
nodes=read.table(NodeFile,sep="\t",header=FALSE)
colnames(nodes)=c("PngName","Group","node","color","size")
edges=read.table(EdgeFile,sep="\t",header=FALSE)
colnames(edges)=c("PngName","Group","node1","node2","color","size")
Unique_ENSG=unique(nodes$PngName)
for(i in 1:length(Unique_ENSG)){
EachNode=nodes[which(nodes$PngName==Unique_ENSG[i]),]
EachEdge=edges[which(edges$PngName==Unique_ENSG[i]),]
Groups=unique(EachNode$Group)
png(file=paste0(outpath,Unique_ENSG[i],".png"),width=500*min(3,length(Groups)),height=500*(floor((length(Groups)-1)/3)+1))
par(mfrow=c((floor((length(Groups)-1)/3)+1),min(3,length(Groups))))
for(j in 1:length(Groups)){
EachGroupNode=EachNode[which(EachNode$Group==Groups[j]),]
EachGroupEdge=EachEdge[which(EachEdge$Group==Groups[j]),]
graph <- graph_from_data_frame(EachGroupEdge[,3:4], directed = FALSE, vertices=EachGroupNode[,3])
set.seed(50)
l<-layout_as_star(graph)
V(graph)$label.cex=1
V(graph)$label.dist=0.5
V(graph)$size=EachGroupNode$size
V(graph)$color=EachGroupNode$color
E(graph)$color=EachGroupEdge$color
E(graph)$width=EachGroupEdge$size
print(plot(graph,main=Groups[j]))
}
dev.off()
}
