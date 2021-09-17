#!/usr/bin/R
#Input File:GTF file
#Output File: UTR3 file

library(GeneStructureTools)#conda install -c bioconda bioconductor-genestructuretools

args=commandArgs(T)
gtfFile=args[1]
UTR3File=args[2]
gtf <- rtracklayer::import(gtfFile)
gtf <- as.data.frame(UTR2UTR53(gtf))
UTR3index=which(gtf$type=="UTR3")
write.table(gtf[UTR3index,],UTR3File,quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)