#!/usr/bin/R
#Input: CancerType, outputname
library(TCGAbiolinks)
library(plyr)
library(limma)
library(biomaRt)
library(SummarizedExperiment)
args=commandArgs(T)
cancer_type = args[1]
outputname = args[2]
Type=args[3]

if(Type == "Gene"){
#gene cancer_type = 'TCGA-BRCA'
query <- GDCquery(project = paste0("TCGA-",cancer_type),data.category = "Transcriptome Profiling",data.type = "Gene Expression Quantification",workflow.type = "HTSeq - FPKM")
GDCdownload(query = query)
dataPrep1 <- GDCprepare(query = query, save = TRUE, save.filename = paste0(outputname,cancer_type,"_case1.rda"))
expdat<-assay(dataPrep1)
write.csv(expdat,file = paste0(outputname,cancer_type,"_expression.csv"),quote = FALSE)
dataPrep2 <- TCGAanalyze_Preprocessing(object = dataPrep1, cor.cut = 0.3, datatype = "HTSeq - FPKM")
write.csv(dataPrep2,file = paste0(outputname,cancer_type,"_expression_preprocessing.csv"),quote = FALSE)
}

#if(Type == "isoform"){
#query <- GDCquery(project = paste0("TCGA-",cancer_type),data.category = "Transcriptome Profiling",data.type = "Isoform Expression Quantification",workflow.type = "HTSeq - FPKM")
#GDCdownload(query = query)
#dataPrep1 <- GDCprepare(query = query, save = TRUE, save.filename = paste0(outputname,cancer_type,"_case1.rda"))
#expdat<-assay(dataPrep1)
#write.csv(expdat,file = paste0(outputname,cancer_type,"_expression.csv"),quote = FALSE)
#}




