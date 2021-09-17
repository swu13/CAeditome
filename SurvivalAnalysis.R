#!/usr/bin/R
#parameter1: RNA editing file, row: RNA editing site, col:tumor samples with header
#parameter2: clinical file, row: Tumor samples, col: RNA editing samples, RNA editing sample column, Clinical samples, survival time and survival status without header
#parameter3: Index of survival time in the clinical file
#parameter4: Index of survival status in the clinical file
#parameter5: survival plots folder
#parameter6: survival results file
#librarys
library("survival")
library("survminer")
#parameters
args=commandArgs(T)
RawFile=args[1]
ClinicalFile=args[2]
SurvivalTimeIndex=as.numeric(args[3])
SurvivalStatusIndex=as.numeric(args[4])
OutputFold=args[5]
OutputFile=args[6]
#Read files
RawData=read.table(RawFile,sep="\t",header=TRUE)
SurvivalData=read.table(ClinicalFile,sep="\t",header=FALSE)
SurvialTime=as.numeric(as.vector(SurvivalData[,SurvivalTimeIndex]))
SurvivalStatusCha=as.character(as.vector(SurvivalData[,SurvivalStatusIndex]))
SurvivalStatusNum=ifelse(SurvivalStatusCha == "alive",0,ifelse(SurvivalStatusCha == "dead",1,NA))
Outputdata=RawData[,1:10]
colnames(Outputdata)=c("EditingInformation","Pkm","Pcox_discrete","HR_discrete","HR_lower_discrete","HR_upper_discrete","Pcox_continuous","HR_continuous","HR_lower_continuous","HR_upper_continuous")
#survival analysis
for(i in 1:dim(RawData)[1]){
EachData=as.double(as.vector(RawData[i,2:dim(RawData)[2]]))
InformativeIndex=setdiff(1:length(EachData),union(union(which(is.na(EachData)=="TRUE"),which(is.na(SurvivalStatusNum)=="TRUE")),union(intersect(which(SurvivalStatusNum==1),which(is.na(SurvialTime)=="TRUE")),intersect(which(SurvivalStatusNum==0),which(is.na(SurvialTime)=="FALSE")))))
if(length(InformativeIndex)==0 || all(EachData[InformativeIndex] == EachData[InformativeIndex[1]])){
Outputdata[i,2:10]=c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
}else{
InformativeEditing=EachData[InformativeIndex]
InformativeStatus=SurvivalStatusNum[InformativeIndex]
InformativeTime=SurvialTime[InformativeIndex]
TranferEditing=ifelse(InformativeEditing > median(InformativeEditing),1,0)
survival_data <- data.frame(group = TranferEditing,status = InformativeStatus,time = InformativeTime,stringsAsFactors = F)
#KM
HighIndex=which(TranferEditing==1)
LowIndex=which(TranferEditing==0)
HighIndexDead=which(unique(survival_data[HighIndex,])$status==1)
LowIndexDead=which(unique(survival_data[LowIndex,])$status==1)
if(dim(unique(survival_data))[1]<2 || dim(unique(survival_data[HighIndex,]))[1]<2 || dim(unique(survival_data[LowIndex,]))[1]<2 || (length(HighIndexDead)<2 && length(LowIndexDead)<2 && unique(survival_data[HighIndex,])$time[HighIndexDead]==unique(survival_data[LowIndex,])$time[LowIndexDead])){
Outputdata[i,2]=NA
}else{
kmfit<-survfit(Surv(time, status) ~ group, data =  survival_data)
diff=survdiff(Surv(time, status)~group,data=survival_data)
Outputdata[i,2]=1-pchisq(diff$chisq,df=length(diff$n)-1)#pvalue
}
#COX--high/low
if(length(which(is.na(survival_data$time)=="FALSE"))<2){
Outputdata[i,3:6]=c(NA,NA,NA,NA)
}else{
coxfit <- coxph(Surv(time, status) ~ group,data =  survival_data)
beta=coef(coxfit)
se=sqrt(diag(vcov(coxfit)))
Outputdata[i,3]=1 - pchisq((beta/se)^2, 1) #p.value=1 - pchisq((beta/se)^2, 1),beta=coef(coxfit),se=sqrt(diag(vcov(coxfit)))
Outputdata[i,4]=exp(beta)  #HR
Outputdata[i,5] = exp(beta - qnorm(.975, 0, 1) * se) # HR lower 95%
Outputdata[i,6] = exp(beta + qnorm(.975, 0, 1) * se) #HR upper 95%
}
#COX--continue
survival_data <- data.frame(group = InformativeEditing,status = InformativeStatus,time = InformativeTime,stringsAsFactors = F)
if(length(which(is.na(survival_data$time)=="FALSE"))<2){
Outputdata[i,7:10]=c(NA,NA,NA,NA)
}else{
coxfit <- coxph(Surv(time, status) ~ group,data =  survival_data)
summary_Results=summary(coxfit)
beta=signif(summary_Results$coef[1])
Outputdata[i,7]=signif(summary_Results$wald["pvalue"])
Outputdata[i,8]=signif(summary_Results$coef[2])  #HR
Outputdata[i,9] = signif(summary_Results$conf.int[,"lower .95"], 2) # HR lower 95%
Outputdata[i,10] = signif(summary_Results$conf.int[,"upper .95"],2) #HR upper 95%
}
}
if((!is.na(Outputdata[i,2])) && Outputdata[i,2]<0.05){
png(file=paste0(OutputFold,"/",RawData[i,1],".png"),width=300,height=300)
pval_value=ifelse(Outputdata[i,2] < 0.001, "Pkm<0.001", paste0("Pkm=",round(Outputdata[i,2],3)))
if(Outputdata[i,4]>1){
lengend_name=c("LowEditing","HighEditing")
color_name=c("blue4","red3")
}else{
lengend_name=c("LowEditing","HighEditing")
color_name=c("red3", "blue4")
}
p=ggsurvplot(kmfit,pval = FALSE,conf.int = FALSE,risk.table = FALSE,linetype = "strata",ggtheme = theme_bw(),palette = color_name,xlab = "Follow up time",legend = c(0.8,0.75), legend.title = "",legend.labs = lengend_name,font.x = c(10, "bold", "black"),font.y = c(10, "bold", "black"),font.tickslab = c(10, "bold", "black"),font.legend=c(10, "bold", "black"))
p$plot <- p$plot +annotate("text", x = 0, y = 0, label = pval_value, cex=3.8, vjust="bottom", hjust = "left", fontface=2)
print(p)
dev.off()
}
}
write.table(Outputdata,OutputFile,quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)