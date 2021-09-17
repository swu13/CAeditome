#For the manuscript
cancertypes=("BLCA"  "LGG"  "BRCA"  "LIHC"  "CESC"  "LUAD"  "COAD"  "LUSC"  "GBM"  "PRAD"  "HNSC"  "STAD"  "KICH"  "THCA"  "KIRC"  "UCEC"  "KIRP"  "ACC"  "DLBC"  "ESCA"  "LAML"  "OV"  "PAAD"  "PCPG"  "SKCM"  "UCS"  "UVM"  "CHOL"  "MESO"  "READ"  "SARC"  "TGCT"  "THYM")
StatisticalFold="/data9/swu13/CAeditome_hg38/Statistic/"

#0. Editing events
function EditingEventsAnalysis(){
mkdir $StatisticalFold"EditingEvents/"
rm -f $StatisticalFold"EditingEvents/EditingEvents"
cut -f1 /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation | tail -n +2 | sort | uniq | wc -l > $StatisticalFold"EditingEvents/EditingEvents" # 1530875
cut -f3 /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation | tail -n +2 | sort | uniq | wc -l >> $StatisticalFold"EditingEvents/EditingEvents" # 18414
awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
            ARGIND==2 && (FNR>1) && ($5 in Tumor){print $5}'\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation | sort | uniq | wc -l
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
awk -F"\t" -v cancer="$cancer" '{if(NR==1){num=0}else{num=num+$5+$6}}END{print cancer"\t"num}'\
            "/data9/swu13/CAeditome_hg38/Differential/Files/"$cancer".DifferentialResults" >> $StatisticalFold"EditingEvents/EditingEvents"
done
}

#0.1 gene for each cancer
function GeneCancer(){
echo -e "Cancer\tENSG\tGeneType\tGeneName" > $StatisticalFold"EditingEvents/GeneCancer"
#for((i=0;i<=32;i++))
#do
#cancer=${cancertypes[$i]}
#awk -F"\t" -v cancer="$cancer" '{print cancer"\t"$2"\t"$3"\t"$4}'\
#            "/data9/swu13/CAeditome_hg38/GeneralAnnotation/"$cancer".Table" | sort | uniq >> $StatisticalFold"EditingEvents/GeneCancer"
#done
awk -F"\t" 'ARGIND==1 && (FNR>1){split($3,a,".");GeneType[a[1]]=a[1]"\t"$4"\t"$5}\
            ARGIND==2 && (($5=="NA" && $8>=50) || ($5<0.0001)){split($3,a,".");print $12"\t"GeneType[a[1]]}'\
            ../GeneralAnnotation/GeneralInformation\
            $StatisticalFold"Specific/TumorSpecific" | sort | uniq > $StatisticalFold"EditingEvents/GeneCancer.1"
awk -F"\t" 'ARGIND==1 && (FNR>1){split($3,a,".");GeneType[a[1]]=a[1]"\t"$4"\t"$5}\
            ARGIND==2 && (FNR>1) && ($5<0.001){split($3,a,".");print $7"\t"GeneType[a[1]]}'\
            ../GeneralAnnotation/GeneralInformation\
            $StatisticalFold"Stage/StageAssociation" | sort | uniq >> $StatisticalFold"EditingEvents/GeneCancer.1"
awk -F"\t" 'ARGIND==1 && (FNR>1){split($3,a,".");GeneType[a[1]]=a[1]"\t"$4"\t"$5}\
            ARGIND==2 && (FNR>1) && ($5<0.01) && ($6<0.01){split($3,a,".");print $8"\t"GeneType[a[1]]}'\
            ../GeneralAnnotation/GeneralInformation\
            $StatisticalFold"Survival/SurvivalAssociation" | sort | uniq >> $StatisticalFold"EditingEvents/GeneCancer.1"            
awk -F"\t" 'ARGIND==1 && (FNR>1){split($3,a,".");GeneType[a[1]]=a[1]"\t"$4"\t"$5}\
            ARGIND==2 && ($6*$8>0) && ($5<0.01) && ($7<0.01){split($3,a,".");print $9"\t"GeneType[a[1]]}'\
            ../GeneralAnnotation/GeneralInformation\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq >> $StatisticalFold"EditingEvents/GeneCancer.1"
awk -F"\t" 'ARGIND==1 && (FNR>1){split($4,a,".");ENSG[$6]=a[1]"\t"$5"\t"$6}\
            ARGIND==2{print $1"\t"ENSG[$2]}'\
            /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinLolliplot\
            $StatisticalFold"Protein/SequenceCancer" | sort | uniq >> $StatisticalFold"EditingEvents/GeneCancer.1"
awk -F"\t" 'ARGIND==1 && (FNR>1){split($3,a,".");GeneType[a[1]]=a[1]"\t"$4"\t"$5}\
            ARGIND==2 && ($9*$11>0){split($4,a,".");print $1"\t"GeneType[a[1]]}'\
            ../GeneralAnnotation/GeneralInformation\
            $StatisticalFold"SplicingAnalysis/PSI" | sort | uniq >> $StatisticalFold"EditingEvents/GeneCancer.1"
awk -F"\t" 'ARGIND==1 && (FNR>1){split($3,a,".");GeneType[a[1]]=a[1]"\t"$4"\t"$5;if($4=="miRNA"){MiRNA[$2]=a[1]"\t"$4"\t"$5}}\
            ARGIND==2 && ($25=="Y" || $26=="Y"){if($6~"EditedmiRNA_"){if($3 in MiRNA){print $1"\t"MiRNA[$3]}}\
                      else{split($4,a,".");if(a[1] in GeneType){print $1"\t"GeneType[a[1]]}}}'\
            ../GeneralAnnotation/GeneralInformation\
            $StatisticalFold"miRNA/miRNARegulation" | sort | uniq >> $StatisticalFold"EditingEvents/GeneCancer.1"
echo -e "Cancer\tENSG\tGeneType\tGeneName" > $StatisticalFold"EditingEvents/GeneCancer"
cat $StatisticalFold"EditingEvents/GeneCancer.1" |sort | uniq >> $StatisticalFold"EditingEvents/GeneCancer"
}

#0.2 gene for each analysis
function GeneEachAnalysis(){
echo -e "AnalysisType\tENSG\tGeneName" > $StatisticalFold"EditingEvents/GeneEachAnalysis"
awk -F"\t" 'ARGIND==1 && (FNR>1){GeneName[$3]=$5}\
            ARGIND==2 && (($5=="NA" && $8>=100) || ($5<0.00001)){split($3,a,".");\
                          print "Tumor Specific RNA editing events\t"a[1]"\t"GeneName[$3]}'\
            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
            $StatisticalFold"Specific/TumorSpecific" | sort | uniq >> $StatisticalFold"EditingEvents/GeneEachAnalysis"
awk -F"\t" 'ARGIND==1 && (FNR>1){GeneName[$3]=$5}\
            ARGIND==2 && ($5<0.001) && (FNR>1){split($3,a,".");\
                          print "Correlation between tumor stages and RNA editing frequency\t"a[1]"\t"GeneName[$3]}'\
            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
            $StatisticalFold"Stage/StageAssociation" | sort | uniq >> $StatisticalFold"EditingEvents/GeneEachAnalysis"
awk -F"\t" 'ARGIND==1 && (FNR>1){GeneName[$3]=$5}\
            ARGIND==2 && ($5<0.05) && ($6<0.05) && (FNR>1){split($3,a,".");\
                          print "Correlation between tumor survival and RNA editing frequency\t"a[1]"\t"GeneName[$3]}'\
            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
            $StatisticalFold"Survival/SurvivalAssociation" | sort | uniq >> $StatisticalFold"EditingEvents/GeneEachAnalysis"
awk -F"\t" '($5<0.001) && ($6*$8>0){print "Correlation between gene expression and RNA editing frequency\t"$3"\t"$4}'\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq >> $StatisticalFold"EditingEvents/GeneEachAnalysis"
awk -F"\t" '(NR>1) && ($9=="D" || $10=="D" || $11=="D" || $12=="D"){split($4,a,".");print "Protein coding RNA editing(s)\t"a[1]"\t"$5}'\
            /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinFunction | sort | uniq >> $StatisticalFold"EditingEvents/GeneEachAnalysis"
awk -F"\t" '{split($4,a,".");print "RNA A-to-I editing(s) in the alternative splicing sites\t"a[1]"\t"$5}'\
            $StatisticalFold"SplicingAnalysis/PSI" | sort | uniq >> $StatisticalFold"EditingEvents/GeneEachAnalysis"
            
awk -F"\t" 'ARGIND==1 && (FNR>1){split($3,a,".");if($4=="miRNA"){MiRNA[$2]=a[1]"\t"$5}}\
            ARGIND==2{split($4,a,".");split($9,b,".");if($6~"EditedmiRNA_"){if($3 in MiRNA){\
                           print "The influence of RNA editing events on miRNA binding and regulation\t"MiRNA[$3]}}\
                      else{split($4,a,".");print "The influence of RNA editing events on miRNA binding and regulation\t"a[1]"\t"$5}}'\
            ../GeneralAnnotation/GeneralInformation\
            $StatisticalFold"miRNA/miRNARegulation" | sort| uniq  >> $StatisticalFold"EditingEvents/GeneEachAnalysis"
awk -F"\t" 'ARGIND==1 && (FNR>1){split($3,a,".");ENSG[$6]=a[1]"\t"$5}\
            ARGIND==2 && (FNR>1) && ($6-$5<-6){print "RNA editings affecting the stability of the RNA structures (from Minimum free energy (MFE))\t"ENSG[$2]}'\
            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
            /data9/swu13/CAeditome_hg38/Structure/Files/RNAstructureInformation | sort | uniq >> $StatisticalFold"EditingEvents/GeneEachAnalysis"
awk -F"\t" '($6>0.5 || $8>0.5 || $10>0.5){split($4,a,".");print "Correlation between ADAR gene expression and RNA editing frequency\t"a[1]"\t"$11}'\
            $StatisticalFold"ADAR/ADARInformation" | sort | uniq >> $StatisticalFold"EditingEvents/GeneEachAnalysis"            
}


#0.3 png files and title
function Title(){
echo -e "PngFileName\tTitleName" > $StatisticalFold"EditingEvents/ENSG.Title"
awk -F"\t" '(NR>1){print $3".png\t"$3"."$5}' /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation >> $StatisticalFold"EditingEvents/ENSG.Title"
echo -e "PngFileName\tTitleName" > $StatisticalFold"EditingEvents/Editing.Title"
awk -F"\t" '(NR>1){print $2".png\t"$1}' /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation >> $StatisticalFold"EditingEvents/Editing.Title"
echo -e "PngFileName\tTitleName" > $StatisticalFold"EditingEvents/EditingENSG.Title"
awk -F"\t" '(NR>1){split($3,a,".");print $2"."a[1]".png\t"$1"."$3"."$5}' /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation >> $StatisticalFold"EditingEvents/EditingENSG.Title"
echo -e "PngFileName\tTitleName" > $StatisticalFold"EditingEvents/ENSG.noversion.Title"
awk -F"\t" '(NR>1){split($3,a,".");print a[1]".png\t"$3"."$5}' /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation >> $StatisticalFold"EditingEvents/ENSG.noversion.Title"
echo -e "PngFileName\tTitleName" > $StatisticalFold"EditingEvents/ADAR.Title"
awk -F"\t" '(NR>1){print "ADAR1/"$3".png\t"$3"."$5".ADAR1";\
                   print "ADAR2/"$3".png\t"$3"."$5".ADAR2";\
                   print "ADAR3/"$3".png\t"$3"."$5".ADAR3";}'\
                   /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation >> $StatisticalFold"EditingEvents/ADAR.Title"
echo -e "PngFileName\tTitleName" > $StatisticalFold"EditingEvents/Splicing.Title"
awk -F"\t" '(NR>1){print $2"."$3".png\t"$1"."$3}' /data9/swu13/CAeditome_hg38/Splicing/Files/SplicingInformation >> $StatisticalFold"EditingEvents/Splicing.Title" 
}

#0.4 sample number figures
function SampleNumber(){
rm -f $StatisticalFold"EditingEvents/SampleNumber.1"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
#sample number
awk -F"\t" -v cancer="$cancer" '(NR==1){Tumor=0;Normal=0;for(i=2;i<=NF;i++){split($i,a,"-");\
                                             if(substr(a[4],1,1)==0){Tumor=Tumor+1}\
                                             else if(substr(a[4],1,1)==1){Normal=Normal+1}}}END{print cancer"\tTumor\t"Tumor"\t"cancer"_Tumor";\
                                                                                                print cancer"\tNormal\t"Normal"\t"cancer"_Normal"}'\
           "/data9/swu13/CAeditome_hg38/RawData/RNAediting/"$cancer"_Editing.txt" >> $StatisticalFold"EditingEvents/SampleNumber.1"
done
echo -e "Cancer\tGroup\tSampleNum\tColorGroup" > $StatisticalFold"EditingEvents/SampleNumber"
cat $StatisticalFold"EditingEvents/SampleNumber.1" | sort | uniq >> $StatisticalFold"EditingEvents/SampleNumber"
Rscript BarPlot.R $StatisticalFold"EditingEvents/SampleNumber" $StatisticalFold"EditingEvents/SampleNumber.png"
}

#1. Specific analysis
function SpecifcEditingAnalysis(){
mkdir $StatisticalFold"Specific/"
rm -f $StatisticalFold"Specific/TumorSpecific"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && (FNR>1){CAeditome[$2]=$1;if(!(Gene[$2]~$3"|"$5)){Gene[$2]=Gene[$2]";"$3"|"$5}}\
                                ARGIND==2{Tumor[$1]=$1}\
                                ARGIND==3 && (($2=="NA" && $5>=5 && $6==0&& $8>=5) || ($2<0.05 && $5+$6>=50)) && (FNR>1){split(Gene[$1],a,";");\
                                            for(i=2;i<=length(a);i++){split(a[i],b,"|");out=CAeditome[$1]"\t"$1"\t"b[1]"\t"b[2]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"cancer;\
                                            if(b[2] in Tumor){out=out"\tY"}else{out=out"\tN"};print out}}'\
                               /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
                               /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
                              "/data9/swu13/CAeditome_hg38/Differential/Files/"$cancer".DifferentialResults" >> $StatisticalFold"Specific/TumorSpecific"
done
#statistical ALL
awk -F"\t" '{print $1}' $StatisticalFold"Specific/TumorSpecific" \
                   | sort | uniq | wc -l > $StatisticalFold"Specific/Statistical"
awk -F"\t" '{print $3}' $StatisticalFold"Specific/TumorSpecific" \
                   | sort | uniq | wc -l >> $StatisticalFold"Specific/Statistical"
awk -F"\t" '($13=="Y"){print $1}' $StatisticalFold"Specific/TumorSpecific" \
                   | sort | uniq | wc -l >> $StatisticalFold"Specific/Statistical"
awk -F"\t" '($13=="Y"){print $3}' $StatisticalFold"Specific/TumorSpecific" \
                   | sort | uniq | wc -l >> $StatisticalFold"Specific/Statistical"
#awk -F"\t" '(($6>$7 && $5<0.05) || ($5=="NA" && $8>=5 && $9==0)){print $1}' ../Statistic/Specific/TumorSpecific | sort | uniq | wc -l
#awk -F"\t" '(($6<$7 && $5<0.05) || ($5=="NA" && $8==0 && $9>=5)){print $1}' ../Statistic/Specific/TumorSpecific | sort | uniq | wc -l
#statistical each cancer
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
#editing
num1=`awk -F"\t" -v cancer="$cancer" '($12==cancer){print $1}' $StatisticalFold"Specific/TumorSpecific" | sort | uniq | wc -l`
#Gene
num2=`awk -F"\t" -v  cancer="$cancer" '($12==cancer){print $3}' $StatisticalFold"Specific/TumorSpecific" | sort | uniq | wc -l`
#editing in tumor-reltaed genes
num3=`awk -F"\t" -v cancer="$cancer" '($12==cancer) && ($13=="Y"){print $1}' $StatisticalFold"Specific/TumorSpecific" | sort | uniq | wc -l`
#tumor-related genes with editing
num4=`awk -F"\t" -v cancer="$cancer" '($12==cancer) && ($13=="Y"){print $3}' $StatisticalFold"Specific/TumorSpecific" | sort | uniq | wc -l`
echo -e $cancer"\t"$num1"\t"$num2"\t"$num3"\t"$num4 >> $StatisticalFold"Specific/Statistical"
done
#Histogram
awk -F"\t" '($5=="NA" && $9==0){print $1"\t"$12"\t"$8"\t"$4"\t"$13}' $StatisticalFold"Specific/TumorSpecific" | sort | uniq > $StatisticalFold"Specific/Statistical.histogram.data"
cut -f1-3 $StatisticalFold"Specific/Statistical.histogram.data" | sort | uniq > $StatisticalFold"Specific/Statistical.histogram.data.figure"
Rscript Histogram.R $StatisticalFold"Specific/Statistical.histogram.data.figure" $StatisticalFold"Specific/Statistical.histogram.png" 6 3
#Bubble plot
awk -F"\t" '($5=="NA" && $9==0 && $13=="Y"){print $1"|"$4"\t"$8"\t"$8"\t"$8/$10*2"\t"$12}' $StatisticalFold"Specific/TumorSpecific" | sort | uniq > $StatisticalFold"Specific/Statistical.bubble.data"
Rscript BubblePlot.R $StatisticalFold"Specific/Statistical.bubble.data" 0 20 $StatisticalFold"Specific/Statistical.bubble.png" 6 3 
cancers=`cut -f2 $StatisticalFold"Specific/Statistical.histogram.data" | sort | uniq`
num=`echo $cancers | awk -F" " '{print NF}'`
for((i=1;i<=$num;i++))
do
cancer=`echo $cancers | cut -d" " -f$i`
awk -F"\t" -v cancer="$cancer" '($5=="Y") && ($2==cancer){print $0}' $StatisticalFold"Specific/Statistical.histogram.data" | sort -nrk 3 | head -n 1
done
#Volcano plot
awk -F"\t" '($5<0.05){print $12"\t"$1"\t"$5"\t"$6-$7}' $StatisticalFold"Specific/TumorSpecific" | sort | uniq > $StatisticalFold"Specific/Statistical.volcano.data"
Rscript VolcanoPlot.R $StatisticalFold"Specific/Statistical.volcano.data" $StatisticalFold"Specific/Statistical.volcano.png" 0.01 0.1 5 4 10
#Bubble plot
awk -F"\t" '($5<0.05) && ($13=="Y"){if($6>$7){size=$6-$7}else{size=$7-$6};print $1"|"$4"\t"log($5)"\t"$6-$7"\t"size"\t"$12}' $StatisticalFold"Specific/TumorSpecific" | sort | uniq > $StatisticalFold"Specific/Statistical.volcano.bubble.data"
Rscript BubblePlot.R $StatisticalFold"Specific/Statistical.volcano.bubble.data" 0 20 $StatisticalFold"Specific/Statistical.volcano.bubble.png" 5 4
#select more significant cases for Database
awk -F"\t" '($5<0.05){print "/data9/swu13/CAeditome_hg38/Differential/Boxplot/"$12"/"$2".png"}' $StatisticalFold"Specific/TumorSpecific" | sort | uniq > $StatisticalFold"Specific/SelectedBoxplot"
}

#2. Stage analysis
function StageEditingAnalysis(){
mkdir $StatisticalFold"Stage/"
echo -e "CAeditomeID\tRNAediting\tENSG\tGeneName\tPvalue\tCorrelationR\tCancer\tStageType\tTumorRelated" > $StatisticalFold"Stage/StageAssociation"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
if [[ -f "/data9/swu13/CAeditome_hg38/StageAssociation/Files/"$cancer".PathStage.CorResults" ]]
then
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && (FNR>1){CAeditome[$2]=$1;if(!(Gene[$2]~$3"|"$5)){Gene[$2]=Gene[$2]";"$3"|"$5}}\
                                ARGIND==2{Tumor[$1]=$1}\
                                ARGIND==3 && ($5>=50) && (FNR>1){Specific[$1]=$1}\
                                ARGIND==4 && ($2<0.05) && ($1 in Specific) && (FNR>1){split(Gene[$1],a,";");\
                                                                for(i=2;i<=length(a);i++){split(a[i],b,"|");out=CAeditome[$1]"\t"$1"\t"b[1]"\t"b[2]"\t"$2"\t"$3"\t"cancer"\tpathologic_stage";\
                                                                if(b[2] in Tumor){out=out"\tY"}else{out=out"\tN"};print out}}'\
                                /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
                                /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
                                "/data9/swu13/CAeditome_hg38/Differential/Files/"$cancer".DifferentialResults"\
                                "/data9/swu13/CAeditome_hg38/StageAssociation/Files/"$cancer".PathStage.CorResults" >> $StatisticalFold"Stage/StageAssociation"
fi
if [[ -f "/data9/swu13/CAeditome_hg38/StageAssociation/Files/"$cancer".CliStage.CorResults" ]]
then
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && (FNR>1){CAeditome[$2]=$1;if(!(Gene[$2]~$3"|"$5)){Gene[$2]=Gene[$2]";"$3"|"$5}}\
                                ARGIND==2{Tumor[$1]=$1}\
                                ARGIND==3 && ($5>=50) && (FNR>1){Specific[$1]=$1}\
                                ARGIND==4 && ($2<0.05) && ($1 in Specific) && (FNR>1){split(Gene[$1],a,";");\
                                                                for(i=2;i<=length(a);i++){split(a[i],b,"|");out=CAeditome[$1]"\t"$1"\t"b[1]"\t"b[2]"\t"$2"\t"$3"\t"cancer"\tclinical_stage";\
                                                                if(b[2] in Tumor){out=out"\tY"}else{out=out"\tN"};print out}}'\
                                /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
                                /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
                                "/data9/swu13/CAeditome_hg38/Differential/Files/"$cancer".DifferentialResults"\
                                "/data9/swu13/CAeditome_hg38/StageAssociation/Files/"$cancer".CliStage.CorResults">> $StatisticalFold"Stage/StageAssociation"
fi
done
echo -e "Cancer\tStageType\tCAeditomeID\tRNAediting\tENSG\tGeneName\tPvalue\tCorrelationR" > $StatisticalFold"Stage/StageEditingAssociation"
awk -F"\t" '(NR>1){print $7"\t"$8"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $StatisticalFold"Stage/StageAssociation" | sort | uniq >> $StatisticalFold"Stage/StageEditingAssociation"
#All statistical
awk -F"\t" '($5<0.05) && (NR>1){print $1}' $StatisticalFold"Stage/StageAssociation" \
                   | sort | uniq | wc -l > $StatisticalFold"Stage/Statistical"
awk -F"\t" '($5<0.05) && (NR>1){print $3}' $StatisticalFold"Stage/StageAssociation" \
                   | sort | uniq | wc -l >> $StatisticalFold"Stage/Statistical"
awk -F"\t" '($5<0.05) && (NR>1) && ($9=="Y"){print $1}' $StatisticalFold"Stage/StageAssociation" \
                   | sort | uniq | wc -l >> $StatisticalFold"Stage/Statistical"
awk -F"\t" '($5<0.05) && (NR>1) && ($9=="Y"){print $3}' $StatisticalFold"Stage/StageAssociation" \
                   | sort | uniq | wc -l >> $StatisticalFold"Stage/Statistical"
awk -F"\t" '($6>0 && NR>1){print $1}' ../Statistic/Stage/StageAssociation | sort |uniq | wc -l
awk -F"\t" '($6<0 && NR>1){print $1}' ../Statistic/Stage/StageAssociation | sort |uniq | wc -l
#statistical each cancer
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
#editing
num2=`awk -F"\t" -v cancer="$cancer" '($7==cancer) && (NR>1) && ($5<0.05){print $1}' $StatisticalFold"Stage/StageAssociation" | sort | uniq | wc -l`
#Gene
num3=`awk -F"\t" -v  cancer="$cancer" '($7==cancer) && (NR>1)&& ($5<0.05){print $3}' $StatisticalFold"Stage/StageAssociation" | sort | uniq | wc -l`
#editing in tumor-reltaed genes
num4=`awk -F"\t" -v cancer="$cancer" '($7==cancer) && (NR>1) && ($9=="Y") && ($5<0.05){print $1}' $StatisticalFold"Stage/StageAssociation" | sort | uniq | wc -l`
#tumor-related genes with editing
num5=`awk -F"\t" -v cancer="$cancer" '($7==cancer) && (NR>1) && ($9=="Y") && ($5<0.05){print $3}' $StatisticalFold"Stage/StageAssociation" | sort | uniq | wc -l`
echo -e $cancer"\t"$num2"\t"$num3"\t"$num4"\t"$num5 >> $StatisticalFold"Stage/Statistical"
done
#Volcano plot
awk -F"\t" '($5<0.05){print $7"\t"$1"\t"$5"\t"$6}' $StatisticalFold"Stage/StageAssociation" | sort | uniq > $StatisticalFold"Stage/Statistical.volcano.data"
Rscript VolcanoPlot.R $StatisticalFold"Stage/Statistical.volcano.data" $StatisticalFold"Stage/Statistical.volcano.png" 0.01 0.1 6 4 5
#Bubble plot
awk -F"\t" '($5<0.05) && ($8=="pathologic_stage") && ($9=="Y"){if($6>0){size=$6}else{size=0-$6};print $1"|"$4"\t"log($5)"\t"$6"\t"size"\t"$7}' $StatisticalFold"Stage/StageAssociation"\
    | sort | uniq > $StatisticalFold"Stage/Statistical.volcano.bubble.data"
unset cancers
cancers=`cut -f5 $StatisticalFold"Stage/Statistical.volcano.bubble.data" | sort | uniq`
awk -F"\t" -v cancers="${cancers[*]}" '($5<0.05) && ($8=="clinical_stage") && ($9=="Y"){split(cancers,mm," ");for(i in mm){nn[mm[i]]=nn[mm[i]]};\
             if(!($7 in nn)){if($6>0){size=$6}else{size=0-$6};print $1"|"$4"\t"log($5)"\t"$6"\t"size"\t"$7}}'\
             $StatisticalFold"Stage/StageAssociation" | sort | uniq >> $StatisticalFold"Stage/Statistical.volcano.bubble.data"
Rscript BubblePlot.R $StatisticalFold"Stage/Statistical.volcano.bubble.data" 0 20 $StatisticalFold"Stage/Statistical.volcano.bubble.png" 6 4
}

#2.1 stage association with gene expression
function StageExpressionAnalysis(){
echo -e "Cancer\tStageType\tENSG\tGeneName\tPvalue\tCorrelationR" > $StatisticalFold"Stage/StageExpressionAssociation"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
if [[ -f "/data9/swu13/CAeditome_hg38/StageExpressionAssociation/Files/"$cancer".PathStage.CorResults" ]]
then
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && ($3=="gene"){gsub(/"| /,"",$9);split($9,a,";");ENSG=substr(a[1],length("gene_id")+1);split(ENSG,b,".");\
                                                          GeneName=substr(a[4],length("gene_name")+1);Gene[b[1]]=GeneName}\
                                ARGIND==2 && ($2<0.05) && ($3>0.3 || $3<-0.3) && (FNR>1){print cancer"\tpathologic_stage\t"$1"\t"Gene[$1]"\t"$2"\t"$3}'\
                                /data9/swu13/CAeditome_hg38/RawData/Reference/gencode.v22.annotation.gtf\
                                "/data9/swu13/CAeditome_hg38/StageExpressionAssociation/Files/"$cancer".PathStage.CorResults" >> $StatisticalFold"Stage/StageExpressionAssociation"
fi
if [[ -f "/data9/swu13/CAeditome_hg38/StageExpressionAssociation/Files/"$cancer".CliStage.CorResults" ]]
then
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && ($3=="gene"){gsub(/"| /,"",$9);split($9,a,";");ENSG=substr(a[1],length("gene_id")+1);split(ENSG,b,".");\
                                                          GeneName=substr(a[4],length("gene_name")+1);Gene[b[1]]=GeneName}\
                                ARGIND==2 && ($2<0.05) && ($3>0.3 || $3<-0.3) && (FNR>1){print cancer"\tclinical_stage\t"$1"\t"Gene[$1]"\t"$2"\t"$3}'\
                                /data9/swu13/CAeditome_hg38/RawData/Reference/gencode.v22.annotation.gtf\
                                "/data9/swu13/CAeditome_hg38/StageExpressionAssociation/Files/"$cancer".CliStage.CorResults">> $StatisticalFold"Stage/StageExpressionAssociation"
fi
done
}


#3. Survival analysis
function SurvivalEditingAnalysis(){
mkdir $StatisticalFold"Survival/"
echo -e "CAeditomeID\tRNAediting\tENSG\tGeneName\tKMpvalue\tCoxpvalue\tCoxHR\tCancer\tTumorRelated" >  $StatisticalFold"Survival/SurvivalAssociation"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && (FNR>1){CAeditome[$2]=$1;if(!(Gene[$2]~$3"|"$5)){Gene[$2]=Gene[$2]";"$3"|"$5}}\
                                ARGIND==2{Tumor[$1]=$1}\
                                ARGIND==3 && ($5>=50) && (FNR>1){Specific[$1]=$1}\
                                ARGIND==4 && ($2<0.05) && ($7<0.05) && ($1 in Specific) && (FNR>1){split(Gene[$1],a,";");\
                                                                for(i=2;i<=length(a);i++){split(a[i],b,"|");out=CAeditome[$1]"\t"$1"\t"b[1]"\t"b[2]"\t"$2"\t"$7"\t"$8"\t"cancer;\
                                                                if(b[2] in Tumor){out=out"\tY"}else{out=out"\tN"};print out}}'\
                                /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
                                /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
                                "/data9/swu13/CAeditome_hg38/Differential/Files/"$cancer".DifferentialResults"\
                                "/data9/swu13/CAeditome_hg38/SurvivalAssociation/Files/"$cancer".survival.results" >> $StatisticalFold"Survival/SurvivalAssociation"
done
echo -e "Cancer\tCAeditomeID\tRNAediting\tENSG\tGeneName\tKMpvalue\tCoxpvalue\tCoxHR" > $StatisticalFold"Survival/SurvivalEditingAssociation"
awk -F"\t" '(NR>1){out=$8;for(i=1;i<=7;i++){out=out"\t"$i};print out}' $StatisticalFold"Survival/SurvivalAssociation" >> $StatisticalFold"Survival/SurvivalEditingAssociation"
#All statistical
awk -F"\t" '($5<0.05) && (NR>1) && ($6<0.05){print $1}' $StatisticalFold"Survival/SurvivalAssociation" \
                   | sort | uniq | wc -l > $StatisticalFold"Survival/Statistical"
awk -F"\t" '($5<0.05) && (NR>1) && ($6<0.05){print $3}' $StatisticalFold"Survival/SurvivalAssociation" \
                   | sort | uniq | wc -l >> $StatisticalFold"Survival/Statistical"
awk -F"\t" '($5<0.05) && (NR>1) && ($6<0.05) && ($9=="Y"){print $1}' $StatisticalFold"Survival/SurvivalAssociation" \
                   | sort | uniq | wc -l >> $StatisticalFold"Survival/Statistical"
awk -F"\t" '($5<0.05) && (NR>1) && ($6<0.05) && ($9=="Y"){print $3}' $StatisticalFold"Survival/SurvivalAssociation" \
                   | sort | uniq | wc -l >> $StatisticalFold"Survival/Statistical"
#awk -F"\t" '($7>1 && NR>1){print $1}' ../Statistic/Survival/SurvivalAssociation | sort | uniq | wc -l
#awk -F"\t" '($7<1 && NR>1){print $1}' ../Statistic/Survival/SurvivalAssociation | sort | uniq | wc -l
#statistical each cancer
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
#editing
num2=`awk -F"\t" -v cancer="$cancer" '($8==cancer) && (NR>1) && ($5<0.05) && ($6<0.05){print $1}' $StatisticalFold"Survival/SurvivalAssociation" | sort | uniq | wc -l`
#Gene
num3=`awk -F"\t" -v  cancer="$cancer" '($8==cancer) && (NR>1) && ($5<0.05) && ($6<0.05){print $3}' $StatisticalFold"Survival/SurvivalAssociation" | sort | uniq | wc -l`
#editing in tumor-reltaed genes
num4=`awk -F"\t" -v cancer="$cancer" '($8==cancer) && (NR>1) && ($9=="Y") && ($5<0.05) && ($6<0.05){print $1}' $StatisticalFold"Survival/SurvivalAssociation" | sort | uniq | wc -l`
#tumor-related genes with editing
num5=`awk -F"\t" -v cancer="$cancer" '($8==cancer) && (NR>1) && ($9=="Y") && ($5<0.05) && ($6<0.05){print $3}' $StatisticalFold"Survival/SurvivalAssociation" | sort | uniq | wc -l`
echo -e $cancer"\t"$num2"\t"$num3"\t"$num4"\t"$num5 >> $StatisticalFold"Survival/Statistical"
done
#Volcano plot
awk -F"\t" '($5<0.05)&& ($6<0.05){if($7>20){HR=20}else{HR=$7};print $8"\t"$1"\t"$6"\t"HR-1}' $StatisticalFold"Survival/SurvivalAssociation" | sort | uniq > $StatisticalFold"Survival/Statistical.volcano.data"
Rscript VolcanoPlot.R $StatisticalFold"Survival/Statistical.volcano.data" $StatisticalFold"Survival/Statistical.volcano.png" 0.01 0.5 6 5
#Bubble plot
awk -F"\t" '($5<0.05)&& ($6<0.05)&& ($9=="Y"){if($7>20){HR=20}else{HR=$7};print $1"|"$4"\t"log($6)"\t"HR"\t"(-log($6))"\t"$8}' $StatisticalFold"Survival/SurvivalAssociation" | sort | uniq > $StatisticalFold"Survival/Statistical.volcano.bubble.data"
Rscript BubblePlot.R $StatisticalFold"Survival/Statistical.volcano.bubble.data" 1 20 $StatisticalFold"Survival/Statistical.volcano.bubble.png" 6 4
#select more significant cases for Database
awk -F"\t" '{print "/data9/swu13/CAeditome_hg38/SurvivalAssociation/KMplot/"$8"/"$2".png"}' $StatisticalFold"Survival/SurvivalAssociation" | sort | uniq  > $StatisticalFold"Survival/SelectedKMplot"
#insert GRIA2 sample
echo -e "UCEC\tCAediting_390714\tchr4_157360142_+\tENSG00000120251.17\tGRIA2\t0.0482098696638136\t0.00962948\t2329240" > /data9/swu13/CAeditome_hg38/Statistic/Survival/SurvivalGRIA2Association
}

#3.1 survival for gene expression
function SurvivalExpressionAnalysis(){
echo -e "Cancer\tENSG\tGeneName\tKMpvalue\tCoxpvalue\tCoxHR" >  $StatisticalFold"Survival/SurvivalExpressionAssociation"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && ($3=="gene"){gsub(/"| /,"",$9);split($9,a,";");ENSG=substr(a[1],length("gene_id")+1);split(ENSG,b,".");\
                                                          GeneName=substr(a[4],length("gene_name")+1);Gene[b[1]]=GeneName}\
                                ARGIND==2 && ($2<0.05) && ($7<0.05) && (FNR>1){print cancer"\t"$1"\t"Gene[$1]"\t"$2"\t"$7"\t"$8}'\
                                /data9/swu13/CAeditome_hg38/RawData/Reference/gencode.v22.annotation.gtf\
                                "/data9/swu13/CAeditome_hg38/SurvivalExpressionAssociation/Files/"$cancer".survival.results" >> $StatisticalFold"Survival/SurvivalExpressionAssociation"
done
}

#4. DEG analysis
function DEGAnalysis(){
mkdir $StatisticalFold"DEG/"
rm -f $StatisticalFold"DEG/EditedNonEditedDEGs"
rm -f $StatisticalFold"DEG/TumorNormalDEGs" 
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && (FNR>1){CAeditome[$2]=$1;split($3,a,".");Gene[a[1]]=$5}\
                                ARGIND==2{Tumor[$1]=$1}\
                                ARGIND==3 && ($7>=50) && ($5>=25) && ($7-$5>=25) && (FNR>1){Specific[$1]=$1}\
                                ARGIND==4 && ($3<0.05) && ($7<0.05) &&($4>1 || $5>1) && ($1 in Specific) && (FNR>1)  && ($6*$8>0)\
                                                         {out=CAeditome[$1]"\t"$1"\t"$2"\t"Gene[$2]"\t"$3"\t"$6"\t"$7"\t"$8"\t"cancer;\
                                                                if(Gene[$2] in Tumor){out=out"\tY"}else{out=out"\tN"};print out}'\
                                /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
                                /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
                                "/data9/swu13/CAeditome_hg38/Differential/Files/"$cancer".DifferentialResults"\
                                "/data9/swu13/CAeditome_hg38/ExpressionCorrelation/Files/"$cancer".Results" >> $StatisticalFold"DEG/EditedNonEditedDEGs"
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && ($3=="gene"){gsub(/"| /,"",$9);split($9,a,";");ENSG=substr(a[1],length("gene_id")+1);split(ENSG,b,".");\
                                                          GeneName=substr(a[4],length("gene_name")+1);Gene[b[1]]=GeneName}\
                                ARGIND==2{Tumor[$1]=$1}\
                                ARGIND==3 && ($2<0.05) && (FNR>1) && ($6+$7>=50) && ($3>1 || $4>1){out=$1"\t"Gene[$1]"\t"$2"\t"$5"\t"cancer;\
                                                                if(Gene[$1] in Tumor){out=out"\tY"}else{out=out"\tN"};print out}'\
                                /data9/swu13/CAeditome_hg38/RawData/Reference/gencode.v22.annotation.gtf\
                                /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
                                "/data9/swu13/CAeditome_hg38/DifferentialExpression/Files/"$cancer".DifferentialResults" >> $StatisticalFold"DEG/TumorNormalDEGs" 
done
#Statistical
awk -F"\t" '($5<0.05) && ($6>0.3) {print $1}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | wc -l > $StatisticalFold"DEG/Statistic"
awk -F"\t" '($5<0.05) && ($6>0.3) {print $3}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | wc -l >> $StatisticalFold"DEG/Statistic"
awk -F"\t" '($5<0.05) && ($6<-0.3){print $1}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | wc -l >> $StatisticalFold"DEG/Statistic"
awk -F"\t" '($5<0.05) && ($6<-0.3){print $3}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | wc -l >> $StatisticalFold"DEG/Statistic"
#editing
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4>0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6>0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) {print $1}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4>0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6<-0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) {print $1}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4<-0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6>0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) {print $1}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4<-0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6<-0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) {print $1}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4>0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6>0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) && ($10=="Y"){print $1}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4>0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6<-0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) && ($10=="Y") {print $1}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4<-0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6>0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) && ($10=="Y") {print $1}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4<-0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6<-0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) && ($10=="Y") {print $1}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
#gene
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4>0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6>0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) {print $3}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4>0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6<-0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) {print $3}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4<-0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6>0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) {print $3}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4<-0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6<-0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) {print $3}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4>0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6>0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) && ($10=="Y"){print $3}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4>0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6<-0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) && ($10=="Y") {print $3}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4<-0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6>0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) && ($10=="Y") {print $3}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4<-0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6<-0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) && ($10=="Y") {print $3}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l  >> $StatisticalFold"DEG/Statistic"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4>0.3 || $4<-0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6>0.3 || $6<-0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) {print $1}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4>0.3 || $4<-0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6>0.3 || $6<-0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) {print $3}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4>0.3 || $4<-0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6>0.3 || $6<-0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) && ($10=="Y") {print $1}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4>0.3 || $4<-0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6>0.3 || $6<-0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) && ($10=="Y") {print $3}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | sort | uniq | wc -l
#ALL statistical for table 1
awk -F"\t" '($5<0.05) && ($6>0.3 || $6<-0.3) && ($6*$8>0){print $1}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | wc -l
awk -F"\t" '($5<0.05) && ($6>0.3 || $6<-0.3) && ($6*$8>0){print $3}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | wc -l
awk -F"\t" '($5<0.05) && ($6>0.3 || $6<-0.3) && ($6*$8>0) && ($10=="Y"){print $1}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | wc -l
awk -F"\t" '($5<0.05) && ($6>0.3 || $6<-0.3) && ($6*$8>0) && ($10=="Y"){print $3}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | wc -l
#Statistical for each cancer
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
num1=`awk -F"\t" -v cancer="$cancer" '($5<0.05) && ($6>0.3 || $6<-0.3) && ($6*$8>0) && ($9==cancer){print $1}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | wc -l`
num2=`awk -F"\t" -v cancer="$cancer" '($5<0.05) && ($6>0.3 || $6<-0.3) && ($6*$8>0) && ($9==cancer){print $3}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | wc -l`
num3=`awk -F"\t" -v cancer="$cancer" '($5<0.05) && ($6>0.3 || $6<-0.3) && ($6*$8>0)  && ($9==cancer) && ($10=="Y"){print $1}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | wc -l`
num4=`awk -F"\t" -v cancer="$cancer" '($5<0.05) && ($6>0.3 || $6<-0.3) && ($6*$8>0)  && ($9==cancer) && ($10=="Y"){print $3}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  | wc -l`
echo -e $cancer"\t"$num1"\t"$num2"\t"$num3"\t"$num4 >> $StatisticalFold"DEG/Statistic"
done
#volcano plot
awk -F"\t" '($5<0.05) && ($6*$8>0){print $9"\t"$1"|"$4"\t"$5"\t"$6}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq > $StatisticalFold"DEG/Statistical.volcano.data"
Rscript VolcanoPlot.R $StatisticalFold"DEG/Statistical.volcano.data" $StatisticalFold"DEG/Statistical.volcano.png" 0.01 0.3 8 4 10
#Bubble plot
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4<-0.3 || $4>0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6<-0.3 || $6>0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) && ($10=="Y") {if($6>0){size=$6}else{size=0-$6};\
                                                                                                                   print $1"|"$4"\t"log($5)"\t"$6"\t"size"\t"$9}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq > $StatisticalFold"DEG/Statistical.volcano.bubble.data"
Rscript BubblePlot.R $StatisticalFold"DEG/Statistical.volcano.bubble.data" 0 20 $StatisticalFold"DEG/Statistical.volcano.bubble.png" 5 4
#enrichment
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($4>0.3 || $4<-0.3){TumorNormal[$1"\t"$5]=$1"\t"$5}\
            ARGIND==2 && ($5<0.05) && ($6<-0.3 || $6>0.3) && ($6*$8>0) && ($3"\t"$9 in TumorNormal) {print $4"\t"$9}'\
            $StatisticalFold"DEG/TumorNormalDEGs"\
            $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq > $StatisticalFold"DEG/Statistical.enrichment.data"
Rscript Enrichment.R $StatisticalFold"DEG/Statistical.enrichment.data" $StatisticalFold"DEG/Statistical.enrichment.KEGG" $StatisticalFold"DEG/Statistical.enrichment.KEGG" 4 0.05 0.2 4 4 
#select more significant cases for Database
awk -F"\t" '($6<-0.3 || $6>0.3){print "/data9/swu13/CAeditome_hg38/ExpressionCorrelation/Boxplot/"$9"/"$2"."$3".png"}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  > $StatisticalFold"DEG/SelectedBoxplot"
#awk -F"\t" '($8<-0.3 || $8>0.3){print "/data9/swu13/CAeditome_hg38/ExpressionCorrelation/Correlation/"$9"/"$2"."$3".png"}' $StatisticalFold"DEG/EditedNonEditedDEGs" | sort | uniq  > $StatisticalFold"DEG/SelectedCorrelation"
#Title for each figure
awk -F"/" 'ARGIND==1 && (FNR>1){split($0,m,"\t");split(m[3],a,".");GeneName[a[1]]=m[5];CAeditomeID[m[2]]=m[1]}\
           ARGIND==2{split($NF,a,".");print $0"\t"$(NF-1)"."CAeditomeID[a[1]]"."a[2]"."GeneName[a[2]]}'\
           /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation $StatisticalFold"DEG/SelectedBoxplot" > $StatisticalFold"DEG/Title"
num=`wc -l $StatisticalFold"DEG/Title" | cut -d" " -f1`
for((i=1;i<=$num;i++))
do
originalPath=`head -n $i $StatisticalFold"DEG/Title" | tail -n 1 | cut -f1`
filename=`head -n $i $StatisticalFold"DEG/Title" | tail -n 1 | cut -f2`
cp $originalPath "/data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedBoxplot/"$filename".png"
done
#extract GRIA2 figures for Database
#grep ENSG00000120251 ../Statistic/DEG/EditedNonEditedDEGs
#grep chr4_157360142_+ /data9/swu13/CAeditome_hg38/ExpressionCorrelation/Tmp/PCPG.RNAediting.tumor > /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/chr4_157360142_+.editing
#grep ENSG00000120251 /data9/swu13/CAeditome_hg38/ExpressionCorrelation/Tmp/PCPG.Expression.tumor | head -n 1 >/data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/GRIA2.expression
#Rscript Correlation.R /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/chr4_157360142_+.editing /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/GRIA2.expression /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/boxplot /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/correlation "FALSE"
#cp /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/boxplot/* /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/"PCPG.CAediting_390714.ENSG00000120251.GRIA2.png"
#rm -f /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/boxplot/*
#grep chr4_157336984_+ /data9/swu13/CAeditome_hg38/ExpressionCorrelation/Tmp/LGG.RNAediting.tumor > /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/chr4_157336984_+.editing
#grep ENSG00000120251 /data9/swu13/CAeditome_hg38/ExpressionCorrelation/Tmp/LGG.Expression.tumor | head -n 1 > /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/GRIA2.expression
#Rscript Correlation.R /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/chr4_157336984_+.editing /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/GRIA2.expression /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/boxplot /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/correlation "FALSE"
#cp /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/boxplot/chr4_157336984_* /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/"LGG.CAediting_390710.ENSG00000120251.GRIA2.png"
#rm -f /data9/swu13/CAeditome_hg38/ExpressionCorrelation/SelectedGRIA2/boxplot/*

}

#5. Protein analysis
function ProteinAnalysis(){
mkdir $StatisticalFold"Protein/"
#All statistical
awk -F"\t" '($3=="nonsynonymous SNV"){print $1}' /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinLolliplot  |sort | uniq | wc -l > $StatisticalFold"Protein/Sequence"
awk -F"\t" '($3=="synonymous SNV"){print $1}' /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinLolliplot  |sort | uniq | wc -l >> $StatisticalFold"Protein/Sequence"
awk -F"\t" '($3=="stoploss"){print $1}' /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinLolliplot  |sort | uniq | wc -l >> $StatisticalFold"Protein/Sequence"
awk -F"\t" '($3=="nonsynonymous SNV") || ($3=="stoploss"){print $6}' /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinLolliplot  |sort | uniq | wc -l >> $StatisticalFold"Protein/Sequence"
#awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
#            ARGIND==2 && (($3=="nonsynonymous SNV")) && ($6 in Tumor){print $1}'\
#            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
#            /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinLolliplot  |sort | uniq | wc -l
#awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
#            ARGIND==2 && (($3=="stoploss")) && ($6 in Tumor){print $1}'\
#            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
#            /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinLolliplot  |sort | uniq | wc -l
#awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
#            ARGIND==2 && (($3=="nonsynonymous SNV") || ($3=="stoploss")) && ($6 in Tumor){print $4}'\
#            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
#            /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinLolliplot  |sort | uniq | wc -l
##function
awk -F"\t" '($9=="D"){print $1"\t"$4"\tSIFT"}' /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinFunction | sort | uniq > $StatisticalFold"Protein/Function"
awk -F"\t" '($10=="D" || $11=="D"){print $1"\t"$4"\tPolyphen2"}' /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinFunction | sort | uniq >> $StatisticalFold"Protein/Function"
awk -F"\t" '($12=="D"){print $1"\t"$4"\tPROVEAN"}' /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinFunction | sort | uniq >> $StatisticalFold"Protein/Function"
cut -f1 $StatisticalFold"Protein/Function" | sort | uniq | wc -l
cut -f2 $StatisticalFold"Protein/Function" | sort | uniq | wc -l
#awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
#            ARGIND==2 && ($9=="D" || $10=="D" || $11=="D" || $12=="D") && ($5 in Tumor){print $1}'\
#            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
#            /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinFunction | sort | uniq | wc -l
#awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
#            ARGIND==2 && ($9=="D" || $10=="D" || $11=="D" || $12=="D") && ($5 in Tumor){print $4}'\
#            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
#            /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinFunction | sort | uniq | wc -l
Rscript VennPlot.R $StatisticalFold"Protein/Function" $StatisticalFold"Protein/Function.editing.venn.png" "TRUE"
awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
            ARGIND==2 && ($9=="D") && ($5 in Tumor){print $1"\t"$4"\tSIFT"}'\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinFunction | sort | uniq > $StatisticalFold"Protein/Function"
awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
            ARGIND==2 && ($10=="D" || $11=="D") && ($5 in Tumor){print $1"\t"$4"\tPolyphen2"}'\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinFunction | sort | uniq >> $StatisticalFold"Protein/Function"
awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
            ARGIND==2 && ($12=="D") && ($5 in Tumor){print $1"\t"$4"\tPROVEAN"}'\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinFunction | sort | uniq >> $StatisticalFold"Protein/Function"
Rscript VennPlot.R $StatisticalFold"Protein/Function" $StatisticalFold"Protein/Function.editing.venn.tumor.png" "TRUE"
cut -f1 $StatisticalFold"Protein/Function" | sort | uniq | wc -l
cut -f2 $StatisticalFold"Protein/Function" | sort | uniq | wc -l
awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
            ARGIND==2 && ($9=="D") && ($10=="D" || $11=="D") && ($12=="D") && ($5 in Tumor){print $5}'\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinFunction | sort | uniq 
#Each cancer for sequence
echo -e "Cancer\tGeneName\tUniprotName\tCAeditomeID\tRNAediting\tEditingType\tEditingProteinPosition" > $StatisticalFold"Protein/SequenceCancer"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && (FNR>1){Editing[$1]=$1}\
                                ARGIND==2 && (FNR>1) && ($2 in Editing) && ($10 == $16){print cancer"\t"$6"\t"$15"\t"$1"\t"$2"\t"$3"\t"substr($14,3,length($14)-2)}'\
                                "/data9/swu13/CAeditome_hg38/GeneralAnnotation/"$cancer".Table"\
                                /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinLolliplot >> $StatisticalFold"Protein/SequenceCancer"
done
head -n 1 $StatisticalFold"Protein/SequenceCancer" | cut -f1 --complement  > $StatisticalFold"Protein/SequenceAll"
tail -n +2 $StatisticalFold"Protein/SequenceCancer" | cut -f1 --complement | sort | uniq >> $StatisticalFold"Protein/SequenceAll"
#File for lolliplot
rm -f $StatisticalFold"Protein/EditedNum"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
awk -F"\t" '(NR>1){print $1"\t"$5+$6}' "/data9/swu13/CAeditome_hg38/Differential/Files/"$cancer".DifferentialResults" >> $StatisticalFold"Protein/EditedNum"
done
echo -e "GeneName\tAAchange\tEditedNumber" > $StatisticalFold"Protein/SequenceAll.lolliplot"
awk -F"\t" 'ARGIND==1{Ed[$1]=Ed[$1]+$2}\
            ARGIND==2 && (FNR>1){print $1"\t"$6"\t"Ed[$4]}'\
            $StatisticalFold"Protein/EditedNum"\
            $StatisticalFold"Protein/SequenceAll" >> $StatisticalFold"Protein/SequenceAll.lolliplot"
#Each cancer for function
head -n 1 /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinFunction | awk '{print "Cancer\t"$0"\tTumorRelatedGene"}' > $StatisticalFold"Protein/FunctionCancer"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && (FNR>1){Editing[$1]=$1}\
                                ARGIND==2{Tumor[$1]=$1}\
                                ARGIND==3 && (FNR>1) && ($2 in Editing) && (($9=="D") || ($10=="D") || ($11=="D") || ($12=="D"))\
                                             {out=cancer"\t"$0;if($5 in Tumor){out=out"\tY"}else{out=out"\tN"};print out}'\
                                "/data9/swu13/CAeditome_hg38/GeneralAnnotation/"$cancer".Table"\
                                /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
                                /data9/swu13/CAeditome_hg38/ProteinAnnotation/ProteinFunction >> $StatisticalFold"Protein/FunctionCancer"
done 
head -n 1 $StatisticalFold"Protein/FunctionCancer" | cut -f1,14 --complement > $StatisticalFold"Protein/FunctionAll"
tail -n +2 $StatisticalFold"Protein/FunctionCancer" | cut -f1,14 --complement | sort | uniq >> $StatisticalFold"Protein/FunctionAll"
#Statistical
rm -f $StatisticalFold"Protein/Statistical"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
num1=`awk -F"\t" -v cancer="$cancer" '($1==cancer){print $2}' $StatisticalFold"Protein/FunctionCancer" | sort | uniq | wc -l`
num2=`awk -F"\t" -v cancer="$cancer" '($1==cancer){print $5}' $StatisticalFold"Protein/FunctionCancer" | sort | uniq | wc -l`
num3=`awk -F"\t" -v cancer="$cancer" '($1==cancer) && ($14=="Y"){print $2}' $StatisticalFold"Protein/FunctionCancer" | sort | uniq | wc -l`
num4=`awk -F"\t" -v cancer="$cancer" '($1==cancer) && ($14=="Y"){print $5}' $StatisticalFold"Protein/FunctionCancer" | sort | uniq | wc -l`
echo -e $cancer"\t"$num1"\t"$num2"\t"$num3"\t"$num4 >> $StatisticalFold"Protein/Statistical"
done
awk -F"\t" '(NR>1){print $2}' $StatisticalFold"Protein/FunctionCancer" | sort | uniq | wc -l
awk -F"\t" '(NR>1){print $5}' $StatisticalFold"Protein/FunctionCancer" | sort | uniq | wc -l
awk -F"\t" '(NR>1) && ($14=="Y"){print $2}' $StatisticalFold"Protein/FunctionCancer" | sort | uniq | wc -l
awk -F"\t" '(NR>1) && ($14=="Y"){print $5}' $StatisticalFold"Protein/FunctionCancer" | sort | uniq | wc -l
}

#6. Region/Repeat analysis
function DistributionAnalysis(){
mkdir $StatisticalFold"Distribution/"
#Type distribution
Types=`cut -f4 /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation | tail -n +2 | sort | uniq`
Types=("protein_coding" "lincRNA" "antisense" "processed_transcript" "transcribed_unprocessed_pseudogene")
Allnum=${#Types[@]}
rm -f $StatisticalFold"Distribution/previoustype"
touch $StatisticalFold"Distribution/previoustype"
for((i=0;i<$Allnum;i++))
do
type=${Types[$i]}
awk -F"\t" -v type="$type" 'ARGIND==1{PreviousType[$1]=$1}\
                            ARGIND==2 && ($4==type) && (!($1 in PreviousType)){print $1}'\
                            $StatisticalFold"Distribution/previoustype"\
                            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation | sort | uniq > $StatisticalFold"Distribution/currenttype"
num=`wc -l $StatisticalFold"Distribution/currenttype" |cut -d" " -f1`
echo -e $type"\t"$num >> $StatisticalFold"Distribution/Types"
cat $StatisticalFold"Distribution/currenttype" >> $StatisticalFold"Distribution/previoustype"
done
awk -F"\t" -v type="$type" 'ARGIND==1{PreviousType[$1]=$1}\
                            ARGIND==2 && (!($1 in PreviousType)){print $1}'\
                            $StatisticalFold"Distribution/previoustype"\
                            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation | sort | uniq | wc -l 
#Region distribution
Types=`cut -f10 /data9/swu13/CAeditome_hg38/GeneAnnotaion/RegionInformation | tail -n +2 | sort | uniq`
Allnum=`echo $Types | awk -F" " '{print NF}'`
rm -f $StatisticalFold"Distribution/Regions"
for((i=1;i<=$Allnum;i++))
do
type=`echo $Types |cut -d" " -f$i`
num=`awk -F"\t" -v type="$type" '($10==type){print $1}' /data9/swu13/CAeditome_hg38/GeneAnnotaion/RegionInformation | sort | uniq | wc -l`
echo -e $type"\t"$num >> $StatisticalFold"Distribution/Regions"
done
#Repeat distribution
Types=`cut -f10 /data9/swu13/CAeditome_hg38/RepeatAnnotation/RepeatInformation | tail -n +2 | sort | uniq`
Types=("Alu" "NoRepeat" "L1" "ERVL-MaLR")
Allnum=${#Types[@]}
rm -f $StatisticalFold"Distribution/Repeats"
for((i=0;i<$Allnum;i++))
do
type=${Types[$i]}
num=`awk -F"\t" -v type="$type" '($10==type){print $1}' /data9/swu13/CAeditome_hg38/RepeatAnnotation/RepeatInformation | sort | uniq | wc -l`
echo -e $type"\t"$num >> $StatisticalFold"Distribution/Repeats"
done
awk -F"\t" -v Types="${Types[*]}" '(NR>1){split(Types,a," ");for(i=1;i<=length(a);i++){type[a[i]]=a[i]};if(!($10 in type)){print $1}}' /data9/swu13/CAeditome_hg38/RepeatAnnotation/RepeatInformation | sort | uniq | wc -l
}

#7. Splicing analysis
function SplicingAnalysis(){
mkdir $StatisticalFold"SplicingAnalysis/"
mv /data9/swu13/CAeditome_hg38/Splicing/Files/SplicingInformation /data9/swu13/CAeditome_hg38/Splicing/Files/SplicingInformation.previous
head -n 1 /data9/swu13/CAeditome_hg38/Splicing/Files/SplicingInformation.previous > /data9/swu13/CAeditome_hg38/Splicing/Files/SplicingInformation
awk -F"\t" 'ARGIND==1 && (FNR>1){split($3,a,".");if(!(Gene[$1] ~ a[1])){Gene[$1]=Gene[$1]"|"a[1]}}\
            ARGIND==2 && (FNR>1){split(Gene[$1],a,"|");split($8,b,".");for(i=2;i<=length(a);i++){if(a[i]==b[1]){print $0}}}'\
            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
            /data9/swu13/CAeditome_hg38/Splicing/Files/SplicingInformation.previous >> /data9/swu13/CAeditome_hg38/Splicing/Files/SplicingInformation
awk -F"\t" '(NR>1){print $13"_"$14"\t"$1"\t"$16-$12}' /data9/swu13/CAeditome_hg38/Splicing/Files/SplicingInformation  |sort | uniq > $StatisticalFold"SplicingAnalysis/Strength"
awk -F"\t" '($1~"3SS_"){print $2}' $StatisticalFold"SplicingAnalysis/Strength" | sort | uniq | wc -l
awk -F"\t" '($1~"5SS_"){print $2}' $StatisticalFold"SplicingAnalysis/Strength" | sort | uniq | wc -l
awk -F"\t" '(NR>1){split($8,a,".");print a[1]}' /data9/swu13/CAeditome_hg38/Splicing/Files/SplicingInformation | sort | uniq | wc -l
Rscript Boxplot.R $StatisticalFold"SplicingAnalysis/Strength" $StatisticalFold"SplicingAnalysis/Strength.png" 500 300 "TRUE"
rm -f $StatisticalFold"SplicingAnalysis/PSI"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
if [[ -f "/data9/swu13/CAeditome_hg38/Splicing/Files/"$cancer".WT.RNAediting.PSI" ]]
then
awk -v cancer="${cancer}" -F"\t" 'ARGIND==1{CAeditomeID[$2]=$1;if(!(Gene[$2]~$3"|"$5)){Gene[$2]=Gene[$2]";"$3"|"$5}}\
                                  ARGIND==2{Tumor[$1]=$1}\
                                  ARGIND==3 && ($5>=25) && ($7-$5>=25) && (FNR>1){Editing[$1]=$1}
                                  ARGIND==4 && ($16<0.05 && $20<0.05) && ($19*$21>0) && (FNR>1) && ($1 in Editing){split(Gene[$1],a,";");for(i=2;i<=length(a);i++){split(a[i],b,"|");\
                                                     out=cancer"\t"CAeditomeID[$1]"\t"$1"\t"b[1]"\t"b[2]"\t"$2"\t"$12"_"$13"\t"$16"\t"$19"\t"$20"\t"$21;\
                                                     if(b[2] in Tumor){out=out"\tY"}else{out=out"\tN"};print out}}'\
                                  /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
                                  /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
                                  "/data9/swu13/CAeditome_hg38/Differential/Files/"$cancer".DifferentialResults"\
                                  "/data9/swu13/CAeditome_hg38/Splicing/Files/"$cancer".WT.RNAediting.PSI" >>  $StatisticalFold"SplicingAnalysis/PSI"
fi
done
awk -F"\t" 'ARGIND==1{Known[$1"\t"$3]=$1"\t"$3}\
            ARGIND==2 && ($2"\t"$6 in Known){print $0}'\
            /data9/swu13/CAeditome_hg38/Splicing/Files/SplicingInformation\
            $StatisticalFold"SplicingAnalysis/PSI" > $StatisticalFold"SplicingAnalysis/PSI.1"
mv $StatisticalFold"SplicingAnalysis/PSI.1" $StatisticalFold"SplicingAnalysis/PSI"
#All
awk -F"\t" '($9*$11>0) && ($7=="3SS_2i" || $7=="5SS_2e"){print $2}' $StatisticalFold"SplicingAnalysis/PSI" | sort | uniq | wc -l
awk -F"\t" '($9*$11>0) && ($7=="3SS_2i" || $7=="5SS_2e"){print $5}' $StatisticalFold"SplicingAnalysis/PSI" | sort | uniq 
awk -F"\t" '($9*$11>0) && ($7=="3SS_2i" || $7=="5SS_2e") && ($12=="Y"){print $2}' $StatisticalFold"SplicingAnalysis/PSI" | sort | uniq | wc -l
#Statistic
rm -f $StatisticalFold"SplicingAnalysis/Statistical"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
num1=`awk -v cancer="${cancer}" -F"\t" '($1==cancer){print $2}' $StatisticalFold"SplicingAnalysis/PSI" |  sort | uniq | wc -l`
num2=`awk -v cancer="${cancer}" -F"\t" '($1==cancer){print $4}' $StatisticalFold"SplicingAnalysis/PSI" |  sort | uniq | wc -l`
num3=`awk -v cancer="${cancer}" -F"\t" '($1==cancer) && ($12=="Y"){print $2}' $StatisticalFold"SplicingAnalysis/PSI" |  sort | uniq | wc -l`
num4=`awk -v cancer="${cancer}" -F"\t" '($1==cancer) && ($12=="Y"){print $4}' $StatisticalFold"SplicingAnalysis/PSI" |  sort | uniq | wc -l`
echo -e $cancer"\t"$num1"\t"$num2"\t"$num3"\t"$num4 >> $StatisticalFold"SplicingAnalysis/Statistical"
done
awk -F"\t" '{print $2}' $StatisticalFold"SplicingAnalysis/PSI" |  sort | uniq | wc -l
awk -F"\t" '{print $4}' $StatisticalFold"SplicingAnalysis/PSI" |  sort | uniq | wc -l
awk -F"\t" '($12=="Y"){print $2}' $StatisticalFold"SplicingAnalysis/PSI" |  sort | uniq | wc -l
awk -F"\t" '($12=="Y"){print $4}' $StatisticalFold"SplicingAnalysis/PSI" |  sort | uniq | wc -l
#Selected for database
#awk -F"\t" '($9<-0.1 || $9>0.1){print "/data9/swu13/CAeditome_hg38/Splicing/Boxplot/"$1"/"$3"."$6".png"}' $StatisticalFold"SplicingAnalysis/PSI" | sort | uniq  > $StatisticalFold"SplicingAnalysis/SelectedBoxplot"
#awk -F"\t" '($11<-0.2 || $11>0.2){print "/data9/swu13/CAeditome_hg38/Splicing/Correlation/"$1"/"$3"."$6".png"}' $StatisticalFold"SplicingAnalysis/PSI" | sort | uniq  > $StatisticalFold"SplicingAnalysis/SelectedCorrelation"
##match ENSG different version to v22
#echo -e "SplicingENSG\tEditingENSG\tGeneName" > $StatisticalFold"SplicingAnalysis/SplicingENSGMatch"
#awk -F"\t" 'ARGIND==1 && (FNR>1){split($3,a,".");ENSG[a[1]]=$3"\t"$5;if(!(Editing[$2]~$3)){Editing[$2]=Editing[$2]";"$3"|"$5}}\
#            ARGIND==2 && (FNR>1){split($8,a,".");if(a[1] in ENSG){print $8"\t"ENSG[a[1]]}\
#                                                 else{split(Editing[$2],m,";");for(i=2;i<=length(m);i++){split(m[i],n,"|");print $8"\t"n[1]"\t"n[2]}}}'\
#           /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
#           /data9/swu13/CAeditome_hg38/Splicing/Files/SplicingInformation  | sort | uniq >> $StatisticalFold"SplicingAnalysis/SplicingENSGMatch"
}

#8. miRNA bind
function miRNABindAnalysis(){
mkdir $StatisticalFold"miRNA/"
#split files for Database
Original=("EditedmRNA|Gain" "EditedmRNA|Loss" "EditedLncRNA|Gain" "EditedLncRNA|Loss" "EditedmiRNA_BindmRNA|Gain" "EditedmiRNA_BindmRNA|Loss" "EditedmiRNA_BindlncRNA|Gain" "EditedmiRNA_BindlncRNA|Loss")
Present=("mRNA3UTR_miRNA_gain.txt" "mRNA3UTR_miRNA_loss.txt" "lncRNA_miRNA_gain.txt" "lncRNA_miRNA_loss.txt" "miRNA_mRNA3UTR_gain.txt" "miRNA_mRNA3UTR_loss.txt" "miRNA_lncRNA_gain.txt" "miRNA_lncRNA_loss.txt")
for((i=0;i<8;i++))
do
Edtype=`echo ${Original[$i]} | cut -d"|" -f1`
GainLoss=`echo ${Original[$i]} | cut -d"|" -f2`
file=${Present[$i]}
awk -F"\t" -v Edtype="$Edtype" -v GainLoss="$GainLoss" '(($14==Edtype) && ($15==GainLoss)) || (NR==1)\
           {out=$1"\t"$2"\t"$4"\t"$5"\t"$3;for(i=6;i<=13;i++){out=out"\t"$i};print out}'\
    /data9/swu13/CAeditome_hg38/miRNA/Files/RNAeditingMiRNABindingInformation > "/data9/swu13/CAeditome_hg38/miRNA/Files/"$file
done
#statistics
awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
            ARGIND==2 && (FNR>1){split($3,a,".");split($4,b,".");\
                                 out=$1"\t"$2"\t"a[1]"\t"b[1]"\t"$5"\t"$6"\t"$14"\t"$15;if($5 in Tumor){out=out"\tY"}else{out=out"\tN"};print out}'\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            /data9/swu13/CAeditome_hg38/miRNA/Files/RNAeditingMiRNABindingInformation > $StatisticalFold"miRNA/miRNABind"            
#All statistic
awk -F"\t" '{print $1}' $StatisticalFold"miRNA/miRNABind" | sort | uniq | wc -l
awk -F"\t" '(!($7~"EditedmiRNA")){print $3}' $StatisticalFold"miRNA/miRNABind" | sort | uniq | wc -l
awk -F"\t" '(($7~"EditedmiRNA")){print $6}' $StatisticalFold"miRNA/miRNABind" | sort | uniq | wc -l
awk -F"\t" '($8=="Gain"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l
awk -F"\t" '($8=="Loss"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l
awk -F"\t" '($7=="EditedLncRNA"){print $3}' $StatisticalFold"miRNA/miRNABind" | sort | uniq | wc -l > $StatisticalFold"miRNA/Statistic"   
awk -F"\t" '($7=="EditedLncRNA"){print $1}' $StatisticalFold"miRNA/miRNABind" | sort | uniq | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7=="EditedmRNA"){print $3}' $StatisticalFold"miRNA/miRNABind" | sort | uniq | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7=="EditedmRNA"){print $1}' $StatisticalFold"miRNA/miRNABind" | sort | uniq | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7~"EditedmiRNA"){print $6}' $StatisticalFold"miRNA/miRNABind" | sort | uniq | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7~"EditedmiRNA"){print $1}' $StatisticalFold"miRNA/miRNABind" | sort | uniq | wc -l >> $StatisticalFold"miRNA/Statistic"
#
awk -F"\t" '($7=="EditedLncRNA") && ($9=="Y"){print $3}' $StatisticalFold"miRNA/miRNABind" | sort | uniq | wc -l > $StatisticalFold"miRNA/Statistic"   
awk -F"\t" '($7=="EditedLncRNA") && ($9=="Y"){print $1}' $StatisticalFold"miRNA/miRNABind" | sort | uniq | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7=="EditedmRNA") && ($9=="Y"){print $3}' $StatisticalFold"miRNA/miRNABind" | sort | uniq | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7=="EditedmRNA") && ($9=="Y"){print $1}' $StatisticalFold"miRNA/miRNABind" | sort | uniq | wc -l >> $StatisticalFold"miRNA/Statistic"
#
awk -F"\t" '($7=="EditedLncRNA") && ($8=="Gain"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l > $StatisticalFold"miRNA/Statistic" 
awk -F"\t" '($7=="EditedLncRNA") && ($8=="Loss"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7=="EditedLncRNA") && ($8=="Gain") && ($9=="Y"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7=="EditedLncRNA") && ($8=="Loss") && ($9=="Y"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l >> $StatisticalFold"miRNA/Statistic"
#
awk -F"\t" '($7=="EditedmRNA") && ($8=="Gain"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l > $StatisticalFold"miRNA/Statistic" 
awk -F"\t" '($7=="EditedmRNA") && ($8=="Loss"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7=="EditedmRNA") && ($8=="Gain") && ($9=="Y"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7=="EditedmRNA") && ($8=="Loss") && ($9=="Y"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l >> $StatisticalFold"miRNA/Statistic"
#
awk -F"\t" '($7=="EditedmiRNA_BindlncRNA") && ($8=="Gain"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l > $StatisticalFold"miRNA/Statistic" 
awk -F"\t" '($7=="EditedmiRNA_BindlncRNA") && ($8=="Loss"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7=="EditedmiRNA_BindlncRNA") && ($8=="Gain") && ($9=="Y"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7=="EditedmiRNA_BindlncRNA") && ($8=="Loss") && ($9=="Y"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l >> $StatisticalFold"miRNA/Statistic"
#
awk -F"\t" '($7=="EditedmiRNA_BindmRNA") && ($8=="Gain"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l > $StatisticalFold"miRNA/Statistic" 
awk -F"\t" '($7=="EditedmiRNA_BindmRNA") && ($8=="Loss"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7=="EditedmiRNA_BindmRNA") && ($8=="Gain") && ($9=="Y"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l >> $StatisticalFold"miRNA/Statistic"
awk -F"\t" '($7=="EditedmiRNA_BindmRNA") && ($8=="Loss") && ($9=="Y"){print $0}' $StatisticalFold"miRNA/miRNABind" | wc -l >> $StatisticalFold"miRNA/Statistic"
#Each Cancer
rm -f $StatisticalFold"miRNA/miRNABind.Cancers"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
awk -v cancer="${cancer}" -F"\t" 'ARGIND==1{Editing[$1]=$1}\
                                  ARGIND==2{Tumor[$1]=$1}\
                                  ARGIND==3 && (FNR>1) && ($2 in Editing){out=cancer"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$14"\t"$15;\
                                                     if($5 in Tumor){out=out"\tY"}else{out=out"\tN"};print out}'\
                                  "/data9/swu13/CAeditome_hg38/GeneralAnnotation/"$cancer".Table"\
                                  /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
                                  /data9/swu13/CAeditome_hg38/miRNA/Files/RNAeditingMiRNABindingInformation >>  $StatisticalFold"miRNA/miRNABind.Cancers"
done
#Statistical
rm -f $StatisticalFold"miRNA/Statistic"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
num1=`awk -v cancer="${cancer}" -F"\t" '($1==cancer){print $2}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l`
num20=`awk -v cancer="${cancer}" -F"\t" '($1==cancer) && (!($8~"EditedmiRNA")){print $5}' $StatisticalFold"miRNA/miRNABind.Cancers"  | sort | uniq |wc -l`
num21=`awk -v cancer="${cancer}" -F"\t" '($1==cancer) && (($8~"EditedmiRNA")){print $7}' $StatisticalFold"miRNA/miRNABind.Cancers"  | sort | uniq |wc -l`
num3=`awk -v cancer="${cancer}" -F"\t" '($1==cancer) && ($10=="Y"){print $2}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l`
num4=`awk -v cancer="${cancer}" -F"\t" '($1==cancer) && ($10=="Y") && (!($8~"EditedmiRNA")){print $5}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l`
num2=`expr $num20 + $num21`
echo -e $cancer"\t"$num1"\t"$num2"\t"$num3"\t"$num4"\t"$num20"\t"$num21 >> $StatisticalFold"miRNA/Statistic"
done
#All
awk -F"\t" '{print $2}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l
awk -F"\t" '(!($8~"EditedmiRNA")){print $5}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l
awk -F"\t" '(($8~"EditedmiRNA")){print $7}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l
awk -F"\t" '($10=="Y"){print $2}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l
awk -F"\t" '($10=="Y") && (!($8~"EditedmiRNA")){print $5}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l
#Statistical
rm -f $StatisticalFold"miRNA/Statistic"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
num1=`awk -v cancer="${cancer}" -F"\t" '($1==cancer){print $2}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l`
num2=`awk -v cancer="${cancer}" -F"\t" '($1==cancer){print $5}' $StatisticalFold"miRNA/miRNABind.Cancers"  | sort | uniq |wc -l`
num3=`awk -v cancer="${cancer}" -F"\t" '($1==cancer) && ($10=="Y"){print $2}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l`
num4=`awk -v cancer="${cancer}" -F"\t" '($1==cancer) && ($10=="Y"){print $5}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l`
echo -e $cancer"\t"$num1"\t"$num2"\t"$num3"\t"$num4 >> $StatisticalFold"miRNA/Statistic"
done
#All
awk -F"\t" '{print $2}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l
awk -F"\t" '{print $5}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l
awk -F"\t" '($10=="Y"){print $2}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l
awk -F"\t" '($10=="Y"){print $5}' $StatisticalFold"miRNA/miRNABind.Cancers" | sort | uniq |wc -l
}

#9. miRNA regulation
function miRNARegulationAnalysis(){
mkdir $StatisticalFold"miRNA/"
awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
            ARGIND==2 && (FNR>1){out=$0;if($5 in Tumor){out=out"\tY"}else{out=out"\tN"};\
                                 if($10 in Tumor){out=out"\tY"}else{out=out"\tN"};print out}'\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            /data9/swu13/CAeditome_hg38/miRNA/Files/miRNARegulation > $StatisticalFold"miRNA/miRNARegulation"
#All Statistic
cut -f2 $StatisticalFold"miRNA/miRNARegulation" | sort | uniq | wc -l
awk -F"\t" '{split($4,a,".");split($9,b,".");if($6~"miRNA"){print a[1]"\n"b[1]}else{print b[1]}}' $StatisticalFold"miRNA/miRNARegulation" | sort | uniq | wc -l
awk -F"\t" '{split($4,a,".");split($9,b,".");if($6~"miRNA"){if($25=="Y" || $26=="Y"){print $2}}\
                                             else{if($26=="Y"){print $2}}}' $StatisticalFold"miRNA/miRNARegulation" | sort | uniq | wc -l
awk -F"\t" '{split($4,a,".");split($9,b,".");if($6~"miRNA"){if($25=="Y"){print a[1]};if($26=="Y"){print b[1]}}\
                                             else{if($26=="Y"){print b[1]}}}' $StatisticalFold"miRNA/miRNARegulation" | sort | uniq | wc -l
#Statistic
rm -f $StatisticalFold"miRNA/Statistic"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
num1=`awk -F"\t" -v cancer="$cancer" '($1==cancer){print $2}' $StatisticalFold"miRNA/miRNARegulation" | sort | uniq | wc -l`
num2=`awk -F"\t" -v cancer="$cancer" '($1==cancer){split($4,a,".");split($9,b,".");if($6~"miRNA"){print a[1]"\n"b[1]}else{print b[1]}}'\
                                        $StatisticalFold"miRNA/miRNARegulation" | sort | uniq | wc -l`
num3=`awk -F"\t" -v cancer="$cancer" '($1==cancer){split($4,a,".");split($9,b,".");if($6~"miRNA"){if($25=="Y" || $26=="Y"){print $2}}\
                                                   else{if($26=="Y"){print $2}}}' $StatisticalFold"miRNA/miRNARegulation" | sort | uniq | wc -l`
num4=`awk -F"\t" -v cancer="$cancer" '($1==cancer){split($4,a,".");split($9,b,".");if($6~"miRNA"){if($25=="Y"){print a[1]};if($26=="Y"){print b[1]}}\
                                                   else{if($26=="Y"){print b[1]}}}' $StatisticalFold"miRNA/miRNARegulation" | sort | uniq | wc -l`
echo -e $cancer"\t"$num1"\t"$num2"\t"$num3"\t"$num4 >> $StatisticalFold"miRNA/Statistic"
done
awk -F"\t" '(NR>1){split($4,a,".");split($9,b,".");if($6~"miRNA"){if(!(a[1]"\t"$8"\t"b[1] in MiRNA) && !(b[1]"\t"$8"\t"a[1] in MiRNA)){\
                                                                        MiRNA[a[1]"\t"$8"\t"b[1]]=a[1]"\t"$8"\t"b[1];print a[1]"\t"$8"\t"b[1]}}\
                                                   else{print a[1]"\t"$8"\t"b[1]}}' $StatisticalFold"miRNA/miRNARegulation" | sort | uniq | wc -l
} 

#10. Structure
function StructureAnalysis(){
mkdir $StatisticalFold"Structure/"
awk -F"\t" 'ARGIND==1 && (substr($1,1,1)==">"){Without[substr($1,2)]=substr($1,2)}\
            ARGIND==2 && (substr($1,1,1)==">") && (substr($1,2) in Without)\
              {out="/data9/swu13/CAeditome_hg38/Structure/RNAstructure/withoutEditingResizePng/"substr($1,2)".png";\
               print out"\t/data9/swu13/CAeditome_hg38/Structure/RNAstructure/withmultipleEditingPng/"substr($1,2)".png"}'\
               /data9/swu13/CAeditome_hg38/Structure/Files/WithoutEditing.fold\
               /data9/swu13/CAeditome_hg38/Structure/Files/WithMultipleEditing.fold > $StatisticalFold"Structure/SelectedStructure"
#Boxplot
awk -F"\t" '(NR>1){print "WideType\t.\t"$5;print "EditedType\t.\t"$4}' /data9/swu13/CAeditome_hg38/Structure/Files/RNAstructureInformation > $StatisticalFold"Structure/Structure.boxplot.data"
Rscript Boxplot.R $StatisticalFold"Structure/Structure.boxplot.data" $StatisticalFold"Structure/Structure.boxplot.png" 200 300 "FALSE"
awk -F"\t" '(NR>1){print "EditedType-WideType\t.\t"$4-$5}' /data9/swu13/CAeditome_hg38/Structure/Files/RNAstructureInformation > $StatisticalFold"Structure/Structure.boxplot.data"
Rscript Boxplot.R $StatisticalFold"Structure/Structure.boxplot.data" $StatisticalFold"Structure/Structure.boxplot.png" 200 300 "FALSE" -2 1
Rscript Violin.R $StatisticalFold"Structure/Structure.boxplot.data" $StatisticalFold"Structure/Structure.violin.png" 200 300
#Statistic
awk -F"\t" '($6<$5) && (NR>1){print $2}' /data9/swu13/CAeditome_hg38/Structure/Files/RNAstructureInformation | sort | uniq | wc -l >>  $StatisticalFold"Structure/Statistic"
cut -f2 /data9/swu13/CAeditome_hg38/Structure/Files/RNAstructureInformation | sort | uniq | wc -l | awk '{print $0-1}' >>  $StatisticalFold"Structure/Statistic"
#gene
awk -F"\t" 'ARGIND==1 && (FNR>1){ENSG[$6]=$3}\
            ARGIND==2 && (FNR>1){print ENSG[$2]}'\
            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
            /data9/swu13/CAeditome_hg38/Structure/Files/RNAstructureInformation | sort | uniq | wc -l
awk -F"\t" 'ARGIND==1 && (FNR>1){ENSG[$6]=$3}\
            ARGIND==2 && (FNR>1) && ($6<$5){print ENSG[$2]}'\
            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
            /data9/swu13/CAeditome_hg38/Structure/Files/RNAstructureInformation | sort | uniq | wc -l
awk -F"\t" 'ARGIND==1 && (FNR>1){Gene[$6]=$5}\
            ARGIND==2 {Tumor[$1]=$1}\
            ARGIND==3 && (FNR>1) && ($4<$5) && (Gene[$2] in Tumor){print $1}'\
            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            /data9/swu13/CAeditome_hg38/Structure/Files/RNAstructureInformation | sort | uniq | wc -l
}

#11. ADAR correlation
function ADARanalysis(){
mkdir $StatisticalFold"ADAR/"
awk -F"\t" 'ARGIND==1 && (FNR>1){Gene[$3]=$5}\
            ARGIND==2{Tumor[$1]=$1}\
            ARGIND==3 && (($5<0.05 && $6>0.3) || ($7<0.05 && $8>0.3) || ($9<0.05 && $10>0.3))\
            {out=$0"\t"Gene[$4];if(Gene[$4 in Tumor]){print out"\tY"}else{print out"\tN"}}'\
            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            /data9/swu13/CAeditome_hg38/ADARCorrelation/Files/ADARInformation > $StatisticalFold"ADAR/ADARInformation"
awk -F"\t" '{if($5<0.05 && $6>0.3){print $2"\t.\tADAR1"};\
             if($7<0.05 && $8>0.3){print $2"\t.\tADAR2"};\
             if($9<0.05 && $10>0.3){print $2"\t.\tADAR3"}}'\
            $StatisticalFold"ADAR/ADARInformation" | sort | uniq > $StatisticalFold"ADAR/vennplot.data"
Rscript VennPlot.R $StatisticalFold"ADAR/vennplot.data" $StatisticalFold"ADAR/vennplot.data.png" "FALSE" 
#
rm -f $StatisticalFold"ADAR/ADARInformation.none"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
ExpressionFile="/data9/swu13/CAeditome_hg38/ADARCorrelation/Files/"$cancer".ADAR.Results"
DifferentialFile="/data9/swu13/CAeditome_hg38/Differential/Files/"$cancer".DifferentialResults"
awk -v cancer="$cancer" 'ARGIND==1 && (FNR>1) && ($5>=50){EditingNum[$1]=$1}\
                         ARGIND==2 && (FNR>1){Known[$3]=$3}
                         ARGIND==3 && (FNR>1) && ($1 in EditingNum) && (($7=="NA") || ($7>=0.05)) && (!($1 in Known)){\
                                          All[cancer"\t"$1]=cancer"\t"$1;\
                                          if($2=="ENSG00000160710"){ADAR1[cancer"\t"$1]=$7"\t"$8}\
                                          else if($2=="ENSG00000197381"){ADAR2[cancer"\t"$1]=$7"\t"$8}\
                                          else if($2=="ENSG00000185736"){ADAR3[cancer"\t"$1]=$7"\t"$8}}\
                                   END{for(i in All){out=i;if((i in ADAR1) && (i in ADAR2) && (i in ADAR3)){\
                                                               out=out"\t"ADAR1[i]"\t"ADAR2[i]"\t"ADAR3[i];;print out}}}'\
                         $DifferentialFile\
                         /data9/swu13/CAeditome_hg38/ADARCorrelation/Files/ADARInformation\
                         $ExpressionFile >> $StatisticalFold"ADAR/ADARInformation.none"
done
awk -F"\t" 'ARGIND==1{Known[$3]=$3}\
            ARGIND==2 && (!($2 in Known)){print $2}'\
            $StatisticalFold"ADAR/ADARInformation"\
            $StatisticalFold"ADAR/ADARInformation.none" |sort| uniq | wc -l
#
rm -f $StatisticalFold"ADAR/Statistic"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
num1=`awk -F"\t" -v cancer="$cancer" '($1==cancer) && ($5<0.05) && ($6>0.3){print $2}'\
                 $StatisticalFold"ADAR/ADARInformation" | sort | uniq | wc -l`
num2=`awk -F"\t" -v cancer="$cancer" '($1==cancer) && ($7<0.05) && ($8>0.3){print $2}'\
                 $StatisticalFold"ADAR/ADARInformation" | sort | uniq | wc -l`
num3=`awk -F"\t" -v cancer="$cancer" '($1==cancer) && ($9<0.05) && ($10>0.3){print $2}'\
                 $StatisticalFold"ADAR/ADARInformation" | sort | uniq | wc -l`
echo -e $cancer"\t"$num1"\t"$num2"\t"$num3 >> $StatisticalFold"ADAR/Statistic"
done         
awk -F"\t" '($5<0.05) && ($6>0.3){print $2}' $StatisticalFold"ADAR/ADARInformation" | sort | uniq | wc -l
awk -F"\t" '($7<0.05) && ($8>0.3){print $2}' $StatisticalFold"ADAR/ADARInformation" | sort | uniq | wc -l
awk -F"\t" '($9<0.05) && ($10>0.3){print $2}' $StatisticalFold"ADAR/ADARInformation" | sort | uniq | wc -l
awk -F"\t" '{print $4}' $StatisticalFold"ADAR/ADARInformation" | sort | uniq | wc -l
awk -F"\t" 'ARGIND==1 && (FNR>1){Gene[$3]=$5}\
            ARGIND==2{Tumor[$1]=$1}\
            ARGIND==3 &&(Gene[$4] in Tumor){print $4}'\
            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            $StatisticalFold"ADAR/ADARInformation" | sort | uniq | wc -l
}

#12. drug and disease analysis
function DrugDiseaseAnalysis(){
mkdir $StatisticalFold"DrugDisease/"
#Disease
awk -F"\t" 'ARGIND==1 && (FNR>1){Gene[$5]=$5}\
            ARGIND==2 && ($2 in Gene){print $2"\t"$3}'\
            ../GeneralAnnotation/GeneralInformation\
            ../RawData/DrugDisease/caeditome_genes_disease.txt | sort | uniq > $StatisticalFold"DrugDisease/Disease"
cut -f1 $StatisticalFold"DrugDisease/Disease" | sort | uniq | wc -l
cut -f2 $StatisticalFold"DrugDisease/Disease" | sort | uniq | wc -l
#Drug
awk -F"\t" 'ARGIND==1 && (FNR>1){Protein[$15]=$6}\
            ARGIND==2 && ($5 in Protein){print Protein[$5]"\t"$2}'\
            ../ProteinAnnotation/ProteinLolliplot\
            ../RawData/DrugDisease/caeditome_drug.txt | sort | uniq > $StatisticalFold"DrugDisease/Drug"
cut -f1 $StatisticalFold"DrugDisease/Drug" | sort | uniq | wc -l
cut -f2 $StatisticalFold"DrugDisease/Drug" | sort | uniq | wc -l
}


