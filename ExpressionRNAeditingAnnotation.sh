cancertypes=("BLCA"  "LGG"  "BRCA"  "LIHC"  "CESC"  "LUAD"  "COAD"  "LUSC"  "GBM"  "PRAD"  "HNSC"  "STAD"  "KICH"  "THCA"  "KIRC"  "UCEC"  "KIRP"  "ACC"  "DLBC"  "ESCA"  "LAML"  "OV"  "PAAD"  "PCPG"  "SKCM"  "UCS"  "UVM"  "CHOL"  "MESO"  "READ"  "SARC"  "TGCT"  "THYM")
DatabaseFold="/data9/swu13/CAeditome_hg38/"

function MkdirPath(){
mkdir $DatabaseFold"ADARCorrelation/"
mkdir $DatabaseFold"ExpressionCorrelation/"
mkdir $DatabaseFold"DifferentialExpression/"
mkdir $DatabaseFold"StageExpressionAssociation/"
mkdir $DatabaseFold"SurvivalExpressionAssociation/"
}

function ExpressionAnnotation(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
DownloadRNAeditingFold=$DownloadFold"RNAediting/"
DownloadExpressionFold=$DownloadFold"GeneExpression/"
RNAeditingFile=$DownloadRNAeditingFold""$cancer"_Editing.txt"
ExpressionFile=$DownloadExpressionFold""$cancer"_expression.csv"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
ExpressionFold=$DatabaseFold"ExpressionCorrelation/"
ExpressionFileFold=$ExpressionFold"Files/"
ExpressionTmpFold=$ExpressionFold"Tmp/"
ExpressionCorrelationFold=$ExpressionFold"Correlation/"
ExpressionBoxplotFold=$ExpressionFold"Boxplot/"
mkdir $ExpressionFileFold
mkdir $ExpressionTmpFold
mkdir $ExpressionCorrelationFold
mkdir $ExpressionBoxplotFold
head -n 1 $RNAeditingFile | tr '\t' '\n' | awk '(NR>1){print $1"\t"NR}' > $ExpressionTmpFold""$cancer".RNAediting.samples"
head -n 1 $ExpressionFile | tr ',' '\n' | awk '(NR>1){print $1"\t"NR}' > $ExpressionTmpFold""$cancer".Expression.samples"
awk -F"\t" 'ARGIND==1{split($1,a,"-");if(substr(a[4],1,1)=="0"){ENSG[a[1]"-"a[2]"-"a[3]"-"a[4]]=ENSG[a[1]"-"a[2]"-"a[3]"-"a[4]]";"$2}}\
            ARGIND==2{if($1 in ENSG){print $0"\t"$1"\t"substr(ENSG[$1],2)}}'\
            $ExpressionTmpFold""$cancer".Expression.samples"\
            $ExpressionTmpFold""$cancer".RNAediting.samples" > $ExpressionTmpFold""$cancer".RNAediting.Expression.tumor.samples"
RNAeditingColumn=`cut -f2 $ExpressionTmpFold""$cancer".RNAediting.Expression.tumor.samples"`
ExpressionColumn=`cut -f4 $ExpressionTmpFold""$cancer".RNAediting.Expression.tumor.samples"`
awk -F"\t" -v RNAeditingColumn="${RNAeditingColumn[*]}" 'ARGIND==1{split($2,a,".");ENSG[$1]=a[1]}\
            ARGIND==2{split($0,a,",");ENSGexpression[a[1]]=a[1]}\
            ARGIND==3 && ($1 in ENSG) && (ENSG[$1] in ENSGexpression){split(RNAeditingColumn,col," ");out=$1;for(i=1;i<=length(col);i++){out=out"\t"$col[i]};print out}'\
            $GeneAnnotationFold""$cancer".Table"\
            $ExpressionFile\
            $RNAeditingFile > $ExpressionTmpFold""$cancer".RNAediting.tumor"
awk -F"\t" -v ExpressionColumn="${ExpressionColumn[*]}" 'ARGIND==1{split($2,a,".");ENSG[$1]=a[1]}\
            ARGIND==2{split($0,a,",");split(ExpressionColumn,col," ");out=a[1];for(i=1;i<=length(col);i++){split(col[i],b,";");m=0;for(j=1;j<=length(b);j++){m=m+a[b[j]]};out=out"\t"m/length(b)};ENSG_value[a[1]]=out}\
            ARGIND==3{print ENSG_value[ENSG[$1]]}'\
            $GeneAnnotationFold""$cancer".Table"\
            $ExpressionFile\
            $ExpressionTmpFold""$cancer".RNAediting.tumor" > $ExpressionTmpFold""$cancer".Expression.tumor"
mkdir $ExpressionBoxplotFold""$cancer"/"
mkdir $ExpressionCorrelationFold""$cancer"/"
Rscript ExpressionCorrelation.R $ExpressionTmpFold""$cancer".RNAediting.tumor" $ExpressionTmpFold""$cancer".Expression.tumor" $ExpressionBoxplotFold""$cancer"/" $ExpressionCorrelationFold""$cancer"/" $ExpressionFileFold""$cancer".Results"
}


function ExpressionHeatmap(){
ExpressionFold=$DatabaseFold"ExpressionCorrelation/"
ExpressionHeatmapFold=$ExpressionFold"Heatmap/"
ExpressionTmpFold=$ExpressionFold"Tmp/"
mkdir $ExpressionHeatmapFold
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
rm -f $ExpressionTmpFold"Heatmap.data"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
ExpressionFile=$ExpressionFold"Files/"$cancer".Results"
DifferentialFile=$DatabaseFold"Differential/Files/"$cancer".DifferentialResults"
awk -v cancer="$cancer" 'ARGIND==1{Editing[$2]=$1}\
                         ARGIND==2 && (FNR>1) && ($5>=25) && ($7-$5>=25){EditedNum[$1]=$1}
                         ARGIND==3 && (FNR>1) && ($1 in EditedNum) &&($3!="NA") && ($3<0.05) && ($7!="NA") && ($7<0.05){\
                                       print cancer"\t"Editing[$1]"\t"$1"\t"$2"\t"$6}'\
                         $GeneAnnotationFold"GeneralInformation"\
                         $DifferentialFile\
                         $ExpressionFile >> $ExpressionTmpFold"Heatmap.lofFC.data"
awk -v cancer="$cancer" 'ARGIND==1{Editing[$2]=$1}\
                         ARGIND==2 && (FNR>1) && ($5>=25){EditedNum[$1]=$1}
                         ARGIND==3 && (FNR>1) && ($1 in EditedNum) &&($3!="NA") && ($3<0.05) && ($7!="NA") && ($7<0.05){\
                                       print cancer"\t"Editing[$1]"\t"$1"\t"$2"\t"$8}'\
                         $GeneAnnotationFold"GeneralInformation"\
                         $DifferentialFile\
                         $ExpressionFile >> $ExpressionTmpFold"Heatmap.R.data"
done
mkdir $ExpressionHeatmapFold"/LogFC/"
Rscript Heatmap.R $ExpressionTmpFold"Heatmap.lofFC.data" 1 4 2 5 0 99 $ExpressionHeatmapFold"/LogFC/" #cancer\tGene\tEditing\tvalue
mkdir $ExpressionHeatmapFold"/CorR/"
Rscript Heatmap.R $ExpressionTmpFold"Heatmap.R.data" 1 4 2 5 0 99 $ExpressionHeatmapFold"/CorR/" #cancer\tGene\tEditing\tvalue
}

function MergeExpression(){
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
ExpressionFold=$DatabaseFold"ExpressionCorrelation/"
echo -e "Cancer\tCAeditomeID\tEditingInformation\tGene\tGeneName\tTtestP\tMeanExpressionEdited\tMeanExpressionNonedited\tlogFC\tCorreP\tCorreR" > $ExpressionFold"Files/RNAeditingExpressionInformation"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && (FNR>1){CAeditomeID[$2]=$1;split($3,a,".");GeneName[a[1]]=$5}\
                                ARGIND==2 && (FNR>1){out=cancer"\t"CAeditomeID[$1]"\t"$1"\t"$2"\t"GeneName[$2];\
                                                     for(i=3;i<=NF;i++){out=out"\t"$i};print out}'\
                                $GeneAnnotationFold"GeneralInformation"\
                                $ExpressionFold"Files/"$cancer".Results" >> $ExpressionFold"Files/RNAeditingExpressionInformation"  
done          
}

function ADARAnnotation(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
DownloadRNAeditingFold=$DownloadFold"RNAediting/"
DownloadExpressionFold=$DownloadFold"GeneExpression/"
RNAeditingFile=$DownloadRNAeditingFold""$cancer"_Editing.txt"
ExpressionFile=$DownloadExpressionFold""$cancer"_expression.csv"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
ExpressionFold=$DatabaseFold"ADARCorrelation/"
ExpressionFileFold=$ExpressionFold"Files/"
ExpressionTmpFold=$ExpressionFold"Tmp/"
ExpressionCorrelationFold=$ExpressionFold"Correlation/"
ExpressionBoxplotFold=$ExpressionFold"Boxplot/"
mkdir $ExpressionFileFold
mkdir $ExpressionTmpFold
mkdir $ExpressionCorrelationFold
mkdir $ExpressionBoxplotFold
head -n 1 $RNAeditingFile | tr '\t' '\n' | awk '(NR>1){print $1"\t"NR}' > $ExpressionTmpFold""$cancer".RNAediting.samples"
head -n 1 $ExpressionFile | tr ',' '\n' | awk '(NR>1){print $1"\t"NR}' > $ExpressionTmpFold""$cancer".Expression.samples"
awk -F"\t" 'ARGIND==1{split($1,a,"-");if(substr(a[4],1,1)=="0"){ENSG[a[1]"-"a[2]"-"a[3]"-"a[4]]=ENSG[a[1]"-"a[2]"-"a[3]"-"a[4]]";"$2}}\
            ARGIND==2{if($1 in ENSG){print $0"\t"$1"\t"substr(ENSG[$1],2)}}'\
            $ExpressionTmpFold""$cancer".Expression.samples"\
            $ExpressionTmpFold""$cancer".RNAediting.samples" > $ExpressionTmpFold""$cancer".RNAediting.Expression.tumor.samples"
RNAeditingColumn=`cut -f2 $ExpressionTmpFold""$cancer".RNAediting.Expression.tumor.samples"`
ExpressionColumn=`cut -f4 $ExpressionTmpFold""$cancer".RNAediting.Expression.tumor.samples"`
awk -F"\t" -v RNAeditingColumn="${RNAeditingColumn[*]}" 'ARGIND==1{RNAediting[$1]=$1}\
            ARGIND==2 && ($1 in RNAediting){split(RNAeditingColumn,col," ");out=$1;for(i=1;i<=length(col);i++){out=out"\t"$col[i]};print out"\n"out"\n"out}'\
            $GeneAnnotationFold""$cancer".Table"\
            $RNAeditingFile > $ExpressionTmpFold""$cancer".RNAediting.tumor"
awk -F"\t" -v ExpressionColumn="${ExpressionColumn[*]}" 'ARGIND==1{split($0,a,",");if(a[1]=="ENSG00000160710" || a[1]=="ENSG00000197381" || a[1]=="ENSG00000185736"){split(ExpressionColumn,col," ");out=a[1];for(i=1;i<=length(col);i++){split(col[i],b,";");m=0;for(j=1;j<=length(b);j++){m=m+a[b[j]]};out=out"\t"m/length(b)};ENSG_value[a[1]]=out}}\
            ARGIND==2 && (FNR % 3==1){for(i in ENSG_value){print ENSG_value[i]}}'\
            $ExpressionFile\
            $ExpressionTmpFold""$cancer".RNAediting.tumor" > $ExpressionTmpFold""$cancer".ADAR.tumor"
mkdir $ExpressionBoxplotFold""$cancer"/"
mkdir $ExpressionCorrelationFold""$cancer"/"
Rscript ADARCorrelation.R $ExpressionTmpFold""$cancer".RNAediting.tumor" $ExpressionTmpFold""$cancer".ADAR.tumor" $ExpressionBoxplotFold""$cancer"/" $ExpressionCorrelationFold""$cancer"/" $ExpressionFileFold""$cancer".ADAR.Results"
}

function ADARHeatmap(){
ExpressionFold=$DatabaseFold"ADARCorrelation/"
ExpressionHeatmapFold=$ExpressionFold"Heatmap/"
ExpressionTmpFold=$ExpressionFold"Tmp/"
mkdir $ExpressionHeatmapFold
mkdir $ExpressionTmpFold
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
rm -f $ExpressionTmpFold"ADAR1.Heatmap.data"
rm -f $ExpressionTmpFold"ADAR2.Heatmap.data"
rm -f $ExpressionTmpFold"ADAR3.Heatmap.data"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
ExpressionFile=$ExpressionFold"Files/"$cancer".ADAR.Results"
awk -v cancer="$cancer" 'ARGIND==1{Editing[$2]=$1; if(!(Gene[$2]~$3)){Gene[$2]=Gene[$2]";"$3}}\
                         ARGIND==2 && (FNR>1) && ($3!="NA") && ($3<0.05) && ($7!="NA") && ($7<0.05) && ($8>0) && ($2=="ENSG00000160710"){\
                                   split(Gene[$1],a,";");for(i=2;i<=length(a);i++){print cancer"\t"Editing[$1]"\t"$1"\t"a[i]"\t"$7"\t"$8}}'\
                         $GeneAnnotationFold"GeneralInformation"\
                         $ExpressionFile >> $ExpressionTmpFold"ADAR1.Heatmap.data"
awk -v cancer="$cancer" 'ARGIND==1{Editing[$2]=$1; if(!(Gene[$2]~$3)){Gene[$2]=Gene[$2]";"$3}}\
                         ARGIND==2 && (FNR>1) && ($3!="NA") && ($3<0.05) && ($7!="NA") && ($7<0.05) && ($8>0) && ($2=="ENSG00000197381"){\
                                   split(Gene[$1],a,";");for(i=2;i<=length(a);i++){print cancer"\t"Editing[$1]"\t"$1"\t"a[i]"\t"$7"\t"$8}}'\
                         $GeneAnnotationFold"GeneralInformation"\
                         $ExpressionFile >> $ExpressionTmpFold"ADAR2.Heatmap.data"
awk -v cancer="$cancer" 'ARGIND==1{Editing[$2]=$1; if(!(Gene[$2]~$3)){Gene[$2]=Gene[$2]";"$3}}\
                         ARGIND==2 && (FNR>1) && ($3!="NA") && ($3<0.05) && ($7!="NA") && ($7<0.05) && ($8>0) && ($2=="ENSG00000185736"){\
                                   split(Gene[$1],a,";");for(i=2;i<=length(a);i++){print cancer"\t"Editing[$1]"\t"$1"\t"a[i]"\t"$7"\t"$8}}'\
                         $GeneAnnotationFold"GeneralInformation"\
                         $ExpressionFile >> $ExpressionTmpFold"ADAR3.Heatmap.data"
done
mkdir $ExpressionHeatmapFold"ADAR1/"
Rscript Heatmap.R $ExpressionTmpFold"ADAR1.Heatmap.data" 1 4 2 6 0 99 $ExpressionHeatmapFold"ADAR1/"
mkdir $ExpressionHeatmapFold"ADAR2/"
Rscript Heatmap.R $ExpressionTmpFold"ADAR2.Heatmap.data" 1 4 2 6 0 99 $ExpressionHeatmapFold"ADAR2/"
mkdir $ExpressionHeatmapFold"ADAR3/"
Rscript Heatmap.R $ExpressionTmpFold"ADAR3.Heatmap.data" 1 4 2 6 0 99 $ExpressionHeatmapFold"ADAR3/"
}

function MergeADAR(){
ExpressionFold=$DatabaseFold"ADARCorrelation/"
ExpressionTmpFold=$ExpressionFold"Tmp/"
echo -e "Cancer\tCAeditomeID\tRNAediting\tGene\tADAR1_CorP\tADAR1_CorR\tADAR2_CorP\tADAR2_CorR\tADAR3_CorP\tADAR3_CorR" > $ExpressionFold"Files/ADARInformation"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
ExpressionFile=$ExpressionFold"Files/"$cancer".ADAR.Results"
DifferentialFile=$DatabaseFold"Differential/Files/"$cancer".DifferentialResults"
awk -v cancer="$cancer" 'ARGIND==1{Editing[$2]=$1; if(!(Gene[$2]~$3)){Gene[$2]=Gene[$2]";"$3}}\
                         ARGIND==2 && (FNR>1) && ($5>=50){EditingNum[$1]=$1}\
                         ARGIND==3 && (FNR>1) && ($1 in EditingNum) && ($7!="NA") && ($7<0.05) && ($8>0){\
                                   split(Gene[$1],a,";");for(i=2;i<=length(a);i++){All[cancer"\t"Editing[$1]"\t"$1"\t"a[i]]=cancer"\t"Editing[$1]"\t"$1"\t"a[i];\
                                                          if($2=="ENSG00000160710"){ADAR1[cancer"\t"Editing[$1]"\t"$1"\t"a[i]]=$7"\t"$8}\
                                                          else if($2=="ENSG00000197381"){ADAR2[cancer"\t"Editing[$1]"\t"$1"\t"a[i]]=$7"\t"$8}\
                                                          else if($2=="ENSG00000185736"){ADAR3[cancer"\t"Editing[$1]"\t"$1"\t"a[i]]=$7"\t"$8}}}\
                                   END{for(i in All){out=i;if(i in ADAR1){out=out"\t"ADAR1[i]}else{out=out"\t.\t."};\
                                                           if(i in ADAR2){out=out"\t"ADAR2[i]}else{out=out"\t.\t."};\
                                                           if(i in ADAR3){out=out"\t"ADAR3[i]}else{out=out"\t.\t."};print out}}'\
                         $GeneAnnotationFold"GeneralInformation"\
                         $DifferentialFile\
                         $ExpressionFile >> $ExpressionFold"Files/ADARInformation"
done        
}

function DifferentialExpressionAnnotation(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
DownloadExpressionFold=$DownloadFold"GeneExpression/"
ExpressionFile=$DownloadExpressionFold""$cancer"_expression.csv"
ExpressionFold=$DatabaseFold"DifferentialExpression/"
ExpressionFileFold=$ExpressionFold"Files/"
ExpressionTmpFold=$ExpressionFold"Tmp/"
ExpressionBoxplotFold=$ExpressionFold"Boxplot/"
mkdir $ExpressionFileFold
mkdir $ExpressionTmpFold
mkdir $ExpressionBoxplotFold
head -n 1 $ExpressionFile | tr ',' '\n' | awk '(NR>1){split($1,a,"-");if(substr(a[4],1,1)==0){print "Tumor"}else if(substr(a[4],1,1)==1){print "Normal"}}' > $ExpressionTmpFold""$cancer".group"
#Files without header
mkdir $ExpressionBoxplotFold""$cancer"/"
Rscript DifferentialExpression.R $ExpressionFile $ExpressionTmpFold""$cancer".group" $ExpressionFileFold""$cancer".DifferentialResults" $ExpressionBoxplotFold""$cancer"/"
}

function StageAssociationExpression(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
DownloadExpressionFold=$DownloadFold"GeneExpression/"
ExpressionFile=$DownloadExpressionFold""$cancer"_expression.csv"
DownloadClinicalFold=$DownloadFold"Clinical/"
ClinicalFile=$DownloadClinicalFold""$cancer".txt"
StageAssociationFold=$DatabaseFold"StageExpressionAssociation/"
StageAssociationFileFold=$StageAssociationFold"Files/"
StageAssociationTmpFold=$StageAssociationFold"Tmp/"
StageAssociationCorrelationFold=$StageAssociationFold"Correlation/"
mkdir $StageAssociationFileFold
mkdir $StageAssociationTmpFold
mkdir $StageAssociationCorrelationFold
awk -F"," '{if(FNR==1){out=$1;for(i=2;i<=NF;i++){split($i,a,"-");if(substr(a[4],1,1)==0){Group[i]="Tumor";out=out"\t"$i}else if(substr(a[4],1,1)==1){Group[i]="Normal"}};print out}\
             else{out=$1;for(i=2;i<=NF;i++){if(Group[i]=="Tumor"){out=out"\t"$i}};print out}}'\
            $ExpressionFile > $StageAssociationTmpFold""$cancer".Expression" #contain tumor samples and RNAediting in Table file
head -n 1 $StageAssociationTmpFold""$cancer".Expression" | tr '\t' '\n' | awk '(NR>1){print $1"\t"NR}' > $StageAssociationTmpFold""$cancer".Expression.samples"
PathologicalStageColumn=`head -n 1 $ClinicalFile | tr '\t' '\n' | grep -n "pathologic_stage" | cut -d":" -f1`
if [[ -n $PathologicalStageColumn ]]
then
awk -F"\t" -v stageindex="$PathologicalStageColumn" 'ARGIND==1{gsub(/stage |a|b|c|1|2|3/,"",$stageindex);Clinical[$1]=$stageindex}\
                                                     ARGIND==2{split($1,a,"-");sample=a[1]"-"a[2]"-"a[3];\
                                                               if(sample in Clinical){print $0"\t"sample"\t"Clinical[sample]}else{print $0"\tNA\tNA"}}'\
                                                     $ClinicalFile\
                                                     $StageAssociationTmpFold""$cancer".Expression.samples" > $StageAssociationTmpFold""$cancer".Expression.PathStage"
Rscript Roman2Num.R $StageAssociationTmpFold""$cancer".Expression.PathStage" 4 $StageAssociationTmpFold""$cancer".Expression.PathStage.roman2num"
mkdir $StageAssociationCorrelationFold"/Pathological/"
mkdir $StageAssociationCorrelationFold"/Pathological/"$cancer"/"
Rscript StageExpressionCorrelation.R $StageAssociationTmpFold""$cancer".Expression" $StageAssociationTmpFold""$cancer".Expression.PathStage.roman2num" 5 $StageAssociationCorrelationFold"/Pathological/"$cancer"/" $StageAssociationFileFold"/"$cancer".PathStage.CorResults"
fi
ClinicalStageColumn=`head -n 1 $ClinicalFile | tr '\t' '\n' | grep -n "clinical_stage" | cut -d":" -f1`
if [[ -n $ClinicalStageColumn ]]
then
awk -F"\t" -v stageindex="$ClinicalStageColumn" 'ARGIND==1{gsub(/stage |a|b|c|1|2|3/,"",$stageindex);Clinical[$1]=$stageindex}\
                                                 ARGIND==2{split($1,a,"-");sample=a[1]"-"a[2]"-"a[3];\
                                                           if(sample in Clinical){print $0"\t"sample"\t"Clinical[sample]}else{print $0"\tNA\tNA"}}'\
                                                 $ClinicalFile\
                                                 $StageAssociationTmpFold""$cancer".Expression.samples" > $StageAssociationTmpFold""$cancer".Expression.CliStage"
Rscript Roman2Num.R $StageAssociationTmpFold""$cancer".Expression.CliStage" 4 $StageAssociationTmpFold""$cancer".Expression.CliStage.roman2num"
mkdir $StageAssociationCorrelationFold"/ClinicalStage/"
mkdir $StageAssociationCorrelationFold"/ClinicalStage/"$cancer"/"
Rscript StageExpressionCorrelation.R $StageAssociationTmpFold""$cancer".Expression" $StageAssociationTmpFold""$cancer".Expression.CliStage.roman2num" 5 $StageAssociationCorrelationFold"/ClinicalStage/"$cancer"/"  $StageAssociationFileFold"/"$cancer".CliStage.CorResults"
fi
}


function SurvivalExpressionAssociation(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
DownloadExpressionFold=$DownloadFold"GeneExpression/"
ExpressionFile=$DownloadExpressionFold""$cancer"_expression.csv"
DownloadClinicalFold=$DownloadFold"Clinical/"
ClinicalFile=$DownloadClinicalFold""$cancer".txt"
SurvivalAssociationFold=$DatabaseFold"SurvivalExpressionAssociation/"
SurvivalAssociationFileFold=$SurvivalAssociationFold"Files/"
SurvivalAssociationTmpFold=$SurvivalAssociationFold"Tmp/"
SurvivalAssociationKMFold=$SurvivalAssociationFold"KMplot/"
mkdir $SurvivalAssociationFileFold
mkdir $SurvivalAssociationTmpFold
mkdir $SurvivalAssociationKMFold
mkdir $SurvivalAssociationKMFold""$cancer"/"
awk -F"," '{if(FNR==1){out=$1;for(i=2;i<=NF;i++){split($i,a,"-");if(substr(a[4],1,1)==0){Group[i]="Tumor";out=out"\t"$i}else if(substr(a[4],1,1)==1){Group[i]="Normal"}};print out}\
             else{out=$1;for(i=2;i<=NF;i++){if(Group[i]=="Tumor"){out=out"\t"$i}};print out}}'\
            $ExpressionFile > $SurvivalAssociationTmpFold""$cancer".Expression" #contain tumor samples and RNAediting in Table file
head -n 1 $SurvivalAssociationTmpFold""$cancer".Expression" | tr '\t' '\n' | awk '(NR>1){print $1"\t"NR}' > $SurvivalAssociationTmpFold""$cancer".Expression.samples"          
#Survival Information
OS_time=`head -n 1 $ClinicalFile | tr '\t' '\n' | grep -n "days_to_death" | cut -d":" -f1`
OS_status=`head -n 1 $ClinicalFile | tr '\t' '\n' | grep -n "vital_status" | cut -d":" -f1`
awk -F"\t" -v Time="$OS_time" -v Status="$OS_status" 'ARGIND==1{Clinical[$1]=$Time"\t"$Status}\
                                                     ARGIND==2{split($1,a,"-");sample=a[1]"-"a[2]"-"a[3];\
                                                              if(sample in Clinical){print $0"\t"sample"\t"Clinical[sample]}else{print $0"\tNA\tNA\tNA"}}'\
                                                     $ClinicalFile\
                                                     $SurvivalAssociationTmpFold""$cancer".Expression.samples" > $SurvivalAssociationTmpFold""$cancer".Expression.survival"
Rscript SurvivalExpressionAnalysis.R $SurvivalAssociationTmpFold""$cancer".Expression" $SurvivalAssociationTmpFold""$cancer".Expression.survival" 4 5 $SurvivalAssociationKMFold""$cancer"/" $SurvivalAssociationFileFold"/"$cancer".survival.results"
}

function main(){
index=$1
ExpressionAnnotation $index
ADARAnnotation $index
}
