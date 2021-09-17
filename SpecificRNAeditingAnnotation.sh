cancertypes=("BLCA"  "LGG"  "BRCA"  "LIHC"  "CESC"  "LUAD"  "COAD"  "LUSC"  "GBM"  "PRAD"  "HNSC"  "STAD"  "KICH"  "THCA"  "KIRC"  "UCEC"  "KIRP"  "ACC"  "DLBC"  "ESCA"  "LAML"  "OV"  "PAAD"  "PCPG"  "SKCM"  "UCS"  "UVM"  "CHOL"  "MESO"  "READ"  "SARC"  "TGCT"  "THYM")
DatabaseFold="/data9/swu13/CAeditome_hg38/"

function MkdirPath(){
mkdir $DatabaseFold"Differential/"
mkdir $DatabaseFold"StageAssociation/"
mkdir $DatabaseFold"SurvivalAssociation/"
}


function DifferentialAnnotation(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
DownloadRNAeditingFold=$DownloadFold"RNAediting/"
RNAeditingFile=$DownloadRNAeditingFold""$cancer"_Editing.txt"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
DifferentialFold=$DatabaseFold"Differential/"
DifferentialFileFold=$DifferentialFold"Files/"
DifferentialTmpFold=$DifferentialFold"Tmp/"
DifferentialBoxplotFold=$DifferentialFold"Boxplot/"
mkdir $DifferentialFileFold
mkdir $DifferentialTmpFold
mkdir $DifferentialBoxplotFold
mkdir $DifferentialBoxplotFold""$cancer"/"
head -n 1 $RNAeditingFile | tr '\t' '\n' | awk '(NR>1){split($1,a,"-");if(substr(a[4],1,1)==0){print "Tumor"}else if(substr(a[4],1,1)==1){print "Normal"}}' > $DifferentialTmpFold""$cancer".group"
awk -F"\t" 'ARGIND==1 && (FNR>1){RNAediting[$1]=$1}ARGIND==2 && (($1 in RNAediting) || (FNR==1)){print $0}' $GeneAnnotationFold""$cancer".Table" $RNAeditingFile > $DifferentialTmpFold""$cancer".RNAediting"
#Files without header
Rscript DifferentialEditing.R $DifferentialTmpFold""$cancer".RNAediting" $DifferentialTmpFold""$cancer".group" $DifferentialFileFold""$cancer".DifferentialResults" $DifferentialBoxplotFold""$cancer"/"
}

function DifferentialHeatmap(){
DifferentialFold=$DatabaseFold"Differential/"
DifferentialHeatmapFold=$DifferentialFold"Heatmap/"
DifferentialTmpFold=$DifferentialFold"Tmp/"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
mkdir $DifferentialHeatmapFold
rm -f $DifferentialTmpFold"Heatmap.data"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
DifferentialFile=$DifferentialFold"Files/"$cancer".DifferentialResults"
awk -v cancer="$cancer" 'ARGIND==1{if(!(Gene[$2]~$3)){Gene[$2]=Gene[$2]";"$3};Editing[$2]=$1}\
                         ARGIND==2 && (FNR>1) && ($3!="NA") && ($4!="NA") && ($2<0.05){split(substr(Gene[$1],2),a,";");for(i=1;i<=length(a);i++){print cancer"\t"Editing[$1]"\t"$1"\t"a[i]"\t"$3-$4}}'\
                         $GeneAnnotationFold"GeneralInformation" \
                         $DifferentialFile >> $DifferentialTmpFold"Heatmap.data"
done
Rscript Heatmap.R $DifferentialTmpFold"Heatmap.data" 1 4 2 5 0 99 $DifferentialHeatmapFold #cancer\tGene\tEditing\tvalue
}

function StageAssociationAnnotation(){
#Install ggpubr: conda install -c conda-forge r-ggpubr
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
DownloadRNAeditingFold=$DownloadFold"RNAediting/"
DownloadClinicalFold=$DownloadFold"Clinical/"
RNAeditingFile=$DownloadRNAeditingFold""$cancer"_Editing.txt"
ClinicalFile=$DownloadClinicalFold""$cancer".txt"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
StageAssociationFold=$DatabaseFold"StageAssociation/"
StageAssociationFileFold=$StageAssociationFold"Files/"
StageAssociationTmpFold=$StageAssociationFold"Tmp/"
StageAssociationCorrelationFold=$StageAssociationFold"Correlation/"
mkdir $StageAssociationFileFold
mkdir $StageAssociationTmpFold
mkdir $StageAssociationCorrelationFold
awk -F"\t" 'ARGIND==1 && (FNR>1){RNAediting[$1]=$1}\
            ARGIND==2{if(FNR==1){out=$1;for(i=2;i<=NF;i++){split($i,a,"-");\
                                               if(substr(a[4],1,1)==0){Group[i]="Tumor";out=out"\t"$i}else if(substr(a[4],1,1)==1){Group[i]="Normal"}};print out}\
                      else{if($1 in RNAediting){out=$1;for(i=2;i<=NF;i++){if(Group[i]=="Tumor"){out=out"\t"$i}};print out}}}'\
            $GeneAnnotationFold""$cancer".Table"\
            $RNAeditingFile > $StageAssociationTmpFold""$cancer".RNAediting" #contain tumor samples and RNAediting in Table file
head -n 1 $StageAssociationTmpFold""$cancer".RNAediting" | tr '\t' '\n' | awk '(NR>1){print $1"\t"NR}' > $StageAssociationTmpFold""$cancer".editing.samples"
PathologicalStageColumn=`head -n 1 $ClinicalFile | tr '\t' '\n' | grep -n "pathologic_stage" | cut -d":" -f1`
if [[ -n $PathologicalStageColumn ]]
then
awk -F"\t" -v stageindex="$PathologicalStageColumn" 'ARGIND==1{gsub(/stage |a|b|c|1|2|3/,"",$stageindex);Clinical[$1]=$stageindex}\
                                                     ARGIND==2{split($1,a,"-");sample=a[1]"-"a[2]"-"a[3];\
                                                               if(sample in Clinical){print $0"\t"sample"\t"Clinical[sample]}else{print $0"\tNA\tNA"}}'\
                                                     $ClinicalFile\
                                                     $StageAssociationTmpFold""$cancer".editing.samples" > $StageAssociationTmpFold""$cancer".editing.PathStage"
Rscript Roman2Num.R $StageAssociationTmpFold""$cancer".editing.PathStage" 4 $StageAssociationTmpFold""$cancer".editing.PathStage.roman2num"
mkdir $StageAssociationCorrelationFold"/Pathological/"
mkdir $StageAssociationCorrelationFold"/Pathological/"$cancer"/"
Rscript StageCorrelation.R $StageAssociationTmpFold""$cancer".RNAediting" $StageAssociationTmpFold""$cancer".editing.PathStage.roman2num" 5 $StageAssociationCorrelationFold"/Pathological/"$cancer"/" $StageAssociationFileFold"/"$cancer".PathStage.CorResults"
fi
ClinicalStageColumn=`head -n 1 $ClinicalFile | tr '\t' '\n' | grep -n "clinical_stage" | cut -d":" -f1`
if [[ -n $ClinicalStageColumn ]]
then
awk -F"\t" -v stageindex="$ClinicalStageColumn" 'ARGIND==1{gsub(/stage |a|b|c|1|2|3/,"",$stageindex);Clinical[$1]=$stageindex}\
                                                 ARGIND==2{split($1,a,"-");sample=a[1]"-"a[2]"-"a[3];\
                                                           if(sample in Clinical){print $0"\t"sample"\t"Clinical[sample]}else{print $0"\tNA\tNA"}}'\
                                                 $ClinicalFile\
                                                 $StageAssociationTmpFold""$cancer".editing.samples" > $StageAssociationTmpFold""$cancer".editing.CliStage"
Rscript Roman2Num.R $StageAssociationTmpFold""$cancer".editing.CliStage" 4 $StageAssociationTmpFold""$cancer".editing.CliStage.roman2num"
mkdir $StageAssociationCorrelationFold"/ClinicalStage/"
mkdir $StageAssociationCorrelationFold"/ClinicalStage/"$cancer"/"
Rscript StageCorrelation.R $StageAssociationTmpFold""$cancer".RNAediting" $StageAssociationTmpFold""$cancer".editing.CliStage.roman2num" 5 $StageAssociationCorrelationFold"/ClinicalStage/"$cancer"/"  $StageAssociationFileFold"/"$cancer".CliStage.CorResults"
fi
}

function StageAssociationHeatmap(){
StageAssociationFold=$DatabaseFold"StageAssociation/"
StageAssociationHeatmapFold=$StageAssociationFold"Heatmap/"
StageAssociationTmpFold=$StageAssociationFold"Tmp/"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
mkdir $StageAssociationHeatmapFold
#Pathological stage
rm -f $StageAssociationTmpFold"Heatmap.Path"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
StageAssociationFile=$StageAssociationFold"Files/"$cancer".PathStage.CorResults"
if [[ -f $StageAssociationFile ]]
then
awk -v cancer="$cancer" 'ARGIND==1{if(!(Gene[$2]~$3)){Gene[$2]=Gene[$2]";"$3};Editing[$2]=$1}\
                         ARGIND==2 && (FNR>1) && ($2!="NA") && ($2<0.05){split(substr(Gene[$1],2),a,";");for(i=1;i<=length(a);i++){print cancer"\t"Editing[$1]"\t"$1"\t"a[i]"\t"$3}}'\
                         $GeneAnnotationFold"GeneralInformation"\
                         $StageAssociationFile >> $StageAssociationTmpFold"Heatmap.Path"
fi
done
mkdir $StageAssociationHeatmapFold"Path/"
Rscript Heatmap.R $StageAssociationTmpFold"Heatmap.Path" 1 4 2 5 0 99 $StageAssociationHeatmapFold"Path/" #cancer\tGene\tEditing\tvalue
#Clinical stage
rm -f $StageAssociationTmpFold"Heatmap.Cli"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
StageAssociationFile=$StageAssociationFold"Files/"$cancer".CliStage.CorResults"
if [[ -f $StageAssociationFile ]]
then
awk -v cancer="$cancer" 'ARGIND==1{if(!(Gene[$2]~$3)){Gene[$2]=Gene[$2]";"$3};Editing[$2]=$1}\
                         ARGIND==2 && (FNR>1) && ($2!="NA") && ($2<0.05){split(substr(Gene[$1],2),a,";");for(i=1;i<=length(a);i++){print cancer"\t"Editing[$1]"\t"$1"\t"a[i]"\t"$3}}'\
                         $GeneAnnotationFold"GeneralInformation" \
                         $StageAssociationFile >> $StageAssociationTmpFold"Heatmap.Cli"
fi
done
mkdir $StageAssociationHeatmapFold"Cli/"
Rscript Heatmap.R $StageAssociationTmpFold"Heatmap.Cli" 1 4 2 5 0 99 $StageAssociationHeatmapFold"Cli/" #cancer\tGene\tEditing\tvalue
}

function SurvivalAssociationAnnotation(){
#Install survival, survminer: install.packages("survival"),install.packages("survminer") or  conda install -c conda-forge r-survminer
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
DownloadRNAeditingFold=$DownloadFold"RNAediting/"
DownloadClinicalFold=$DownloadFold"Clinical/"
RNAeditingFile=$DownloadRNAeditingFold""$cancer"_Editing.txt"
ClinicalFile=$DownloadClinicalFold""$cancer".txt"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
SurvivalAssociationFold=$DatabaseFold"SurvivalAssociation/"
SurvivalAssociationFileFold=$SurvivalAssociationFold"Files/"
SurvivalAssociationTmpFold=$SurvivalAssociationFold"Tmp/"
SurvivalAssociationKMFold=$SurvivalAssociationFold"KMplot/"
mkdir $SurvivalAssociationFileFold
mkdir $SurvivalAssociationTmpFold
mkdir $SurvivalAssociationKMFold
mkdir $SurvivalAssociationKMFold""$cancer"/"
awk -F"\t" 'ARGIND==1 && (FNR>1){RNAediting[$1]=$1}\
            ARGIND==2{if(FNR==1){out=$1;for(i=2;i<=NF;i++){split($i,a,"-");\
                                                          if(substr(a[4],1,1)==0){Group[i]="Tumor";out=out"\t"$i}else if(substr(a[4],1,1)==1){Group[i]="Normal"}};print out}\
                      else{if($1 in RNAediting){out=$1;for(i=2;i<=NF;i++){if(Group[i]=="Tumor"){out=out"\t"$i}};print out}}}'\
            $GeneAnnotationFold""$cancer".Table"\
            $RNAeditingFile > $SurvivalAssociationTmpFold""$cancer".RNAediting" #contain tumor samples and RNAediting in Table file
head -n 1 $SurvivalAssociationTmpFold""$cancer".RNAediting" | tr '\t' '\n' | awk '(NR>1){print $1"\t"NR}' > $SurvivalAssociationTmpFold""$cancer".editing.samples"          
#Survival Information
OS_time=`head -n 1 $ClinicalFile | tr '\t' '\n' | grep -n "days_to_death" | cut -d":" -f1`
OS_status=`head -n 1 $ClinicalFile | tr '\t' '\n' | grep -n "vital_status" | cut -d":" -f1`
awk -F"\t" -v Time="$OS_time" -v Status="$OS_status" 'ARGIND==1{Clinical[$1]=$Time"\t"$Status}\
                                                     ARGIND==2{split($1,a,"-");sample=a[1]"-"a[2]"-"a[3];\
                                                              if(sample in Clinical){print $0"\t"sample"\t"Clinical[sample]}else{print $0"\tNA\tNA\tNA"}}'\
                                                     $ClinicalFile\
                                                     $SurvivalAssociationTmpFold""$cancer".editing.samples" > $SurvivalAssociationTmpFold""$cancer".editing.survival"
Rscript SurvivalAnalysis.R $SurvivalAssociationTmpFold""$cancer".RNAediting" $SurvivalAssociationTmpFold""$cancer".editing.survival" 4 5 $SurvivalAssociationKMFold""$cancer"/" $SurvivalAssociationFileFold"/"$cancer".survival.results"
}

function SurvivalAssociationHeatmap(){
SurvivalAssociationFold=$DatabaseFold"SurvivalAssociation/"
SurvivalAssociationHeatmapFold=$SurvivalAssociationFold"Heatmap/"
SurvivalAssociationTmpFold=$SurvivalAssociationFold"Tmp/"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
mkdir $SurvivalAssociationHeatmapFold
#Survival
rm -f $SurvivalAssociationTmpFold"Heatmap.Data"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
SurvivalAssociationFile=$SurvivalAssociationFold"Files/"$cancer".survival.results"
awk -v cancer="$cancer" 'ARGIND==1{if(!(Gene[$2]~$3)){Gene[$2]=Gene[$2]";"$3};Editing[$2]=$1}\
                         ARGIND==2 && (FNR>1) && ($2!="NA") && ($2<0.05) && ($7!="NA") && ($7<0.05){split(substr(Gene[$1],2),a,";");for(i=1;i<=length(a);i++){print cancer"\t"Editing[$1]"\t"$1"\t"a[i]"\t"$8}}'\
                         $GeneAnnotationFold"GeneralInformation" \
                         $SurvivalAssociationFile >> $SurvivalAssociationTmpFold"Heatmap.Data"
done
Rscript Heatmap.R $SurvivalAssociationTmpFold"Heatmap.Data" 1 4 2 5 1 30 $SurvivalAssociationHeatmapFold #cancer\tGene\tEditing\tvalue\tDifferenatialNumber
}

function GeneBrowserAnnotation(){
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
GeneAnnotationTmpFold=$GeneAnnotationFold"Tmp/"
GeneBrowserFold=$GeneAnnotationFold"Bed/"
mkdir $GeneBrowserFold
GeneralFile=$GeneAnnotationFold"GeneralInformation"
DownloadFold=$DatabaseFold"RawData/"
ReferenceFold=$DownloadFold"Reference/"
ReferenceGTFFile=$ReferenceFold"gencode.v22.annotation.gtf"
rm -f $GeneAnnotationFold"Genome.browser.data"
DifferentialFold=$DatabaseFold"Differential/"
StageAssociationFold=$DatabaseFold"StageAssociation/"
SurvivalAssociationFold=$DatabaseFold"SurvivalAssociation/"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
DifferentialFile=$DifferentialFold"Files/"$cancer".DifferentialResults"
awk -F"\t" -v cancer="$cancer" '(($2=="NA" && $5>=5 && $6==0&& $8>=5) || ($2=="NA" && $5==0 && $6>=5 && $7>=5)  || ($2<0.05 && $5+$6>=50)) && (FNR>1){print cancer"\t"$1}' $DifferentialFile >> $GeneAnnotationFold"Genome.browser.data"
StageAssociationFile=$StageAssociationFold"Files/"$cancer".PathStage.CorResults"
if [[ -f $StageAssociationFile ]]
then
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && ($5>=50) && (FNR>1){Specific[$1]=$1}\
                                ARGIND==2 && ($2<0.05) && ($1 in Specific) && (FNR>1){print cancer"\t"$1}'\
                                $DifferentialFile\
                                $StageAssociationFile >> $GeneAnnotationFold"Genome.browser.data"
fi
StageAssociationFile=$StageAssociationFold"Files/"$cancer".CliStage.CorResults"
if [[ -f $StageAssociationFile ]]
then
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && ($5>=50) && (FNR>1){Specific[$1]=$1}\
                                ARGIND==2 && ($2<0.05) && ($1 in Specific) && (FNR>1){print cancer"\t"$1}'\
                                $DifferentialFile\
                                $StageAssociationFile >> $GeneAnnotationFold"Genome.browser.data"
fi
SurvivalAssociationFile=$SurvivalAssociationFold"Files/"$cancer".survival.results"
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && ($5>=50) && (FNR>1){Specific[$1]=$1}\
                                ARGIND==2 && ($2<0.05) && ($7<0.05) && ($1 in Specific) && (FNR>1){print cancer"\t"$1}'\
                                $DifferentialFile\
                                $SurvivalAssociationFile >> $GeneAnnotationFold"Genome.browser.data"
done
awk -F"\t" 'ARGIND==1{if(!(Editing[$2]~$1)){Editing[$2]=Editing[$2]";"$1}}\
            ARGIND==2 && ($3=="gene"){gsub(/"| /,"",$9);split($9,a,";");GeneID=substr(a[1],length("gene_id")+1);Gene[GeneID]=$4"-"$5}\
            ARGIND==3 && ($2 in Editing){split(substr(Editing[$2],2),a,";");for(i=1;i<=length(a);i++){print $1"\t"$2"\t"$3"\t"Gene[$3]"\t"a[i]}}'\
            $GeneAnnotationFold"Genome.browser.data"\
            $ReferenceGTFFile\
            $GeneAnnotationFold"GeneralInformation" | sort | uniq > $GeneAnnotationFold"Genome.browser.data.1"
mv  $GeneAnnotationFold"Genome.browser.data.1" $GeneAnnotationFold"Genome.browser.data"
bash generate_ucsc_bed.sh $GeneAnnotationFold"Genome.browser.data" $GeneBrowserFold $GeneAnnotationTmpFold
}

function MergeAll(){
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
DifferentialFold=$DatabaseFold"Differential/"
StageAssociationFold=$DatabaseFold"StageAssociation/"
SurvivalAssociationFold=$DatabaseFold"SurvivalAssociation/"
head -n 1 $DifferentialFold"Files/BLCA"* | awk -F"\t" '{print "Cancer\tCAeditomeID\t"$0}' > $DifferentialFold"Files/DifferentialResults"
head -n 1 $StageAssociationFold"Files/BLCA"* | awk -F"\t" '{print "Cancer\tStageType\tCAeditomeID\t"$0}' > $StageAssociationFold"Files/Stage.CorResults"
head -n 1 $SurvivalAssociationFold"Files/BLCA"* | awk -F"\t" '{print "Cancer\tCAeditomeID\t"$0}' > $SurvivalAssociationFold"Files/Survival.results"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
DifferentialFile=$DifferentialFold"Files/"$cancer".DifferentialResults"
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && (FNR>1){CAeditomeID[$2]=$1}\
                                ARGIND==2 && (FNR>1){print cancer"\t"CAeditomeID[$1]"\t"$0}'\
                                $GeneAnnotationFold"GeneralInformation" $DifferentialFile >> $DifferentialFold"Files/DifferentialResults"
StageAssociationFile=$StageAssociationFold"Files/"$cancer".PathStage.CorResults"
if [[ -f $StageAssociationFile ]]
then
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && (FNR>1){CAeditomeID[$2]=$1}\
                                ARGIND==2 && (FNR>1){print cancer"\tpathologic_stage\t"CAeditomeID[$1]"\t"$0}'\
                                $GeneAnnotationFold"GeneralInformation" $StageAssociationFile >> $StageAssociationFold"Files/Stage.CorResults"
fi
StageAssociationFile=$StageAssociationFold"Files/"$cancer".CliStage.CorResults"
if [[ -f $StageAssociationFile ]]
then
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && (FNR>1){CAeditomeID[$2]=$1}\
                                ARGIND==2 && (FNR>1){print cancer"\tclinical_stage\t"CAeditomeID[$1]"\t"$0}'\
                                $GeneAnnotationFold"GeneralInformation" $StageAssociationFile >> $StageAssociationFold"Files/Stage.CorResults"
fi
SurvivalAssociationFile=$SurvivalAssociationFold"Files/"$cancer".survival.results"
awk -F"\t" -v cancer="$cancer" 'ARGIND==1 && (FNR>1){CAeditomeID[$2]=$1}\
                                ARGIND==2 && (FNR>1){print cancer"\t"CAeditomeID[$1]"\t"$0}'\
                                $GeneAnnotationFold"GeneralInformation" $SurvivalAssociationFile >> $SurvivalAssociationFold"Files/Survival.results"
done
}

function main(){
index=$1
DifferentialAnnotation $index
StageAssociationAnnotation $index
SurvivalAssociationAnnotation $index
}