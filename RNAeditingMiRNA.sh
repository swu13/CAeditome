cancertypes=("BLCA"  "LGG"  "BRCA"  "LIHC"  "CESC"  "LUAD"  "COAD"  "LUSC"  "GBM"  "PRAD"  "HNSC"  "STAD"  "KICH"  "THCA"  "KIRC"  "UCEC"  "KIRP"  "ACC"  "DLBC"  "ESCA"  "LAML"  "OV"  "PAAD"  "PCPG"  "SKCM"  "UCS"  "UVM"  "CHOL"  "MESO"  "READ"  "SARC"  "TGCT"  "THYM")
DatabaseFold="/data9/swu13/CAeditome_hg38/"

function MkdirPath(){
mkdir $DatabaseFold"miRNA/"
}

function TargetscanBindingPrediction(){
#Install targetscan
#ToolPath=$DatabaseFold"Tool/"
#wget -P $ToolPath http://www.targetscan.org/vert_72/vert_72_data_download/targetscan_70.zip
#mkdir $ToolPath"targetscan/"
#unzip $ToolPath"targetscan_70.zip" -d $ToolPath"targetscan/"
InputFile=$1
miRNAFile=$2
OutputFile=$3
ToolPath=$DatabaseFold"Tool/"
TargetscanPath=$ToolPath"targetscan/"
$TargetscanPath"targetscan_70.pl" $miRNAFile $InputFile $OutputFile
}

function miRandaBindingPrediction(){
#Install miRanda
#The miRanda database is not working. I already download the tool in the synapse and put it in the ToolPath
#ToolPath=$DatabaseFold"Tool/"
#cd $ToolPath
#synapse get syn23546820 
#tar -zxvf $ToolPath"miRanda-aug2010.tar.gz" -C $ToolPath
#cd $ToolPath"miRanda-3.3a/"
#./configure --prefix=$ToolPath"miRanda-3.3a/"
#make
#make install
InputFile=$1
miRNAFile=$2
OutputFile=$3
ToolPath=$DatabaseFold"Tool/"
miRandaPath=$ToolPath"miRanda-3.3a/bin/"
$miRandaPath"miranda" $miRNAFile $InputFile > $OutputFile
}

function BasemiRNABinding(){
DownloadFold=$DatabaseFold"RawData/"
DownloadReference=$DownloadFold"Reference/"
miRNAFold=$DatabaseFold"miRNA/"
miRNAFileFold=$miRNAFold"Files/"
miRNATmpFold=$miRNAFold"Tmp/"
mkdir $miRNAFileFold
mkdir $miRNATmpFold
awk -F"\t" '(NR % 2==1){out=substr($0,2);getline;gsub(/T/,"U",$0);print out"\t9606\t"$0}' $DownloadReference"gencode.v22.UTR3.Really.fa" > $miRNATmpFold"UTR3.Fortargetscan"
TargetscanBindingPrediction $miRNATmpFold"UTR3.Fortargetscan" $DownloadReference"hsa_miRNA_miRBase.txt" $miRNATmpFold"UTR3.targetscan.prediction"  &
awk -F"\t" '(NR % 2==1){out=$0;getline;gsub(/T/,"U",$0);print out"\n"$0}' $DownloadReference"gencode.v22.UTR3.Really.fa" > $miRNATmpFold"UTR3.FormiRanda"
miRandaBindingPrediction $miRNATmpFold"UTR3.FormiRanda" $DownloadReference"hsa_miRNA_miRBase.fa" $miRNATmpFold"UTR3.miRanda.prediction" &
awk -F"\t" '(NR % 2==1){out=substr($0,2);getline;gsub(/T/,"U",$0);print out"\t9606\t"$0}' $DownloadReference"gencode.v22.lncRNA.Really.fa" > $miRNATmpFold"lncRNA.Fortargetscan"
TargetscanBindingPrediction $miRNATmpFold"lncRNA.Fortargetscan" $DownloadReference"hsa_miRNA_miRBase.txt" $miRNATmpFold"lncRNA.targetscan.prediction"  &
awk -F"\t" '(NR % 2==1){out=$0;getline;gsub(/T/,"U",$0);print out"\n"$0}' $DownloadReference"gencode.v22.lncRNA.Really.fa" > $miRNATmpFold"lncRNA.FormiRanda"
miRandaBindingPrediction $miRNATmpFold"lncRNA.FormiRanda" $DownloadReference"hsa_miRNA_miRBase.fa" $miRNATmpFold"lncRNA.miRanda.prediction" &
}

function LncRNAMerge(){
InputName=$1
TargetscanFile=$InputName".targetscan.prediction" #$miRNATmpFold""$cancer".withoutEditing.targetscan"
miRandaFile=$InputName".miRanda.prediction" #$miRNATmpFold""$cancer".withoutEditing.miRanda"
TmpFold=$2
OutName=$3
mkdir $TmpFold
echo -e "ENST\tmiRNA\tStart\tTerminal\tScore" > $TmpFold"targetscan.prediction"
awk -F"\t" '(NR>1){split($1,a,"|");PM=a[length(a)-1];Exs=a[length(a)-2];split(Exs,exons,";");num1=0;num2=0;\
                                  if(PM=="+"){for(j=1;j<=length(exons);j++){split(exons[j],EachExon,"-");\
                                                  if(num1+EachExon[2]-EachExon[1]+1>=$4 && num1<$4){start=EachExon[1]+$4-num1-1};num1=num1+EachExon[2]-EachExon[1]+1;\
                                                  if(num2+EachExon[2]-EachExon[1]+1>=$5 && num2<$5){terminal=EachExon[1]+$5-num2-1};num2=num2+EachExon[2]-EachExon[1]+1;\
                                                  if(num1>$4 && num2>$5){break}}}\
                                  else{for(j=1;j<=length(exons);j++){split(exons[j],EachExon,"-");\
                                                  if(num1+EachExon[2]-EachExon[1]+1>=$4 && num1<$4){terminal=EachExon[2]-($4-num1)+1};num1=num1+EachExon[2]-EachExon[1]+1;\
                                                  if(num2+EachExon[2]-EachExon[1]+1>=$5 && num2<$5){start=EachExon[2]-($5-num2)+1};num2=num2+EachExon[2]-EachExon[1]+1;\
                                                  if(num1>$4 && num2>$5){break}}};\
                                  print $1"\t"$2"\t"start"\t"terminal"\t"$11}'\
                                  $TargetscanFile >> $TmpFold"targetscan.prediction"
echo -e "ENST\tmiRNA\tStart\tTerminal\tScore\tEnergyScore" > $TmpFold"miRanda.prediction"
awk -F"\t" '(substr($1,1,1)==">" && substr($1,1,2)!=">>" ){split($1,a,":");split($2,b,"|");PM=b[length(b)-1];Exs=b[length(b)-2];split(Exs,exons,";");num1=0;num2=0;split($6,m," ");\
                                  if(PM=="+"){for(j=1;j<=length(exons);j++){split(exons[j],EachExon,"-");\
                                                  if(num1+EachExon[2]-EachExon[1]+1>=m[1] && num1<m[1]){start=EachExon[1]+m[1]-num1-1};num1=num1+EachExon[2]-EachExon[1]+1;\
                                                  if(num2+EachExon[2]-EachExon[1]+1>=m[2] && num2<m[2]){terminal=EachExon[1]+m[2]-num2-1};num2=num2+EachExon[2]-EachExon[1]+1;\
                                                  if(num1>m[1] && num2>m[2]){break}}}\
                                  else{for(j=1;j<=length(exons);j++){split(exons[j],EachExon,"-");\
                                                  if(num1+EachExon[2]-EachExon[1]+1>=m[1] && num1<m[1]){terminal=EachExon[2]-(m[1]-num1)+1};num1=num1+EachExon[2]-EachExon[1]+1;\
                                                  if(num2+EachExon[2]-EachExon[1]+1>=m[2] && num2<m[2]){start=EachExon[2]-(m[2]-num2)+1};num2=num2+EachExon[2]-EachExon[1]+1;\
                                                  if(num1>m[1] && num2>m[2]){break}}};\
                                  print $2"\t"substr(a[1],2)"\t"start"\t"terminal"\t"$3"\t"$4}'\
                                  $miRandaFile >> $TmpFold"miRanda.prediction"
cp $TmpFold"targetscan.prediction" $OutName".targetscan.prediction"
cp $TmpFold"miRanda.prediction" $OutName".miRanda.prediction"
rm -r $TmpFold
} 

function UTR3Merge(){
InputName=$1
TargetscanFile=$InputName".targetscan.prediction" #$miRNATmpFold""$cancer".withoutEditing.targetscan"
miRandaFile=$InputName".miRanda.prediction" #$miRNATmpFold""$cancer".withoutEditing.miRanda"
TmpFold=$2
OutName=$3
mkdir $TmpFold
echo -e "ENST\tmiRNA\tStart\tTerminal\tScore" > $TmpFold"targetscan.prediction"
awk -F"\t" '(NR>1){split($1,a,":");split(a[length(a)],b,"(");split(b[1],c,"-");PM=substr(b[2],1,1);\
                                  if(PM=="+"){start=c[1]+$4-1;terminal=c[1]+$5-1}\
                                  else{start=c[2]-$5+1;terminal=c[2]-$4+1};\
                                  print $1"\t"$2"\t"start"\t"terminal"\t"$11}'\
                                  $TargetscanFile >> $TmpFold"targetscan.prediction"
echo -e "ENST\tmiRNA\tStart\tTerminal\tScore\tEnergyScore" > $TmpFold"miRanda.prediction"
awk -F"\t" '(substr($1,1,1)==">" && substr($1,1,2)!=">>" ){split($1,a,":");split($2,n,":");split(n[length(n)],b,"(");split(b[1],c,"-");PM=substr(b[2],1,1);split($6,m," ");\
                                  if(PM=="+"){start=c[1]+m[1]-1;terminal=c[1]+m[2]-1}\
                                  else{start=c[2]-m[2]+1;terminal=c[2]-m[1]+1};\
                                  print $2"\t"substr(a[1],2)"\t"start"\t"terminal"\t"$3"\t"$4}'\
                                  $miRandaFile >> $TmpFold"miRanda.prediction"
cp $TmpFold"targetscan.prediction" $OutName".targetscan.prediction"
cp $TmpFold"miRanda.prediction" $OutName".miRanda.prediction"
rm -r $TmpFold
}

function miRNABinding(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
DownloadFold=$DatabaseFold"RawData/"
DownloadReference=$DownloadFold"Reference/"
miRNAFold=$DatabaseFold"miRNA/"
miRNAFileFold=$miRNAFold"Files/"
miRNATmpFold=$miRNAFold"Tmp/"
#Withone RNAediting
KnownEditingFile=$miRNATmpFold"KnownEditing.Table"
RunningEditingFile=$miRNATmpFold"RunningEditing.Table"
if [[ ! -f $KnownEditingFile ]]
then
cut -f1,5 $GeneAnnotationFold""$cancer".Table" | tail -n +2 | sort | uniq > $RunningEditingFile
cp $RunningEditingFile $KnownEditingFile 
else
awk -F"\t" 'ARGIND==1{Known[$0]=$0}\
            ARGIND==2 && (!($1"\t"$5 in Known)) && (FNR>1){print $1"\t"$5}'\
            $KnownEditingENSTFile\
            $GeneAnnotationFold""$cancer".Table" | sort | uniq > $RunningEditingFile
cat $RunningEditingFile >> $KnownEditingFile 
fi
#mRNA 3UTR
awk -F"\t" '{split($1,a,"_");chromo=a[1];for(i=2;i<=length(a)-2;i++){chromo=chromo"_"a[i]};print chromo"\t"a[length(a)-1]-1"\t"a[length(a)-1]"\t.\t.\t"a[length(a)]}' $RunningEditingFile | sort | uniq > $RunningEditingFile".bed"
bedtools intersect -a $RunningEditingFile".bed" -b $DownloadReference"gencode.v22.UTR3.Really.bed" -s -wa -wb > $miRNATmpFold"RNAediting.UTR3"
awk -F"\t" 'ARGIND==1 && (FNR % 2==1){split($1,a,":");getline;Seq[substr(a[1],2)]=$0}\
            ARGIND==2{WTseq=Seq[$10];if($6=="+"){po=$3-$8}else if($6=="-"){po=$9-$2};\
                                     if(substr(WTseq,po,1)!="N"){print ">"$1"_"$3"_"$6"|"$10"\n"substr(WTseq,1,po-1)"G"substr(WTseq,po+1)}}'\
            $DownloadReference"gencode.v22.UTR3.Really.fa"\
            $miRNATmpFold"RNAediting.UTR3" > $miRNATmpFold"RNAediting.UTR3.fa"
awk -F"\t" '(NR % 2==1){out=substr($0,2);getline;gsub(/T/,"U",$0);print out"\t9606\t"$0}' $miRNATmpFold"RNAediting.UTR3.fa" > $miRNATmpFold"RNAediting.UTR3.Fortargetscan"
TargetscanBindingPrediction $miRNATmpFold"RNAediting.UTR3.Fortargetscan" $DownloadReference"hsa_miRNA_miRBase.txt" $miRNATmpFold"RNAediting.UTR3.targetscan.prediction"  &
awk -F"\t" '(NR % 2==1){out=$0;getline;gsub(/T/,"U",$0);print out"\n"$0}' $miRNATmpFold"RNAediting.UTR3.fa" > $miRNATmpFold"RNAediting.UTR3.FormiRanda"
miRandaBindingPrediction $miRNATmpFold"RNAediting.UTR3.FormiRanda" $DownloadReference"hsa_miRNA_miRBase.fa" $miRNATmpFold"RNAediting.UT3.miRanda.prediction" &
#lncRNA
awk -F"\t" 'ARGIND==1 && (FNR %2 ==1){split($0,a,"|");ENST=substr(a[1],2);strand[ENST]=a[length(a)-1];Exons[ENST]=a[length(a)-2];Chr[ENST]=a[length(a)];getline;Seq[ENST]=$0}\
            ARGIND==2 && ($2 in strand){PM=strand[$2];Exs=Exons[$2];WTseq=Seq[$2];split(Exs,exons,";");num=0;flag="N";split($1,a,"_");Editing=a[length(a)-1];\
                                        if(PM=="+"){for(j=1;j<=length(exons);j++){split(exons[j],EachExon,"-");\
                                                        if(EachExon[2]<Editing){num=num+EachExon[2]-EachExon[1]+1}\
                                                        else if(EachExon[1]<=Editing && EachExon[2]>=Editing){num=num+Editing-EachExon[1]+1;flag="Y";break}\
                                                        else{flag="N";break}}}\
                                        else{for(j=1;j<=length(exons);j++){split(exons[j],EachExon,"-");\
                                                        if(Editing<EachExon[1]){num=num+EachExon[2]-EachExon[1]+1}\
                                                        else if (EachExon[1]<=Editing && EachExon[2]>=Editing){num=num+EachExon[2]-Editing+1;flag="Y";break}\
                                                        else{flag="N";break}}}\
                                        if(flag=="Y"){print ">"$1"|"$2"|"Exs"|"PM"|"Chr[$2]"\n"substr(WTseq,1,num-1)"G"substr(WTseq,num+1)}}'\
            $DownloadReference"gencode.v22.lncRNA.Really.fa"\
            $RunningEditingFile > $miRNATmpFold"RNAediting.lncRNA.fa"
awk -F"\t" '(NR % 2==1){out=substr($0,2);getline;gsub(/T/,"U",$0);print out"\t9606\t"$0}' $miRNATmpFold"RNAediting.lncRNA.fa" > $miRNATmpFold"RNAediting.lncRNA.Fortargetscan"
TargetscanBindingPrediction $miRNATmpFold"RNAediting.lncRNA.Fortargetscan" $DownloadReference"hsa_miRNA_miRBase.txt" $miRNATmpFold"RNAediting.lncRNA.targetscan.prediction"  &
awk -F"\t" '(NR % 2==1){out=$0;getline;gsub(/T/,"U",$0);print out"\n"$0}' $miRNATmpFold"RNAediting.lncRNA.fa" > $miRNATmpFold"RNAediting.lncRNA.FormiRanda"
miRandaBindingPrediction $miRNATmpFold"RNAediting.lncRNA.FormiRanda" $DownloadReference"hsa_miRNA_miRBase.fa" $miRNATmpFold"RNAediting.lncRNA.miRanda.prediction" &
#miRNA
bedtools intersect -a $RunningEditingFile".bed" -b $DownloadReference"hsa.miRNA.bed" -s -wa -wb > $miRNATmpFold"RNAediting.miRNA.pre"
awk -F"\t" 'ARGIND==1{Seed[$1]=$2}\
            ARGIND==2{if($6=="+"){po=$3-$8}else if($6=="-"){po=$9-$2};if(po>=2 && po<=8){print $1"_"$3"_"$6"|"$10"\t"substr(Seed[$10],1,po-2)"G"substr(Seed[$10],po)"\t9606"}}'\
            $DownloadReference"hsa_miRNA_miRBase.txt"\
            $miRNATmpFold"RNAediting.miRNA.pre" > $miRNATmpFold"RNAediting.miRNA.txt"
awk -F"\t" 'ARGIND==1 && (FNR % 2 ==1){split($0,a,":");getline;Seed[substr(a[1],2)]=$0}\
            ARGIND==2{if($6=="+"){po=$3-$8}else if($6=="-"){po=$9-$2};if(po>=2 && po<=8){print ">"$1"_"$3"_"$6"|"$10"\n"substr(Seed[$10],1,po-1)"G"substr(Seed[$10],po+1)}}'\
            $DownloadReference"hsa_miRNA_miRBase.fa"\
            $miRNATmpFold"RNAediting.miRNA.pre" > $miRNATmpFold"RNAediting.miRNA.fa"
TargetscanBindingPrediction $miRNATmpFold"UTR3.Fortargetscan" $miRNATmpFold"RNAediting.miRNA.txt" $miRNATmpFold"UTR3.EditedmiRNA.targetscan.prediction"  &
miRandaBindingPrediction $miRNATmpFold"UTR3.FormiRanda" $miRNATmpFold"RNAediting.miRNA.fa" $miRNATmpFold"UTR3.EditedmiRNA.miRanda.prediction" &
TargetscanBindingPrediction $miRNATmpFold"lncRNA.Fortargetscan" $miRNATmpFold"RNAediting.miRNA.txt" $miRNATmpFold"lncRNA.EditedmiRNA.targetscan.prediction"  &
miRandaBindingPrediction $miRNATmpFold"lncRNA.FormiRanda" $miRNATmpFold"RNAediting.miRNA.fa" $miRNATmpFold"lncRNA.EditedmiRNA.miRanda.prediction" &
#check if it is complete
while true
do
JobNums_UTR_targetscan=`ps -aux | grep targetscan_70 | grep "RNAediting.UTR3.Fortargetscan" | wc -l`
JobNums_lnc_targetscan=`ps -aux | grep targetscan_70 | grep "RNAediting.lncRNA.Fortargetscan" | wc -l`
JobNums_miR_UTR_targetscan=`ps -aux | grep targetscan_70 | grep "UTR3.EditedmiRNA" | wc -l`
JobNums_miR_lnc_targetscan=`ps -aux | grep targetscan_70 | grep "lncRNA.EditedmiRNA" | wc -l`
JobNums_UTR_miRanda=`ps -aux | grep miranda | grep "RNAediting.UTR3.FormiRanda" | wc -l`
JobNums_lnc_miRanda=`ps -aux | grep miranda | grep "RNAediting.lncRNA.FormiRanda" | wc -l`
JobNums_miR_UTR_miRanda=`ps -aux | grep miranda | grep "UTR3.EditedmiRNA" | wc -l`
JobNums_miR_lnc_miRanda=`ps -aux | grep miranda | grep "lncRNA.EditedmiRNA" | wc -l`
if [[ $JobNums_UTR_targetscan -gt 0 ]] || [[ $JobNums_lnc_targetscan -gt 0 ]] || [[ $JobNums_miR_UTR_targetscan -gt 0 ]] || [[ $JobNums_miR_lnc_targetscan -gt 0 ]] || [[ $JobNums_UTR_miRanda -gt 0 ]] || [[ $JobNums_lnc_miRanda -gt 0 ]] || [[ $JobNums_miR_UTR_miRanda -gt 0 ]] || [[ $JobNums_miR_lnc_miRanda -gt 0 ]]
then
continue
else
cat $miRNATmpFold"RNAediting.UTR3.targetscan.prediction" >> $miRNAFileFold"RNAediting.UTR3.targetscan.prediction" &
cat $miRNATmpFold"RNAediting.lncRNA.targetscan.prediction" >> $miRNAFileFold"RNAediting.lncRNA.targetscan.prediction" &
cat $miRNATmpFold"RNAediting.UTR3.miRanda.prediction" >> $miRNAFileFold"RNAediting.UTR3.miRanda.prediction" &
cat $miRNATmpFold"RNAediting.lncRNA.miRanda.prediction" >> $miRNAFileFold"RNAediting.lncRNA.miRanda.prediction" &
cat $miRNATmpFold"UTR3.EditedmiRNA.targetscan.prediction" >> $miRNAFileFold"UTR3.EditedmiRNA.targetscan.prediction" &
cat $miRNATmpFold"lncRNA.EditedmiRNA.targetscan.prediction" >> $miRNAFileFold"lncRNA.EditedmiRNA.targetscan.prediction" &
cat $miRNATmpFold"UTR3.EditedmiRNA.miRanda.prediction" >> $miRNAFileFold"UTR3.EditedmiRNA.miRanda.prediction" &
cat $miRNATmpFold"lncRNA.EditedmiRNA.miRanda.prediction" >> $miRNAFileFold"lncRNA.EditedmiRNA.miRanda.prediction" &

LncRNAMerge $miRNATmpFold"RNAediting.lncRNA" $miRNATmpFold"RNAediting_lncRNA/" $miRNATmpFold"ProcessRNAediting.lncRNA"
cat $miRNATmpFold"ProcessRNAediting.lncRNA.targetscan.prediction" >> $miRNAFileFold"ProcessRNAediting.lncRNA.targetscan.prediction"
cat $miRNATmpFold"ProcessRNAediting.lncRNA.miRanda.prediction" >> $miRNAFileFold"ProcessRNAediting.lncRNA.miRanda.prediction"
UTR3Merge $miRNATmpFold"RNAediting.UTR3" $miRNATmpFold"RNAediting_UTR3/" $miRNATmpFold"ProcessRNAediting.UTR3"
cat $miRNATmpFold"ProcessRNAediting.UTR3.targetscan.prediction" >> $miRNAFileFold"ProcessRNAediting.UTR3.targetscan.prediction"
cat $miRNATmpFold"ProcessRNAediting.UTR3.miRanda.prediction" >> $miRNAFileFold"ProcessRNAediting.UTR3.miRanda.prediction"
LncRNAMerge $miRNATmpFold"lncRNA.EditedmiRNA" $miRNATmpFold"lncRNA_EditedmiRNA/" $miRNATmpFold"ProcesslncRNA.EditedmiRNA"
cat $miRNATmpFold"ProcesslncRNA.EditedmiRNA.targetscan.prediction" >> $miRNAFileFold"ProcesslncRNA.EditedmiRNA.targetscan.prediction"
cat $miRNATmpFold"ProcesslncRNA.EditedmiRNA.miRanda.prediction" >> $miRNAFileFold"ProcesslncRNA.EditedmiRNA.miRanda.prediction"
UTR3Merge $miRNATmpFold"UTR3.EditedmiRNA" $miRNATmpFold"UTR3_EditedmiRNA/" $miRNATmpFold"ProcessUTR3.EditedmiRNA"
cat $miRNATmpFold"ProcessUTR3.EditedmiRNA.targetscan.prediction" >> $miRNAFileFold"ProcessUTR3.EditedmiRNA.targetscan.prediction"
cat $miRNATmpFold"ProcessUTR3.EditedmiRNA.miRanda.prediction" >> $miRNAFileFold"ProcessUTR3.EditedmiRNA.miRanda.prediction"
break
fi
done
}

function BasemiRNAprocessing(){
DownloadFold=$DatabaseFold"RawData/"
DownloadReference=$DownloadFold"Reference/"
miRNAFold=$DatabaseFold"miRNA/"
miRNAFileFold=$miRNAFold"Files/"
miRNATmpFold=$miRNAFold"Tmp/"
LncRNAMerge $miRNATmpFold"lncRNA" $miRNATmpFold"lncRNA/" $miRNAFileFold"ProcesslncRNA"
cp $miRNATmpFold"lncRNA.targetscan.prediction" $miRNAFileFold
cp $miRNATmpFold"lncRNA.miRanda.prediction" $miRNAFileFold
UTR3Merge $miRNATmpFold"UTR3" $miRNATmpFold"UTR3/" $miRNAFileFold"ProcessUTR3"
cp $miRNATmpFold"UTR3.targetscan.prediction" $miRNAFileFold
cp $miRNATmpFold"UTR3.miRanda.prediction" $miRNAFileFold
} 

function EditingmiRNAprocess(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
DownloadFold=$DatabaseFold"RawData/"
DownloadReference=$DownloadFold"Reference/"
miRNAFold=$DatabaseFold"miRNA/"
miRNAFileFold=$miRNAFold"Files/"
miRNATmpFold=$miRNAFold"Tmp/"
cut -f1,5 $GeneAnnotationFold""$cancer".Table" | tail -n +2 | sort | uniq > $miRNATmpFold""$cancer".Table" 
$miRNAFileFold
}


function Compete(){
DownloadFold=$DatabaseFold"RawData/"
DownloadReference=$DownloadFold"Reference/"
miRNAFold=$DatabaseFold"miRNA/"
miRNAFileFold=$miRNAFold"Files/"

#GainLoss
Rscript miRTarMerge.R $miRNAFileFold"ProcesslncRNA" $miRNAFileFold"ProcessRNAediting.lncRNA" $miRNAFileFold"Final.lnRNA" "ENST" "TRUE" $miRNAFileFold"Final.lncRNA.base"
Rscript miRTarMerge.R $miRNAFileFold"ProcesslncRNA" $miRNAFileFold"ProcesslncRNA.EditedmiRNA" $miRNAFileFold"Final.miRNA.lncRNA" "miRNA" "FALSE"
Rscript miRTarMerge.R $miRNAFileFold"ProcessUTR3" $miRNAFileFold"ProcessRNAediting.UTR3" $miRNAFileFold"Final.UTR3" "ENST" "TRUE" $miRNAFileFold"Final.UTR3.base"
Rscript miRTarMerge.R $miRNAFileFold"ProcessUTR3" $miRNAFileFold"ProcessUTR3.EditedmiRNA" $miRNAFileFold"Final.miRNA.UTR3" "miRNA" "FALSE"
#Filter some cases in Final.UTR3.base due to the incorrect UTR prediction in the base file
awk -F"\t" 'ARGIND==1{UTR3[$1]=$1}\
            ARGIND==2{if(FNR==1){print $0}\
                      else{split($1,a,"_");out=a[1];for(i=2;i<=length(a)-1;i++){out=out"_"a[i]};if(out in UTR3){print $0}}}'\
            ../RawData/Reference/gencode.v22.UTR3.Really\
            $miRNAFileFold"Final.UTR3.base" > $miRNAFileFold"Final.UTR3.base.1"
mv $miRNAFileFold"Final.UTR3.base.1" $miRNAFileFold"Final.UTR3.base"
#mRNA-miRNA-lncRNA
echo -e "CAeditomeID\tRNAediting\tmRNAENSG\tmRNAENST\tmiRNA\tlncRNAENSG\tlncRNAENST\tGainLoss\tEditedRNA" > $miRNAFileFold"mRNA_miRNA_lncRNA" 
awk -F"\t" 'ARGIND==1 && ($3=="transcript"){gsub(/"| /,"",$9);split($9,a,";");ENSG=substr(a[1],length("gene_id")+1);ENST=substr(a[2],length("transcript_id")+1);Gene[ENST]=ENSG;\
                               split(ENST,b,".");Gene[b[1]]=ENSG}\
            ARGIND==2 && (FNR>1){split($1,a,"_");if(!(UTRBase[$2]~a[1])){UTRBase[$2]=UTRBase[$2]";"a[1]}}\
            ARGIND==3 && (FNR>1){CAeditomeID[$2]=$1}\
            ARGIND==4 && (FNR>1){split($2,a,"|");split($1,b,"|");split(UTRBase[a[2]],c,";");for(i=2;i<=length(c);i++){print CAeditomeID[a[1]]"\t"a[1]"\t"Gene[c[i]]"\t"c[i]"\t"a[2]"\t"Gene[b[1]]"\t"b[1]"\tGain\tmiRNA"}}\
            ARGIND==5 && (FNR>1){split($2,a,"|");split($1,b,"|");split(UTRBase[a[2]],c,";");for(i=2;i<=length(c);i++){print CAeditomeID[a[1]]"\t"a[1]"\t"Gene[c[i]]"\t"c[i]"\t"a[2]"\t"Gene[b[1]]"\t"b[1]"\tLoss\tmiRNA"}}'\
            $DownloadReference"gencode.v22.annotation.gtf"\
            $miRNAFileFold"Final.UTR3.base"\
            $DatabaseFold"GeneralAnnotation/GeneralInformation"\
            $miRNAFileFold"Final.miRNA.lncRNA.Gain"\
            $miRNAFileFold"Final.miRNA.lncRNA.Loss" >> $miRNAFileFold"mRNA_miRNA_lncRNA"
awk -F"\t" 'ARGIND==1 && ($3=="transcript"){gsub(/"| /,"",$9);split($9,a,";");ENSG=substr(a[1],length("gene_id")+1);ENST=substr(a[2],length("transcript_id")+1);Gene[ENST]=ENSG;\
                               split(ENST,b,".");Gene[b[1]]=ENSG}\
            ARGIND==2 && (FNR>1){split($1,a,"|");if(!(LncBase[$2]~a[1])){LncBase[$2]=LncBase[$2]";"a[1]}}\
            ARGIND==3 && (FNR>1){CAeditomeID[$2]=$1}\
            ARGIND==4 && (FNR>1){split($2,a,"|");split($1,b,"_");split(LncBase[a[2]],c,";");for(i=2;i<=length(c);i++){print CAeditomeID[a[1]]"\t"a[1]"\t"Gene[b[1]]"\t"b[1]"\t"a[2]"\t"Gene[c[i]]"\t"c[i]"\tGain\tmiRNA"}}\
            ARGIND==5 && (FNR>1){split($2,a,"|");split($1,b,"_");split(LncBase[a[2]],c,";");for(i=2;i<=length(c);i++){print CAeditomeID[a[1]]"\t"a[1]"\t"Gene[b[1]]"\t"b[1]"\t"a[2]"\t"Gene[c[i]]"\t"c[i]"\tLoss\tmiRNA"}}'\
            $DownloadReference"gencode.v22.annotation.gtf"\
            $miRNAFileFold"Final.lncRNA.base"\
            $DatabaseFold"GeneralAnnotation/GeneralInformation"\
            $miRNAFileFold"Final.miRNA.UTR3.Gain"\
            $miRNAFileFold"Final.miRNA.UTR3.Loss" >> $miRNAFileFold"mRNA_miRNA_lncRNA"
awk -F"\t" 'ARGIND==1 && ($3=="transcript"){gsub(/"| /,"",$9);split($9,a,";");ENSG=substr(a[1],length("gene_id")+1);ENST=substr(a[2],length("transcript_id")+1);Gene[ENST]=ENSG;\
                               split(ENST,b,".");Gene[b[1]]=ENSG}\
            ARGIND==2 && (FNR>1){split($1,a,"_");if(!(UTRBase[$2]~a[1])){UTRBase[$2]=UTRBase[$2]";"a[1]}}\
            ARGIND==3 && (FNR>1){CAeditomeID[$2]=$1}\
            ARGIND==4 && (FNR>1){split($1,a,"|");split(UTRBase[$2],c,";");for(i=2;i<=length(c);i++){print CAeditomeID[a[1]]"\t"a[1]"\t"Gene[c[i]]"\t"c[i]"\t"$2"\t"Gene[a[2]]"\t"a[2]"\tGain\tlncRNA"}}\
            ARGIND==5 && (FNR>1){split($1,a,"|");split(UTRBase[$2],c,";");for(i=2;i<=length(c);i++){print CAeditomeID[a[1]]"\t"a[1]"\t"Gene[c[i]]"\t"c[i]"\t"$2"\t"Gene[a[2]]"\t"a[2]"\tLoss\tlncRNA"}}'\
            $DownloadReference"gencode.v22.annotation.gtf"\
            $miRNAFileFold"Final.UTR3.base"\
            $DatabaseFold"GeneralAnnotation/GeneralInformation"\
            $miRNAFileFold"Final.lncRNA.Gain"\
            $miRNAFileFold"Final.lncRNA.Loss" >> $miRNAFileFold"mRNA_miRNA_lncRNA"
awk -F"\t" 'ARGIND==1 && ($3=="transcript"){gsub(/"| /,"",$9);split($9,a,";");ENSG=substr(a[1],length("gene_id")+1);ENST=substr(a[2],length("transcript_id")+1);Gene[ENST]=ENSG;\
                               split(ENST,b,".");Gene[b[1]]=ENSG}\
            ARGIND==2 && (FNR>1){split($1,a,"|");if(!(UTRBase[$2]~a[1])){UTRBase[$2]=UTRBase[$2]";"a[1]}}\
            ARGIND==3 && (FNR>1){CAeditomeID[$2]=$1}\
            ARGIND==4 && (FNR>1){split($1,a,"|");split(a[2],b,"_");split(UTRBase[$2],c,";");for(i=2;i<=length(c);i++){print CAeditomeID[a[1]]"\t"a[1]"\t"Gene[b[1]]"\t"b[1]"\t"$2"\t"Gene[c[i]]"\t"c[i]"\tGain\tmRNA"}}\
            ARGIND==5 && (FNR>1){split($1,a,"|");split(a[2],b,"_");split(UTRBase[$2],c,";");for(i=2;i<=length(c);i++){print CAeditomeID[a[1]]"\t"a[1]"\t"Gene[b[1]]"\t"b[1]"\t"$2"\t"Gene[c[i]]"\t"c[i]"\tLoss\tmRNA"}}'\
            $DownloadReference"gencode.v22.annotation.gtf"\
            $miRNAFileFold"Final.lncRNA.base"\
            $DatabaseFold"GeneralAnnotation/GeneralInformation"\
            $miRNAFileFold"Final.UTR3.Gain"\
            $miRNAFileFold"Final.UTR3.Loss" >> $miRNAFileFold"mRNA_miRNA_lncRNA"
}

function CompeteExpression(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
index=$2
Flag=$3
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
DownloadFold=$DatabaseFold"RawData/"
DownloadExpression=$DownloadFold"GeneExpression/"
DownloadRNAediting=$DownloadFold"RNAediting/"
miRNAFold=$DatabaseFold"miRNA/"
miRNAFileFold=$miRNAFold"Files/"
miRNATmpFold=$miRNAFold"TmpExpression/"
mkdir $miRNATmpFold
#Expression
awk -F"," 'ARGIND==1 && (FNR==1){for(i=2;i<=NF;i++){split($i,a,"-");if(substr(a[4],1,1)==0){sam=a[1]"-"a[2]"-"a[3]"-"a[4];Expression[sam]=Expression[sam]";"i}}}\
           ARGIND==2 && (FNR==1){split($0,a,"\t");for(i=2;i<=length(a);i++){if(a[i] in Expression){print a[i]"\t"i"\t"substr(Expression[a[i]],2)}}}'\
           $DownloadExpression""""$cancer"_expression.csv"\
           $DownloadRNAediting""$cancer"_Editing.txt"  > $miRNATmpFold""$cancer".editing.expression.samples"
RNAeditingColumn=`cut -f2 $miRNATmpFold""$cancer".editing.expression.samples"`
ExpressionColumn=`cut -f3 $miRNATmpFold""$cancer".editing.expression.samples"`
awk -F"\t" -v RNAeditingColumn="${RNAeditingColumn[*]}" '{split(RNAeditingColumn,col," ");out=$1;for(i=1;i<=length(col);i++){out=out"\t"$col[i]};print out}'\
            $DownloadRNAediting""$cancer"_Editing.txt" > $miRNATmpFold""$cancer".RNAediting"
awk -F"," -v ExpressionColumn="${ExpressionColumn[*]}" '{split(ExpressionColumn,col," ");\
                                                         if(NR==1){out="Gene";for(i=1;i<=length(col);i++){split(col[i],b,";");out=out"\t"$b[1]};print out}\
                                                         else{out=$1;for(i=1;i<=length(col);i++){\
                                                         split(col[i],b,";");m=0;for(j=1;j<=length(b);j++){m=m+$b[j]};out=out"\t"m/length(b)};print out}}'\
            $DownloadExpression""$cancer"_expression.csv" > $miRNATmpFold""$cancer".Expression"
awk -F"\t" 'ARGIND==1 && (FNR>1){Editing[$1]=$1}\
            ARGIND==2 && (FNR>1){Gene[$1]=$1}\
            ARGIND==3 && (FNR>1) && ($2 in Editing){split($3,a,".");split($6,b,".");if((a[1] in Gene) && (b[1] in Gene)){print $2"\t"a[1]"\t"b[1]}}'\
            $GeneAnnotationFold""$cancer".Table"\
            $miRNATmpFold""$cancer".Expression"\
            $miRNAFileFold"mRNA_miRNA_lncRNA_"$index | sort | uniq > $miRNATmpFold""$cancer".mRNA_miRNA_lncRNA_"$index
awk -F"\t" 'ARGIND==1{Gene[$1]=$0}\
            ARGIND==2{print Gene[$2]}'\
            $miRNATmpFold""$cancer".Expression"\
            $miRNATmpFold""$cancer".mRNA_miRNA_lncRNA_"$index > $miRNATmpFold""$cancer".Expression.mRNA"
awk -F"\t" 'ARGIND==1{Gene[$1]=$0}\
            ARGIND==2{print Gene[$3]}'\
            $miRNATmpFold""$cancer".Expression"\
            $miRNATmpFold""$cancer".mRNA_miRNA_lncRNA_"$index > $miRNATmpFold""$cancer".Expression.lncRNA"
awk -F"\t" 'ARGIND==1{RNAediting[$1]=$0}\
            ARGIND==2 {print RNAediting[$1]}'\
            $miRNATmpFold""$cancer".RNAediting"\
            $miRNATmpFold""$cancer".mRNA_miRNA_lncRNA_"$index > $miRNATmpFold""$cancer".Expression.RNAediting"
mkdir $miRNATmpFold""$cancer
split -l 50000 $miRNATmpFold""$cancer".Expression.mRNA" -d $miRNATmpFold""$cancer"/"$cancer".Expression.mRNA_"
split -l 50000 $miRNATmpFold""$cancer".Expression.lncRNA" -d $miRNATmpFold""$cancer"/"$cancer".Expression.lncRNA_"
split -l 50000 $miRNATmpFold""$cancer".Expression.RNAediting" -d $miRNATmpFold""$cancer"/"$cancer".Expression.RNAediting_"
rm -f $miRNATmpFold""$cancer".Expression."*
if [[ $Flag == "TRUE" ]]
then
for file in `ls $miRNATmpFold""$cancer"/"$cancer".Expression.mRNA_"*`
do
while true
do
num=`ps -aux | grep "CompeteCompute" | wc -l`
if [[ $num -lt 100 ]]
then
Rscript CompeteCompute.R ${file/"Expression.mRNA"/"Expression.RNAediting"} $file ${file/"Expression.mRNA"/"Expression.lncRNA"} ${file/"Expression.mRNA"/"compete.results"} &
break
fi
done
done
fi
}

function CompeteMerge(){
CompeteResultsFile=$1
CompeteReferenceFile=$2
EditedFlag=$3
OutputFile=$4
head -n 1 $CompeteResultsFile"_00" > $CompeteResultsFile
for file in `ls $CompeteResultsFile"_"*`
do
tail -n +2 $file >> $CompeteResultsFile
done
head -n 1 $CompeteReferenceFile | awk -F"\t" '{print $0"\tmRNATtestP\tmRNAMeanExpressionEdited\tmRNAMeanExpressionNonedited\tmRNAlogFC\tmRNACorreP\tmRNACorreR\tlncRNATtestP\tlncRNAMeanExpressionEdited\tlncRNAMeanExpressionNonedited\tlncRNAlogFC\tlncRNACorreP\tlncRNACorreR"}' > $OutputFile
if [[ $EditedFlag == "mRNA" ]]
then
awk -F"\t" 'ARGIND==1 && (FNR>1){Known[$1"\t"$2"\t"$9]=$0}\
            ARGIND==2 && (FNR>1){split($3,a,".");split($6,b,".");if($2"\t"a[1]"\t"b[1] in Known){split(Known[$2"\t"a[1]"\t"b[1]],c,"\t");\
                      if(c[3]<0.05 && c[7]<0.05 && c[10]<0.05 && c[14]<0.05 && c[6]*c[13]<0 && c[8]*c[15]<0 && \
                         (($8=="Gain" && c[6]<0 && c[8]<0) || ($8=="Loss" && c[6]>0 && c[8]>0))){print $0"\t"c[3]"\t"c[4]"\t"c[5]"\t"c[6]"\t"c[7]"\t"c[8]"\t"\
                                                                                                                 c[10]"\t"c[11]"\t"c[12]"\t"c[13]"\t"c[14]"\t"c[15]}}}'\
            $CompeteResultsFile\
            $CompeteReferenceFile >> $OutputFile
elif [[ $EditedFlag == "lncRNA" ]]
then
awk -F"\t" 'ARGIND==1 && (FNR>1){Known[$1"\t"$2"\t"$9]=$0}\
            ARGIND==2 && (FNR>1){split($3,a,".");split($6,b,".");if($2"\t"a[1]"\t"b[1] in Known){split(Known[$2"\t"a[1]"\t"b[1]],c,"\t");\
                      if(c[3]<0.05 && c[7]<0.05 && c[10]<0.05 && c[14]<0.05 && c[6]*c[13]<0 && c[8]*c[15]<0 && \
                         (($8=="Gain" && c[13]<0 && c[15]<0) || ($8=="Loss" && c[13]>0 && c[15]>0))){print $0"\t"c[3]"\t"c[4]"\t"c[5]"\t"c[6]"\t"c[7]"\t"c[8]"\t"\
                                                                                                                 c[10]"\t"c[11]"\t"c[12]"\t"c[13]"\t"c[14]"\t"c[15]}}}'\
            $CompeteResultsFile\
            $CompeteReferenceFile >> $OutputFile
fi
}

function BindMerge(){
FileFold=$1
GeneralFile=$DatabaseFold"GeneralAnnotation/GeneralInformation"
ReferenceFile=$DatabaseFold"RawData/Reference/gencode.v22.annotation.gtf"
echo -e "CAeditomeID\tRNAediting\tENST\tENSG\tGeneName\tmiRNA\tTargetStart\tTargetTerminal\tTargetScore\tmiRandaStart\tmiRandaTerminal\tmiRandaScore\tmiRandaEnergyScore\tEditedRNA_BindRNA\tGainLoss" > $FileFold"RNAeditingMiRNABindingInformation"
awk -F"\t" 'ARGIND==1{Editing[$2]=$1}\
            ARGIND==2 && ($3=="transcript"){gsub(/"| /,"",$9);split($9,a,";");ENSG=substr(a[1],length("gene_id")+1);ENST=substr(a[2],length("transcript_id")+1);\
                                            GeneName=substr(a[5],length("gene_name")+1);Gene[ENST]=ENSG"\t"GeneName;\
                                            split(ENST,b,".");Gene[b[1]]=ENSG"\t"GeneName}\
            ARGIND==3 && (FNR>1){split($1,a,"|");out=Editing[a[1]]"\t"a[1]"\t"a[2]"\t"Gene[a[2]];for(i=2;i<=NF;i++){out=out"\t"$i};out=out"\tEditedLncRNA\tGain";print out}\
            ARGIND==4 && (FNR>1){split($1,a,"|");out=Editing[a[1]]"\t"a[1]"\t"a[2]"\t"Gene[a[2]];for(i=2;i<=NF;i++){out=out"\t"$i};out=out"\tEditedLncRNA\tLoss";print out}\
            ARGIND==5 && (FNR>1){split($1,a,"|");split(a[2],b,"_");out=Editing[a[1]]"\t"a[1]"\t"b[1]"\t"Gene[b[1]];for(i=2;i<=NF;i++){out=out"\t"$i};\
                                 out=out"\tEditedmRNA\tGain";print out}\
            ARGIND==6 && (FNR>1){split($1,a,"|");split(a[2],b,"_");out=Editing[a[1]]"\t"a[1]"\t"b[1]"\t"Gene[b[1]];for(i=2;i<=NF;i++){out=out"\t"$i};\
                                 out=out"\tEditedmRNA\tLoss";print out}\
            ARGIND==7 && (FNR>1){split($1,a,"|");split($2,b,"|");out=Editing[b[1]]"\t"b[1]"\t"a[1]"\t"Gene[a[1]]"\t"b[2];\
                                 for(i=3;i<=NF;i++){out=out"\t"$i};out=out"\tEditedmiRNA_BindlncRNA\tGain";print out}\
            ARGIND==8 && (FNR>1){split($1,a,"|");split($2,b,"|");out=Editing[b[1]]"\t"b[1]"\t"a[1]"\t"Gene[a[1]]"\t"b[2];\
                                 for(i=3;i<=NF;i++){out=out"\t"$i};out=out"\tEditedmiRNA_BindlncRNA\tLoss";print out}\
            ARGIND==9 && (FNR>1){split($1,a,"_");split($2,b,"|");out=Editing[b[1]]"\t"b[1]"\t"a[1]"\t"Gene[a[1]]"\t"b[2];\
                                   for(i=3;i<=NF;i++){out=out"\t"$i};out=out"\tEditedmiRNA_BindmRNA\tGain";print out}\
            ARGIND==10 && (FNR>1){split($1,a,"_");split($2,b,"|");out=Editing[b[1]]"\t"b[1]"\t"a[1]"\t"Gene[a[1]]"\t"b[2];\
                                   for(i=3;i<=NF;i++){out=out"\t"$i};out=out"\tEditedmiRNA_BindmRNA\tLoss";print out}'\
            $GeneralFile\
            $ReferenceFile\
            $FileFold"Final.lncRNA.Gain"\
            $FileFold"Final.lncRNA.Loss"\
            $FileFold"Final.UTR3.Gain"\
            $FileFold"Final.UTR3.Loss"\
            $FileFold"Final.miRNA.lncRNA.Gain"\
            $FileFold"Final.miRNA.lncRNA.Loss"\
            $FileFold"Final.miRNA.UTR3.Gain"\
            $FileFold"Final.miRNA.UTR3.Loss" >> $FileFold"RNAeditingMiRNABindingInformation"
#Filter some information because of the incorrect code for UTR3 download
awk -F"\t" 'ARGIND==1{UTR3[$1]=$1}\
            ARGIND==2{if($14=="EditedmRNA" || $14=="EditedmiRNA_BindmRNA"){if($3 in UTR3){print $0}}\
                          else{print $0}}'\
            ../RawData/Reference/gencode.v22.UTR3.Really\
            $FileFold"RNAeditingMiRNABindingInformation" > $FileFold"RNAeditingMiRNABindingInformation.1"
mv $FileFold"RNAeditingMiRNABindingInformation.1" $FileFold"RNAeditingMiRNABindingInformation"
            
}

function CompeteExpression2(){
CancerIndex=$1
Merge=$2
cancer=${cancertypes[$CancerIndex]}
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
DownloadFold=$DatabaseFold"RawData/"
DownloadExpression=$DownloadFold"GeneExpression/"
DownloadRNAediting=$DownloadFold"RNAediting/"
FileFold="/data9/swu13/CAeditome_hg38/miRNA/Files/"
GainLossFile=$FileFold"RNAeditingMiRNABindingInformation"
BaseLncRNAFile=$FileFold"Final.lncRNA.base"
BaseMRNAFile=$FileFold"Final.UTR3.base.1"
miRNATmpFold="/data9/swu13/CAeditome_hg38/miRNA/TmpMergeExpression/"
if [[ $Merge != "TRUE" ]]
then
#Expression
awk -F"," 'ARGIND==1 && (FNR==1){for(i=2;i<=NF;i++){split($i,a,"-");if(substr(a[4],1,1)==0){sam=a[1]"-"a[2]"-"a[3]"-"a[4];Expression[sam]=Expression[sam]";"i}}}\
           ARGIND==2 && (FNR==1){split($0,a,"\t");for(i=2;i<=length(a);i++){if(a[i] in Expression){print a[i]"\t"i"\t"substr(Expression[a[i]],2)}}}'\
           $DownloadExpression""""$cancer"_expression.csv"\
           $DownloadRNAediting""$cancer"_Editing.txt"  > $miRNATmpFold""$cancer".editing.expression.samples"
RNAeditingColumn=`cut -f2 $miRNATmpFold""$cancer".editing.expression.samples"`
ExpressionColumn=`cut -f3 $miRNATmpFold""$cancer".editing.expression.samples"`
awk -F"\t" -v RNAeditingColumn="${RNAeditingColumn[*]}" '{split(RNAeditingColumn,col," ");out=$1;for(i=1;i<=length(col);i++){out=out"\t"$col[i]};print out}'\
            $DownloadRNAediting""$cancer"_Editing.txt" > $miRNATmpFold""$cancer".RNAediting"
awk -F"," -v ExpressionColumn="${ExpressionColumn[*]}" '{split(ExpressionColumn,col," ");\
                                                         if(NR==1){out="Gene";for(i=1;i<=length(col);i++){split(col[i],b,";");out=out"\t"$b[1]};print out}\
                                                         else{out=$1;for(i=1;i<=length(col);i++){\
                                                         split(col[i],b,";");m=0;for(j=1;j<=length(b);j++){m=m+$b[j]};out=out"\t"m/length(b)};print out}}'\
            $DownloadExpression""$cancer"_expression.csv" > $miRNATmpFold""$cancer".Expression"
#pair
awk -F"\t" 'ARGIND==1 && (FNR>1){Editing[$1]=$1}\
            ARGIND==2 && (FNR>1){Gene[$1]=$1}\
            ARGIND==3 && (FNR>1) && ($2 in Editing){split($4,a,".");if((a[1] in Gene)){print $2"\t"a[1]}}'\
            $GeneAnnotationFold""$cancer".Table"\
            $miRNATmpFold""$cancer".Expression"\
            $GainLossFile | sort | uniq > $miRNATmpFold""$cancer"Edited"
#Compute correlations
awk -F"\t" 'ARGIND==1{Gene[$1]=$0}\
            ARGIND==2{print Gene[$2]}'\
            $miRNATmpFold""$cancer".Expression"\
            $miRNATmpFold""$cancer"Edited" > $miRNATmpFold""$cancer".Expression.RNA"
awk -F"\t" 'ARGIND==1{RNAediting[$1]=$0}\
            ARGIND==2 {print RNAediting[$1]}'\
            $miRNATmpFold""$cancer".RNAediting"\
            $miRNATmpFold""$cancer"Edited" > $miRNATmpFold""$cancer".Expression.RNAediting"
Rscript CompeteComputeTwo.R  $miRNATmpFold""$cancer".Expression.RNAediting" $miRNATmpFold""$cancer".Expression.RNA" $miRNATmpFold""$cancer"Edited.results"
awk -F"\t" 'ARGIND==1 && ($3<0.05) && ($7<0.05){Cor[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}\
            ARGIND==2{split($4,a,".");if($2"\t"a[1] in Cor){split(Cor[$2"\t"a[1]],b,"\t");\
                                      if(($NF=="Gain" && b[2]<b[3] && b[6]<0) || ($NF=="Loss" && b[2]>b[3] && b[6]>0)){
                                      out=$1"\t"$2"\t"$4;for(i=6;i<=NF;i++){out=out"\t"$i};print out"\t"Cor[$2"\t"a[1]]}}}'\
            $miRNATmpFold""$cancer"Edited.results"\
            $GainLossFile | sort | uniq > $miRNATmpFold""$cancer"Edited.significant"
###Base
awk -F"\t" 'ARGIND==1 && (FNR>1){Ed[$1]=$1}\
            ARGIND==2 && ($3=="transcript"){gsub(/"| /,"",$9);split($9,a,";");ENSG=substr(a[1],length("gene_id")+1);ENST=substr(a[2],length("transcript_id")+1);Gene[ENST]=ENSG;\
                                            split(ENST,b,".");Gene[b[1]]=ENSG}\
            ARGIND==3 && (FNR>1){Expression[$1]=$1}\
            ARGIND==4{if($12=="EditedmRNA"){if(!(LncEditing[$4]~$2)){LncEditing[$4]=LncEditing[$4]";"$2}}\
                      if($12=="EditedLncRNA"){if(!(mEditing[$4]~$2)){mEditing[$4]=mEditing[$4]";"$2}}}\
            ARGIND==5 && (FNR>1) && ($2 in LncEditing){split($1,a,"|");split(Gene[a[1]],b,".");split(LncEditing[$2],c,";");\
                                                  if(b[1] in Expression){for(i=2;i<=length(c);i++){if(c[i] in Ed){print c[i]"\t"b[1]"\t"$2}}}}\
            ARGIND==6 && (FNR>1) && ($2 in mEditing){split($1,a,"_");split(Gene[a[1]],b,".");split(mEditing[$2],c,";");\
                                                  if(b[1] in Expression){for(i=2;i<=length(c);i++){if(c[i] in Ed){print c[i]"\t"b[1]"\t"$2}}}}'\
            $GeneAnnotationFold""$cancer".Table"\
            $DownloadFold"Reference/gencode.v22.annotation.gtf"\
            $miRNATmpFold""$cancer".Expression"\
            $miRNATmpFold""$cancer"Edited.significant"\
            $BaseLncRNAFile\
            $BaseMRNAFile > $miRNATmpFold""$cancer"Base"
cut -f1-2 $miRNATmpFold""$cancer"Base" | sort | uniq > $miRNATmpFold""$cancer"Pair"
#Compute correlations
awk -F"\t" 'ARGIND==1{Gene[$1]=$0}\
            ARGIND==2{print Gene[$2]}'\
            $miRNATmpFold""$cancer".Expression"\
            $miRNATmpFold""$cancer"Pair" > $miRNATmpFold""$cancer".Expression.RNA"
awk -F"\t" 'ARGIND==1{RNAediting[$1]=$0}\
            ARGIND==2 {print RNAediting[$1]}'\
            $miRNATmpFold""$cancer".RNAediting"\
            $miRNATmpFold""$cancer"Pair" > $miRNATmpFold""$cancer".Expression.RNAediting"
Rscript CompeteComputeTwo.R  $miRNATmpFold""$cancer".Expression.RNAediting" $miRNATmpFold""$cancer".Expression.RNA" $miRNATmpFold""$cancer"Base.results"
else
echo -e "Cancer\tCAeditome\tRNAediting\tENSG\tGene\tEditedType\tGainLoss\tmiRNA\tCompeteENSG\tCompeteGene\tCompeteEditedType\tCompeteGainLoss\tTtestP1\tMeanExpressionEdited1\tMeanExpressionNonedited1\tlogFC1\tCorreP1\tCorreR1\tTtestP2\tMeanExpressionEdited2\tMeanExpressionNonedited2\tlogFC2\tCorreP2\tCorreR2" > $miRNATmpFold"miRNARegulation"
for((i=0;i<=32;i++))
do
cancer=${cancertypes[$i]}
awk -F"\t" -v cancer="$cancer"\
           'ARGIND==1 && ($3<0.05) && ($7<0.05) && ($6*$8>0) && (FNR>1){Base[$1"\t"$2]=$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}\
            ARGIND==2 && ($1"\t"$2 in Base){if(!(CorrGene[$1"\t"$3] ~$2)){CorrGene[$1"\t"$3]=CorrGene[$1"\t"$3]";"$2}}\
            ARGIND==3 && ($3=="gene"){gsub(/"| /,"",$9);split($9,a,";");ENSG=substr(a[1],length("gene_id")+1);GeneName=substr(a[4],length("gene_name")+1);Gene[ENSG]=GeneName;\
                                            split(ENSG,b,".");Gene[b[1]]=GeneName}\
            ARGIND==4 && ($12=="EditedmRNA" || $12=="EditedLncRNA"){split(CorrGene[$2"\t"$4],a,";");for(i=2;i<=length(a);i++){\
                                                         if($2"\t"a[i] in Base){split(Base[$2"\t"a[i]],c,"\t");\
                                                             if(($13=="Gain" && c[4]>0 && c[6]>0) || ($13=="Loss" && c[4]<0 && c[6]<0)){out=cancer"\t"$1"\t"$2"\t"$3"\t"Gene[$3]"\t"$12"\t"$13"\t"$4"\t"a[i]"\t"Gene[a[i]]"\t.\t.";\
                                                             for(j=14;j<=NF;j++){out=out"\t"$j};print out"\t"Base[$2"\t"a[i]]}}}}'\
            $miRNATmpFold""$cancer"Base.results"\
            $miRNATmpFold""$cancer"Base"\
            $DownloadFold"Reference/gencode.v22.annotation.gtf"\
            $miRNATmpFold""$cancer"Edited.significant" >> $miRNATmpFold"miRNARegulation"
awk -F"\t" -v cancer="$cancer" \
           'ARGIND==1 && ($3=="gene"){gsub(/"| /,"",$9);split($9,a,";");ENSG=substr(a[1],length("gene_id")+1);GeneName=substr(a[4],length("gene_name")+1);\
                                      Gene[ENSG]=GeneName;split(ENSG,b,".");Gene[b[1]]=GeneName}\
            ARGIND==2{if($12=="EditedmiRNA_BindmRNA"){Editing[$1"\t"$4]=$1"\t"$4;\
                                                      if($13=="Gain"){mRNAGain[$1"\t"$4]=mRNAGain[$1"\t"$4]";"$0}\
                                                      else if($13=="Loss"){mRNALoss[$1"\t"$4]=mRNALoss[$1"\t"$4]";"$0}}\
                      else if($12=="EditedmiRNA_BindlncRNA"){Editing[$1"\t"$4]=$1"\t"$4;\
                                                      if($13=="Gain"){lncRNAGain[$1"\t"$4]=lncRNAGain[$1"\t"$4]";"$0}\
                                                      else if($13=="Loss"){lncRNALoss[$1"\t"$4]=lncRNALoss[$1"\t"$4]";"$0}}}\
             END{for(i in Editing){if((i in lncRNALoss) && (i in mRNAGain)){split(lncRNALoss[i],a,";");split(mRNAGain[i],b,";");\
                                              for(i=2;i<=length(a);i++){for(j=2;j<=length(b);j++){split(b[j],m,"\t");split(a[i],n,"\t");\
                                              out=cancer"\t"m[1]"\t"m[2]"\t"m[3]"\t"Gene[m[3]]"\t"m[12]"\t"m[13]"\t"m[4]"\t"n[3]"\t"Gene[n[3]]"\t"n[12]"\t"n[13];\
                                              for(p=14;p<=19;p++){out=out"\t"m[p]};\
                                              for(p=14;p<=19;p++){out=out"\t"n[p]};print out}}}\
                             else if((i in mRNALoss) && (i in lncRNAGain)){split(lncRNAGain[i],a,";");split(mRNALoss[i],b,";");\
                                              for(i=2;i<=length(a);i++){for(j=2;j<=length(b);j++){split(b[j],m,"\t");split(a[i],n,"\t");\
                                              out=cancer"\t"m[1]"\t"m[2]"\t"m[3]"\t"Gene[m[3]]"\t"m[12]"\t"m[13]"\t"m[4]"\t"n[3]"\t"Gene[n[3]]"\t"n[12]"\t"n[13];\
                                              for(p=14;p<=19;p++){out=out"\t"m[p]};\
                                              for(p=14;p<=19;p++){out=out"\t"n[p]};print out}}}}}'\
            $DownloadFold"Reference/gencode.v22.annotation.gtf"\
            $miRNATmpFold""$cancer"Edited.significant" >> $miRNATmpFold"miRNARegulation"
done
fi
}

function RegulationNetwork(){
Flag=$1
miRNAFold=$DatabaseFold"miRNA/"
RegulationFold=$miRNAFold"Files/"
NetworkFold=$miRNAFold"Network"$Flag"/"
mkdir $NetworkFold
mkdir $miRNAFold"TmpNetwork/"
#
if [[ $Flag == "All" ]]
then
rm -f $miRNAFold"TmpNetwork/Node.Data"
rm -f $miRNAFold"TmpNetwork/Edge.Data"
awk -F"\t" '(NR>1) && (!($6~"miRNA")){print $4"\t"$1"\t"$5"("$2")\tpurple\t8";\
                                      print $4"\t"$1"\t"$10"\tnavy\t8";\
                                      print $4"\t"$1"\t"$8"\tblack\t8"}'\
            $RegulationFold"miRNARegulation" >> $miRNAFold"TmpNetwork/Node.Data"
awk -F"\t" '(NR>1) && (($6~"miRNA")){print $8"\t"$1"\t"$5"\tpurple\t8";\
                                     print $8"\t"$1"\t"$10"\tpurple\t8";\
                                     print $8"\t"$1"\t"$8"("$2")\tblack\t8"}'\
            $RegulationFold"miRNARegulation" >> $miRNAFold"TmpNetwork/Node.Data"
awk -F"\t" '(NR>1) && (!($6~"miRNA")){if($7=="Gain"){col1="red";col2="green"}else{col1="green";col2="red"};\
                                      if($16>0){SA=$16}else{SA=0-$16};if($24>0){SB=$24}else{SB=0-$24};\
                                      print $4"\t"$1"\t"$5"("$2")\t"$8"\t"col1"\t"SA;\
                                      print $4"\t"$1"\t"$10"\t"$8"\t"col2"\t"SB}'\
            $RegulationFold"miRNARegulation" >> $miRNAFold"TmpNetwork/Edge.Data"           
awk -F"\t" '(NR>1) && (($6~"miRNA")){if($7=="Gain"){col1="red";col2="green"}else{col1="green";col2="red"};\
                                     if($16>0){SA=$16}else{SA=0-$16};if($24>0){SB=$24}else{SB=0-$24};\
                                     print $8"\t"$1"\t"$5"\t"$8"("$2")\t"col1"\t"SA;\
                                     print $8"\t"$1"\t"$10"\t"$8"("$2")\t"col2"\t"SB}'\
            $RegulationFold"miRNARegulation" >> $miRNAFold"TmpNetwork/Edge.Data" 
cat $miRNAFold"TmpNetwork/Node.Data" | sort | uniq > $miRNAFold"TmpNetwork/Node.Data.1"
mv $miRNAFold"TmpNetwork/Node.Data.1" $miRNAFold"TmpNetwork/Node.Data"
cat $miRNAFold"TmpNetwork/Edge.Data" | sort  |uniq > $miRNAFold"TmpNetwork/Edge.Data.1"
mv $miRNAFold"TmpNetwork/Edge.Data.1" $miRNAFold"TmpNetwork/Edge.Data"
Rscript Network.R $miRNAFold"TmpNetwork/Node.Data" $miRNAFold"TmpNetwork/Edge.Data" $NetworkFold
fi
if [[ $Flag == "Tumor" ]]
then
#
rm -f $miRNAFold"TmpNetwork/Node.Data"
rm -f $miRNAFold"TmpNetwork/Edge.Data"
awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
            ARGIND==2 && (FNR>1) && (!($6~"miRNA"))&& ($10 in Tumor){\
                                      print $4"\t"$1"\t"$5"("$2")\tpurple\t8";\
                                      print $4"\t"$1"\t"$10"\tnavy\t8";\
                                      print $4"\t"$1"\t"$8"\tblack\t8"}'\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            $RegulationFold"miRNARegulation" >> $miRNAFold"TmpNetwork/Node.Data"
awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
            ARGIND==2 && (FNR>1) && (($6~"miRNA")){if($5 in Tumor){print $8"\t"$1"\t"$5"\tpurple\t8"};\
                                                   if($10 in Tumor){print $8"\t"$1"\t"$10"\tpurple\t8"};\
                                                   if(($5 in Tumor) || ($10 in Tumor)){print $8"\t"$1"\t"$8"("$2")\tblack\t8"}}'\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            $RegulationFold"miRNARegulation" >> $miRNAFold"TmpNetwork/Node.Data"
awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
            ARGIND==2 && (FNR>1) && (!($6~"miRNA"))&& ($10 in Tumor){\
                                      if($7=="Gain"){col1="red";col2="green"}else{col1="green";col2="red"};\
                                      if($16>0){SA=$16}else{SA=0-$16};if($22>0){SB=$22}else{SB=0-$22};\
                                      print $4"\t"$1"\t"$5"("$2")\t"$8"\t"col1"\t"SA;\
                                      print $4"\t"$1"\t"$10"\t"$8"\t"col2"\t"SB}'\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            $RegulationFold"miRNARegulation" >> $miRNAFold"TmpNetwork/Edge.Data"           
awk -F"\t" 'ARGIND==1{Tumor[$1]=$1}\
            ARGIND==2 && (FNR>1) && (($6~"miRNA")){\
                                     if($7=="Gain"){col1="red";col2="green"}else{col1="green";col2="red"};\
                                     if($16>0){SA=$16}else{SA=0-$16};if($22>0){SB=$22}else{SB=0-$22};\
                                     if($5 in Tumor){print $8"\t"$1"\t"$5"\t"$8"("$2")\t"col1"\t"SA};\
                                     if($10 in Tumor){print $8"\t"$1"\t"$10"\t"$8"("$2")\t"col2"\t"SB}}'\
            /data9/swu13/CAeditome_hg38/RawData/TumorGene/Genes\
            $RegulationFold"miRNARegulation" >> $miRNAFold"TmpNetwork/Edge.Data" 
cat $miRNAFold"TmpNetwork/Node.Data" | sort | uniq > $miRNAFold"TmpNetwork/Node.Data.1"
mv $miRNAFold"TmpNetwork/Node.Data.1" $miRNAFold"TmpNetwork/Node.Data"
cat $miRNAFold"TmpNetwork/Edge.Data" | sort  |uniq > $miRNAFold"TmpNetwork/Edge.Data.1"
mv $miRNAFold"TmpNetwork/Edge.Data.1" $miRNAFold"TmpNetwork/Edge.Data"
Rscript Network.R $miRNAFold"TmpNetwork/Node.Data" $miRNAFold"TmpNetwork/Edge.Data" $NetworkFold
mv $NetworkFold"hsa-miR-1304-3p.png" $NetworkFold"ENSG00000221170.1.png"
mv $NetworkFold"hsa-miR-200b-3p.png" $NetworkFold"ENSG00000207730.2.png"
mv $NetworkFold"hsa-miR-27a-3p.png" $NetworkFold"ENSG00000207808.1.png"
mv $NetworkFold"hsa-miR-4477b.png" $NetworkFold"ENSG00000266017.1.png"
mv $NetworkFold"hsa-miR-548aa.png" $NetworkFold"ENSG00000207688.2.png"
mv $NetworkFold"hsa-miR-589-3p.png" $NetworkFold"ENSG00000207973.1.png"
mv $NetworkFold"hsa-miR-624-3p.png" $NetworkFold"ENSG00000207952.1.png"
fi
}


function main(){
index=$1
miRNABinding $index
}