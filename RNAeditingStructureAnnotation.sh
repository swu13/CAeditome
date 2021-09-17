cancertypes=("BLCA"  "LGG"  "BRCA"  "LIHC"  "CESC"  "LUAD"  "COAD"  "LUSC"  "GBM"  "PRAD"  "HNSC"  "STAD"  "KICH"  "THCA"  "KIRC"  "UCEC"  "KIRP"  "ACC"  "DLBC"  "ESCA"  "LAML"  "OV"  "PAAD"  "PCPG"  "SKCM"  "UCS"  "UVM"  "CHOL"  "MESO"  "READ"  "SARC"  "TGCT"  "THYM")
DatabaseFold="/data9/swu13/CAeditome_hg38/"

function MkdirPath(){
mkdir $DatabaseFold"Structure/"
}


function Transfer(){
File=$1
pa1=$File"_dp.ps"
pa2=$File"_ss.ps"
pa3=$File"_rss.ps"
ToolFold=$DatabaseFold"Tool/"
ViennaFold=$ToolFold"ViennaRNA-2.4.17/"
if [[ -f $pa1 ]] && [[ -f $pa2 ]] && [[ ! -f ${pa3/_rss.ps/_rss.png} ]]
then
perl $ViennaFold"src/Utils/relplot.pl" $pa2 $pa1 > $pa3
gmt psconvert $pa3 -Tg
convert -resize 1280 ${pa3/_rss.ps/_rss.png} ${pa3/_rss.ps/_rss_resized.png}
fi
#ResizedFile=${pa3/_rss.ps/_rss_resized.png}
#if [[ -f $ResizedFile ]]
#then
#rm -f $pa1 $pa2 $pa3 ${pa3/_rss.ps/_rss.png}
#fi
}

function RNAfoldPrediction(){
#Install RNAfold; GMT (conda config --add channels conda-forge,conda install gmt=6.0.0); imagemagick(conda install imagemagick;there will be confilts, but it is ok)
#ToolFold=$DatabaseFold"Tool/"
#wget -P $ToolFold https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.17.tar.gz
#tar -zxvf $ToolFold"ViennaRNA-2.4.17.tar.gz" -C $ToolFold
#cd $ToolFold"ViennaRNA-2.4.17/"
#export CPATH=/data1/anaconda2/pkgs/perl-5.32.0-h36c2ea0_0/lib/5.32.0/x86_64-linux-thread-multi/CORE/ # if the perl version is 5, you may need to find EXTERN.h first 
#./configure --prefix=$ToolFold"ViennaRNA-2.4.17/ViennaRNA/"
#make
#make install
InputFile=$1 #fastq
OutputFold=$2
OutputFile=$3
TmpFold=$4
Keyword=$5
ToolFold=$DatabaseFold"Tool/"
ViennaFold=$ToolFold"ViennaRNA-2.4.17/"
RNAfoldFold=$ViennaFold"bin/"
CurrentPath=`pwd`
cd $OutputFold
AllSampleNum=`wc -l $InputFile | cut -d" " -f1 | awk '{print $0/2}'`
cp $InputFile $TmpFold""$Keyword".Running.fa"
while true
do
$RNAfoldFold"/RNAfold" -p --jobs=30 < $TmpFold""$Keyword".Running.fa" > $TmpFold""$Keyword".Running.fold"
awk -F" " '(!($0 ~ "WARNING" || $0 ~"scaling factor" || $0 ~"free energy" || $0~ "pf_scale" || $0 ~ "unbalanced brackets" || $0=="")){print $0}' $TmpFold""$Keyword".Running.fold" > $TmpFold""$Keyword".Complete.fold"
cat $TmpFold""$Keyword".Complete.fold" >> $OutputFile
CompleteSampleNum=`awk '(NR % 6 ==1 && ($0 ~ "ENST")){print $0}' $OutputFile | wc -l | cut -d" " -f1`
if [[ $CompleteSampleNum -ne $AllSampleNum ]]
then
awk -F"\t" 'ARGIND==1 && (FNR % 6 ==1) && ($0 ~ "ENST"){Known[$1]=$1}\
            ARGIND==2 && (FNR % 2 ==1) && (!($1 in Known)){out=$1;getline;print out"\n"$0}'\
            $OutputFile\
            $InputFile > $TmpFold""$Keyword".Running.fa"
else
break
fi
done
for((i=1;i<=$AllSampleNum;i++))
do
multipleindex=`echo $i | awk '{print $0 % 30}'`
if [[ $multipleindex -eq 1 ]] || [[ $i -eq $AllSampleNum ]]
then
sleep 60
wait
else
na1=`awk -v i="$i" '(NR==(i-1)*2+1){print substr($0,2)}' $InputFile`
Transfer $na1 &
fi
done
cd $CurrentPath
}

function StructureAnnotation(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
DownloadFold=$DatabaseFold"RawData/"
DownloadReferenceFold=$DownloadFold"Reference/"
StructureFold=$DatabaseFold"Structure/"
StructureFileFold=$StructureFold"Files/"
StructureTmpFold=$StructureFold"Tmp/"
StructureFold=$StructureFold"RNAstructure/"
mkdir $StructureFileFold
mkdir $StructureTmpFold
mkdir $StructureFold
#Without Editing
KnownTableFile=$StructureTmpFold"KnownENST.Table"
RunningTableFile=$StructureTmpFold"RunningENST.Table"
if [[ ! -f $KnownTableFile ]]
then
cut -f5 $GeneAnnotationFold""$cancer".Table" | tail -n +2 > $RunningTableFile
cp $RunningTableFile $KnownTableFile 
else
awk -F"\t" 'ARGIND==1{Known[$1]=$1}\
            ARGIND==2 && (!($5 in Known)) && (FNR>1){print $5}'\
            $KnownTableFile\
            $GeneAnnotationFold""$cancer".Table" > $RunningTableFile
cat $RunningTableFile >> $KnownTableFile 
fi
awk -F"\t" 'ARGIND==1{ENST[$1]=$1}\
            ARGIND==2 && (FNR %2 ==1){split($0,a,"|");if(substr(a[1],2) in ENST){getline;if(length($0)<33679){print a[1]"\n"$0}}}'\
            $RunningTableFile\
            $DownloadReferenceFold"gencode.v22.annotation.transcripts.d1.vd1.fa" > $StructureTmpFold""$cancer".withoutEditing.fa"
seqkit sort -l $StructureTmpFold""$cancer".withoutEditing.fa" > $StructureTmpFold""$cancer".withoutEditing.sorted.fa"
mkdir $StructureFold"withoutEditing/"
RNAfoldPrediction $StructureTmpFold""$cancer".withoutEditing.sorted.fa" $StructureFold"withoutEditing/" $StructureFileFold""$cancer".withoutEditing.fold" $StructureTmpFold $cancer".withoutEditing"
cat $StructureTmpFold""$cancer".withoutEditing.fold" $StructureFileFold"WithoutEditing.fold"
#Withone Editing
KnownEditingENSTFile=$StructureTmpFold"KnownEditingENST.Table"
RunningEditingENSTFile=$StructureTmpFold"RunningEditingENST.Table"
if [[ ! -f $KnownEditingENSTFile ]]
then
cut -f1,5 $GeneAnnotationFold""$cancer".Table" | tail -n +2 > $RunningEditingENSTFile
cp $RunningEditingENSTFile $KnownEditingENSTFile 
else
awk -F"\t" 'ARGIND==1{Known[$0]=$0}\
            ARGIND==2 && (!($1"\t"$5 in Known)) && (FNR>1){print $1"\t"$5}'\
            $KnownEditingENSTFile\
            $GeneAnnotationFold""$cancer".Table" > $RunningEditingENSTFile
cat $RunningEditingENSTFile >> $KnownEditingENSTFile 
fi
awk -F"\t" 'ARGIND==1 && (FNR %2 ==1){split($0,a,"|");ENST=substr(a[1],2);strand[ENST]=a[length(a)-1];Exons[ENST]=a[length(a)-2];getline;Seq[ENST]=$0}\
            ARGIND==2{PM=strand[$2];Exs=Exons[$2];WTseq=Seq[$2];split(Exs,exons,";");num=0;flag="N";split($1,a,"_");Editing=a[length(a)-1];\
                      if(PM=="+"){for(j=1;j<=length(exons);j++){split(exons[j],EachExon,"-");\
                                      if(EachExon[2]<Editing){num=num+EachExon[2]-EachExon[1]+1}\
                                      else if(EachExon[1]<=Editing && EachExon[2]>=Editing){num=num+Editing-EachExon[1]+1;flag="Y";break}\
                                      else{flag="N";break}}}\
                      else{for(j=1;j<=length(exons);j++){split(exons[j],EachExon,"-");\
                                      if(Editing<EachExon[1]){num=num+EachExon[2]-EachExon[1]+1}\
                                      else if (EachExon[1]<=Editing && EachExon[2]>=Editing){num=num+EachExon[2]-Editing+1;flag="Y";break}\
                                      else{flag="N";break}}}\
                      if(flag=="Y"){print ">"$2"."$1"\n"substr(WTseq,1,num-1)"G"substr(WTseq,num+1)}}'\
            $DownloadReferenceFold"gencode.v22.annotation.transcripts.d1.vd1.fa"\
            $RunningEditingENSTFile > $StructureTmpFold""$cancer".withoneEditing.fa"
seqkit sort -l $StructureTmpFold""$cancer".withoneEditing.fa" > $StructureTmpFold""$cancer".withoneEditing.sorted.fa"
mkdir $StructureFold"withoneEditing/"
RNAfoldPrediction $StructureTmpFold""$cancer".withoneEditing.sorted.fa" $StructureFold"withoneEditing/" $StructureFileFold""$cancer".withoneEditing.fold" $StructureTmpFold $cancer".withoneEditing"
#WithMultiple Editing
KnownEditingENSTFile=$StructureTmpFold"KnownMultipleEditingENST.Table"
RunningEditingENSTFile=$StructureTmpFold"RunningMultipleEditingENST.Table"
cut -f2,6 ../GeneralAnnotation/GeneralInformation | tail -n +2 | sort | uniq |\
     awk -F"\t" '{if(!(ENSG[$2] ~ $1)){ENSG[$2]=ENSG[$2]";"$1}}\
              END{for(i in ENSG){split(ENSG[i],a,";");delete Ed;\
                      for(j=2;j<=length(a);j++){split(a[j],b,"_");Ed[b[2]]=b[2]};L=asort(Ed,RankEd);out="";\
                      for(j=1;j<=L;j++){out=out";"RankEd[j]};print i"\t"out}}' > $RunningEditingENSTFile
awk -F"\t" 'ARGIND==1 && (FNR %2 ==1){split($0,a,"|");ENST=substr(a[1],2);strand[ENST]=a[length(a)-1];Exons[ENST]=a[length(a)-2];getline;Seq[ENST]=$0}\
            ARGIND==2{PM=strand[$1];Exs=Exons[$1];WTseq=Seq[$1];split(Exs,exons,";");num=0;split($2,Eds,";");delete new_number;\
                      if(PM=="+"){pre=1;for(j=1;j<=length(exons);j++){split(exons[j],EachExon,"-");\
                                      for(m=pre+1;m<=length(Eds);m++){\
                                        if(EachExon[2]<Eds[m]){num=num+EachExon[2]-EachExon[1]+1;pre=m-1;break}\
                                        else if(EachExon[1]<=Eds[m] && EachExon[2]>=Eds[m]){new_number[m]=num+Eds[m]-EachExon[1]+1;Flag[m]="Y";pre=m}\
                                        else{Flag[m]="N"}\
                                      }\
                                  }\
                      pre=1;out="";\
                      for(m=2;m<=length(Eds);m++){if(Flag[m]=="Y"){out=out""substr(WTseq,pre,new_number[m]-pre)"G";pre=new_number[m]+1}}\
                      }\
                      else{pre=length(Eds)+1;for(j=1;j<=length(exons);j++){split(exons[j],EachExon,"-");\
                                for(m=pre-1;m>=2;m--){\
                                  if(EachExon[1]>Eds[m]){num=num+EachExon[2]-EachExon[1]+1;pre=m+1;break}\
                                  else if(EachExon[1]<=Eds[m] && EachExon[2]>=Eds[m]){new_number[m]=num+EachExon[2]-Eds[m]+1;Flag[m]="Y";pre=m}\
                                  else{Flag[m]="N"}\
                                }\
                            }\
                      pre=1;out="";\
                      for(m=length(Eds);m>=2;m--){if(Flag[m]=="Y"){out=out""substr(WTseq,pre,new_number[m]-pre)"G";pre=new_number[m]+1}}\
                     };if(pre!=1){print ">"$1"\n"out""substr(WTseq,pre)}}'\
            $DownloadReferenceFold"gencode.v22.annotation.transcripts.d1.vd1.fa"\
            $RunningEditingENSTFile | awk -F"\t" '(substr($1,1,1)==">"){out=$0;getline;if(length($0)<=33679){print out"\n"$0}}'> $StructureTmpFold"WithMultipleEditing.fa"
seqkit sort -l $StructureTmpFold""$cancer".withoneEditing.fa" 
}

function Fold2Table(){
WithoutInputFile=$1
WithoneInputFile=$2
WithMultipleFile=$3
TmpFold=$4
awk -F"\t" '(substr($1,1,1)==">"){header=substr($0,2);getline;getline;\
                                  split($0,a,"(");gsub(/ /,"",a[length(a)]);print header"\t"substr(a[length(a)],1,length(a[length(a)])-1)}'\
                                  $WithoutInputFile > $TmpFold"WithoutEditing.table"
awk -F"\t" 'ARGIND==1{Known[$1]=$1}\
            ARGIND==2 && (substr($1,1,1)==">") && (!(substr($1,2) in Known)){NoPS[substr($1,2)]=substr($1,2)}\
            ARGIND==3 && (substr($1,1,1)==">") && (substr($1,2) in NoPS){header=substr($0,2);getline;getline;\
                                  split($0,a,"(");gsub(/ /,"",a[length(a)]);print header"\t"substr(a[length(a)],1,length(a[length(a)])-1)}'\
            $TmpFold"WithoutEditing.table"\
            ${WithoutInputFile/".fold"/".fa"}\
            ${WithoutInputFile/".fold"/".noPS.fold"} >> $TmpFold"WithoutEditing.table"
awk -F"\t" '(substr($1,1,1)==">"){header=substr($0,2);split(header,a,".chr");ENSG=a[1];Editing="chr"a[2];getline;getline;\
                                  split($0,a,"(");gsub(/ /,"",a[length(a)]);print ENSG"\t"Editing"\t"substr(a[length(a)],1,length(a[length(a)])-1)}'\
                                  $WithoneInputFile > $TmpFold"WithoneEditing.table"
awk -F"\t" '(substr($1,1,1)==">"){header=substr($0,2);getline;getline;\
                                  split($0,a,"(");gsub(/ /,"",a[length(a)]);print header"\t"substr(a[length(a)],1,length(a[length(a)])-1)}'\
                                  $WithMultipleFile > $TmpFold"WithmultipleEditing.table"
awk -F"\t" 'ARGIND==1{Known[$1]=$1}\
            ARGIND==2 && (substr($1,1,1)==">") && (!(substr($1,2) in Known)){NoPS[substr($1,2)]=substr($1,2)}\
            ARGIND==3 && (substr($1,1,1)==">") && (substr($1,2) in NoPS){header=substr($0,2);getline;getline;\
                                  split($0,a,"(");gsub(/ /,"",a[length(a)]);print header"\t"substr(a[length(a)],1,length(a[length(a)])-1)}'\
            $TmpFold"WithmultipleEditing.table"\
            ${WithMultipleFile/".fold"/".fa"}\
            ${WithMultipleFile/".fold"/".noPS.fold"} >> $TmpFold"WithmultipleEditing.table"
echo -e "CAeditomeID\tENST\tEditingInformation\tWithoneEditingMFE\tWideTypeMFE\tWithMultipleEditingMFE" > $TmpFold"RNAstructureInformation"
awk -F"\t" 'ARGIND==1 && (FNR>1){CAeditomeID[$2]=$1}\
            ARGIND==2{Without[$1]=$2}\
            ARGIND==3{WithMultiple[$1]=$2}\
            ARGIND==4{print CAeditomeID[$2]"\t"$1"\t"$2"\t"$3"\t"Without[$1]"\t"WithMultiple[$1]}'\
            /data9/swu13/CAeditome_hg38/GeneralAnnotation/GeneralInformation\
            $TmpFold"WithoutEditing.table"\
            $TmpFold"WithmultipleEditing.table"\
            $TmpFold"WithoneEditing.table" >> $TmpFold"RNAstructureInformation"            
}

function main(){
index=$1
StructureAnnotation $index &
}