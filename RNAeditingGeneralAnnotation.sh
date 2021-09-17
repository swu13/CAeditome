#!/bin/bash
cancertypes=("BLCA"  "LGG"  "BRCA"  "LIHC"  "CESC"  "LUAD"  "COAD"  "LUSC"  "GBM"  "PRAD"  "HNSC"  "STAD"  "KICH"  "THCA"  "KIRC"  "UCEC"  "KIRP"  "ACC"  "DLBC"  "ESCA"  "LAML"  "OV"  "PAAD"  "PCPG"  "SKCM"  "UCS"  "UVM"  "CHOL"  "MESO"  "READ"  "SARC"  "TGCT"  "THYM")
DatabaseFold="/data9/swu13/CAeditome_hg38/"

function MkdirPath(){
mkdir $DatabaseFold"GeneralAnnotation/"
mkdir $DatabaseFold"GeneAnnotaion/"
mkdir $DatabaseFold"RepeatAnnotation/"
}

function FilterSNP(){
#SNP file: $SNPfold"PossibleSNP.dbSNP" $SNPfold"GenomeWideSNP_6.na35.annot.RefAlt.hg38" "Chr\tPosition\tStrand\tRSID\tRef\tAlt"
#Input: RNA editing list file,chr_position_+
#Output: RNA editing list file exluding the SNP data
InputFile=$1
OutputFile=$2
DownloadFold=$DatabaseFold"RawData/"
SNPfold=$DownloadFold"SNP/"
awk -F"\t" 'ARGIND==1{if($5=="A" && ($6=="G" || $6=="N")){KnownSNP[$1"_"$2"_+"]=$1"_"$2"_+"}else if($5=="T" && ($6=="C" || $6=="N")){KnownSNP[$1"_"$2"_-"]=$1"_"$2"_-"}}\
            ARGIND==2 && ($5=="A" && $6=="G"){KnownSNP[$1"_"$2"_"$3]=$1"_"$2"_"$3}\
            ARGIND==3 && (!($1 in KnownSNP)){print $0}'\
            $SNPfold"PossibleSNP.dbSNP"\
            $SNPfold"GenomeWideSNP_6.na35.annot.RefAlt.hg38"\
            $InputFile > $OutputFile
}

function RNAediting2GeneTable(){
RNAeditingListFile=$1
GeneTableFile=$2
TmpFold=$3
DownloadFold=$DatabaseFold"RawData/"
ReferenceFold=$DownloadFold"Reference/"
ReferenceGTFFile=$ReferenceFold"gencode.v22.annotation.gtf"
cat $RNAeditingListFile | sort | uniq | awk -F"\t" '{split($0,a,"_");chromo=a[1];for(i=2;i<=length(a)-2;i++){chromo=chromo"_"a[i]};print chromo"\t"a[length(a)-1]-1"\t"a[length(a)-1]"\t.\t.\t"a[length(a)]}' $RNAeditingListFile > $RNAeditingListFile".bed" 
bedtools intersect -a $RNAeditingListFile".bed" -b $ReferenceGTFFile -s -wa -wb  > $RNAeditingListFile".gtf"
awk -F"\t" '($9=="transcript"){gsub(/"| /,"",$15);split($15,Gene,";");print $1"_"$3"_"$6"\t"substr(Gene[1],length("gene_id")+1)"\t"substr(Gene[3],length("gene_type")+1)"\t"substr(Gene[5],length("gene_name")+1)"\t"substr(Gene[2],length("transcript_id")+1)"\t"substr(Gene[6],length("transcript_type")+1)"\t"substr(Gene[8],length("transcript_name")+1)}' $RNAeditingListFile".gtf" > $RNAeditingListFile".transcript"
awk -F"\t" 'ARGIND==1{Known[$1"\t"$2"\t"$3"\t"$4]=$1"\t"$2"\t"$3"\t"$4}\
            ARGIND==2 && ($9=="gene"){gsub(/"| /,"",$15);split($15,Gene,";");Editing=$1"_"$3"_"$6;ENSG=substr(Gene[1],length("gene_id")+1);\
                                      GeneType=substr(Gene[2],length("gene_type")+1);GeneName=substr(Gene[4],length("gene_name")+1);\
                                      if(!(Editing"\t"ENSG"\t"GeneType"\t"GeneName in Known)){print Editing"\t"ENSG"\t"GeneType"\t"GeneName"\t.\t.\t."}}'\
            $RNAeditingListFile".transcript" $RNAeditingListFile".gtf" > $RNAeditingListFile".gene"
echo -e "EditingInformation\tENSG\tGeneType\tGeneName\tENST\tTranscriptType\tTranscriptName" > $GeneTableFile
cat $RNAeditingListFile".transcript" >> $GeneTableFile
cat $RNAeditingListFile".gene" >> $GeneTableFile
}

function GeneTranscriptAnnotation(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
RNAeditingFold=$DownloadFold"RNAediting/"
RNAeditingFile=$RNAeditingFold""$cancer"_Editing.txt"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
TmpFold=$GeneAnnotationFold"Tmp/"
mkdir $TmpFold
cut -f1 $RNAeditingFile | tail -n +2 > $TmpFold""$cancer".RNAeditinglist"
FilterSNP $TmpFold""$cancer".RNAeditinglist" $TmpFold""$cancer".RNAeditinglist.FilterSNP"
RNAediting2GeneTable $TmpFold""$cancer".RNAeditinglist.FilterSNP" $GeneAnnotationFold""$cancer".Table" $TmpFold
#rm -r $TmpFold
}

function RNAediting2GeneRegion(){
AnnovarInput=$1
AnnovarOutput=$2
Keyword=$3
ToolFold=$DatabaseFold"Tool/"
ANNOVARToolFold=$ToolFold"annovar/"
DownloadFold=$DatabaseFold"RawData/"
ReferenceFold=$DownloadFold"Reference/AnnovarReference/GRCh38_v22/"
$ANNOVARToolFold"annotate_variation.pl" -geneanno -out $Keyword -build hg38 $AnnovarInput $ReferenceFold -dbtype refGencode2
echo -e "EditingInforamtion\tRegion\tSpecficRegion" > $AnnovarOutput
#splicing = exonic > ncRNA> UTR3 = UTR5 > intron > upstream = downstream > intergenic
RegionSeq=("splicing" "exonic" "ncRNA_splicing" "ncRNA_exonic" "UTR3" "UTR5" "intronic" "ncRNA_intronic" "upstream" "downstream" "intergenic")
awk -F"\t" '{if($6=="A" && $7=="G"){out=$3"_"$4"_+"}else{out=$3"_"$4"_-"};\
             if(!(a[out]~$1)){a[out]=a[out]";"$1}}\
             END{for(i in a){print i"\t"substr(a[i],2)}}' $Keyword".variant_function" |\
             awk -F"\t" -v RegionSeq="${RegionSeq[*]}" '{split(RegionSeq,b," ");for(j=1;j<=length(b);j++){\
                    if(j==1 || j==2 || j==7){if((!($2~"ncRNA")) && ($2~b[j])){print $0"\t"b[j];break}}\
                    else{if($2~b[j]){print $0"\t"b[j];break}}}}' >> $AnnovarOutput
}

function GeneRegionAnnotation(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
GeneRegionAnnotionFold=$DatabaseFold"GeneAnnotaion/"
TmpFold=$GeneRegionAnnotionFold"Tmp/"
mkdir $TmpFold
awk '(NR>1){split($1,a,"_");chromo=a[1];for(i=2;i<=length(a)-2;i++){chromo=chromo"_"a[i]};out=chromo"\t"a[length(a)-1]"\t"a[length(a)-1];if(a[length(a)]=="+"){out=out"\tA\tG"}else{out=out"\tT\tC"};print out}' $GeneAnnotationFold""$cancer".Table" | sort | uniq > $TmpFold""$cancer".RNAeditinglist.avinput"
RNAediting2GeneRegion $TmpFold""$cancer".RNAeditinglist.avinput" $GeneRegionAnnotionFold""$cancer".Table" $TmpFold""$cancer
}

function RNAediting2Repeat(){
AnnovarInput=$1
AnnovarOutput=$2
Keyword=$3
ToolFold=$DatabaseFold"Tool/"
ANNOVARToolFold=$ToolFold"annovar/"
DownloadFold=$DatabaseFold"RawData/"
ReferenceFold=$DownloadFold"Reference/AnnovarReference/GRCh38_v22/"
$ANNOVARToolFold"annotate_variation.pl" -regionanno -out $Keyword -build hg38 $AnnovarInput $ReferenceFold -dbtype rmsk
echo -e "EditingInforamtion\tRepeatFamily\tRepeatSubFamily\tRepeatName" > $AnnovarOutput
awk -F"\t" 'ARGIND==1{RepeatFamily[$11]=$12"\t"$13}\
            ARGIND==2{split(substr($2,6),m,",");if($6=="A" && $7=="G"){out=$3"_"$4"_+"}else if($6=="T" && $7=="C"){out=$3"_"$4"_-"};\
                      for(i=1;i<=length(m);i++){print out"\t"RepeatFamily[m[i]]"\t"m[i]}}'\
            $ReferenceFold"hg38_rmsk.txt" $Keyword".hg38_rmsk" >> $AnnovarOutput
}

function RepeatAnnotation(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
RepeatAnnotationFold=$DatabaseFold"RepeatAnnotation/"
TmpFold=$RepeatAnnotationFold"Tmp/"
mkdir $TmpFold
awk '(NR>1){split($1,a,"_");chromo=a[1];for(i=2;i<=length(a)-2;i++){chromo=chromo"_"a[i]};out=chromo"\t"a[length(a)-1]"\t"a[length(a)-1];if(a[length(a)]=="+"){out=out"\tA\tG"}else{out=out"\tT\tC"};print out}' $GeneAnnotationFold""$cancer".Table" | sort | uniq > $TmpFold""$cancer".RNAeditinglist.avinput"
RNAediting2Repeat $TmpFold""$cancer".RNAeditinglist.avinput" $RepeatAnnotationFold""$cancer".Table" $TmpFold""$cancer
}


function MergeAnnotation(){
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
if [[ $i -eq 0 ]]
then
cat $GeneAnnotationFold""$cancer".Table" > $GeneAnnotationFold"GeneralInformation"
else
tail -n +2 $GeneAnnotationFold""$cancer".Table" >> $GeneAnnotationFold"GeneralInformation"
fi
done
tail -n +2 $GeneAnnotationFold"GeneralInformation" | sort | uniq > $GeneAnnotationFold"GeneralInformation.1"
CHROM=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
rm -f $GeneAnnotationFold"GeneralInformation.2"
for((i=0;i<=${#CHROM[*]};i++))
do
chr=${CHROM[$i]}
awk -F"\t" -v chr="$chr" '{split($1,a,"_");if(a[1]=="chr"chr){out=$1"\t"a[length(a)-1];for(i=2;i<=NF;i++){out=out"\t"$i};print out}}' $GeneAnnotationFold"GeneralInformation.1" |\
                            sort -n -k 2 -t $'\t' | cut -f2 --complement >> $GeneAnnotationFold"GeneralInformation.2"
done
head -n 1 $GeneAnnotationFold"GeneralInformation" | awk '{print "CAeditomeID\t"$0}' > $GeneAnnotationFold"GeneralInformation.3"
awk 'BEGIN{num=0}{if(!($1 in CAeditomeID)){num=num+1;CAeditomeID[$1]="CAediting_"num};print CAeditomeID[$1]"\t"$0}' $GeneAnnotationFold"GeneralInformation.2" >> $GeneAnnotationFold"GeneralInformation.3"
mv $GeneAnnotationFold"GeneralInformation.3" $GeneAnnotationFold"GeneralInformation"
rm -f $GeneAnnotationFold"GeneralInformation."*
}

function MergeGeneAnnotation(){
GeneRegionAnnotionFold=$DatabaseFold"GeneAnnotaion/"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
TmpFold=$GeneRegionAnnotionFold"Tmp/"
rm -f  $TmpFold"AllTable"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
tail -n +2 $GeneRegionAnnotionFold""$cancer".Table" >> $TmpFold"AllTable"
done
head -n 1 $GeneAnnotationFold"GeneralInformation" | awk -F"\t" '{print $0"\tRegion\tSpecficRegion"}' > $GeneRegionAnnotionFold"RegionInformation"
awk -F"\t" 'ARGIND==1{Region[$1]=$2"\t"$3}\
            ARGIND==2 && (FNR>1){print $0"\t"Region[$2]}'\
            $TmpFold"AllTable"\
            $GeneAnnotationFold"GeneralInformation" >> $GeneRegionAnnotionFold"RegionInformation"
}

function MergeRepeatAnnotation(){
RepeatAnnotationFold=$DatabaseFold"RepeatAnnotation/"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
TmpFold=$RepeatAnnotationFold"Tmp/"
rm -f  $TmpFold"AllTable"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
tail -n +2 $RepeatAnnotationFold""$cancer".Table" >> $TmpFold"AllTable"
done
head -n 1 $GeneAnnotationFold"GeneralInformation" | awk -F"\t" '{print $0"\tRepeatFamily\tRepeatSubFamily\tRepeatName"}' > $RepeatAnnotationFold"RepeatInformation"
awk -F"\t" 'ARGIND==1{Repeat[$1]=$2"\t"$3"\t"$4}\
            ARGIND==2 && (FNR>1) {if($2 in Repeat){print $0"\t"Repeat[$2]}else{print $0"\tNoRepeat\tNoRepeat\tNoRepeat"}}'\
            $TmpFold"AllTable"\
            $GeneAnnotationFold"GeneralInformation" >> $RepeatAnnotationFold"RepeatInformation"
}

function main(){
index=$1
GeneTranscriptAnnotation $index
GeneRegionAnnotation $index
RepeatAnnotation $index
}