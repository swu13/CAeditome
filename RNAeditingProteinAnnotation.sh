#!/bin/bash
cancertypes=("BLCA"  "LGG"  "BRCA"  "LIHC"  "CESC"  "LUAD"  "COAD"  "LUSC"  "GBM"  "PRAD"  "HNSC"  "STAD"  "KICH"  "THCA"  "KIRC"  "UCEC"  "KIRP"  "ACC"  "DLBC"  "ESCA"  "LAML"  "OV"  "PAAD"  "PCPG"  "SKCM"  "UCS"  "UVM"  "CHOL"  "MESO"  "READ"  "SARC"  "TGCT"  "THYM")
DatabaseFold="/data9/swu13/CAeditome_hg38/"

function MkdirPath(){
ProteinAnnotionFold=$DatabaseFold"ProteinAnnotation/"
mkdir $ProteinAnnotionFold
}


function RNAediting2Protein(){
AnnovarInput=$1
AnnovarOutput=$2
keyword=$3
ToolFold=$DatabaseFold"Tool/"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
Reference=$DatabaseFold"RawData/Reference/"
ANNOVARToolFold=$ToolFold"annovar/"
$ANNOVARToolFold"annotate_variation.pl" -geneanno -out $keyword -build hg38 $AnnovarInput $Reference"AnnovarReference/GRCh38_v22/" -dbtype refGencode2
awk -F"\t" '{if($7=="A" && $8=="G"){print $4"_"$5"_+\t"$2"\t"$3};if($7=="T" && $8=="C"){print $4"_"$5"_-\t"$2"\t"$3}}' $keyword".exonic_variant_function" | awk -F"\t" '{split($3,a,",");for(i=1;i<length(a);i++){split(a[i],b,":");print $1"\t"$2"\t"b[1]"\t"b[2]"\t"b[3]"\t"b[4]"\t"b[5]}}' > $keyword".exonic.variant_function.all"  
echo -e "EditingInformation\tMutationType\tENSG\tGeneType\tGeneName\tENST\tTranscriptType\tTranscriptName\tENSPsequence\tENSP\tExon\tNucleotideMutaion\tProteinMutation\tUniprotID\tUniprotSequence\tUniprotMutation" > $AnnovarOutput
awk -F"\t" 'ARGIND==1{Transcript[$5]=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}\
            ARGIND==2 && (substr($1,1,1)==">"){split($1,a,"|");split(a[1],b,".");getline;ENSPsequence[substr(b[1],2)]=$0}\
            ARGIND==3{split($2,a,".");split($3,b,".");ENSPname[a[1]]=b[1];if($6!=""){Uniprot[a[1]]=Uniprot[a[1]]";"$6}else{Uniprot[a[1]]=Uniprot[a[1]]";."}}\
            ARGIND==4{split($1,a,"|");UniprotSequence[a[2]]=$2;UniprotSequence["."]="."}\
            ARGIND==5{split($4,a,".");split(substr(Uniprot[a[1]],2),b,";");\
                      for(i=1;i<=length(b);i++){print $1"\t"$2"\t"Transcript[$4]"\t"ENSPsequence[ENSPname[a[1]]]"\t"ENSPname[a[1]]"\t"$5"\t"$6"\t"$7"\t"b[i]"\t"UniprotSequence[b[i]]}}'\
            $GeneAnnotationFold""$cancer".Table" \
            $Reference"gencode.v22.annotation.translation.d1.vd1.fa"\
            $Reference"ENSG_ENST_Uniprot.txt" \
            $Reference"uniprot_sprot.txt" \
            $keyword".exonic.variant_function.all"  | \
            awk -F"\t" '{if($9==$15){print $0"\t"$13}else{print $0"\t."}}' >> $AnnovarOutput
}

function RNAediting2ProteinFunction(){
AnnovarInput=$1
LolliplotInput=$2
AnnovarOutput=$3
keyword=$4
ToolFold=$DatabaseFold"Tool/"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
Reference=$DatabaseFold"RawData/Reference/"
ANNOVARToolFold=$ToolFold"annovar/"
$ANNOVARToolFold"table_annovar.pl" $AnnovarInput -out $keyword -protocol dbnsfp41a -operation f -build hg38 $Reference"AnnovarReference/GRCh38_v22/" -nastring .
SIFTcol=`head -n 1 $keyword".hg38_multianno.txt" | tr '\t' '\n' | grep -n SIFT_pred | cut -d":" -f1`
Polyphen2_HDIVcol=`head -n 1 $keyword".hg38_multianno.txt" | tr '\t' '\n' | grep -n Polyphen2_HDIV_pred | cut -d":" -f1`
Polyphen2_HVARcol=`head -n 1 $keyword".hg38_multianno.txt" | tr '\t' '\n' | grep -n Polyphen2_HVAR_pred | cut -d":" -f1`
PROVEANcol=`head -n 1 $keyword".hg38_multianno.txt" | tr '\t' '\n' | grep -n PROVEAN_pred | cut -d":" -f1`
awk -F"\t" -v SIFTcol="$SIFTcol" -v Polyphen2_HDIVcol="$Polyphen2_HDIVcol" -v Polyphen2_HVARcol="$Polyphen2_HVARcol" -v PROVEANcol="$PROVEANcol" '(NR>1){if($4=="A" && $5=="G"){strand="+"};if($4=="T" && $5=="C"){strand="-"};print $1"_"$2"_"strand"\t"$SIFTcol"\t"$Polyphen2_HDIVcol"\t"$Polyphen2_HVARcol"\t"$PROVEANcol}' $keyword".hg38_multianno.txt" > $keyword".protein_function.all"  
echo -e "EditingInformation\tMutationType\tENST\tNucleotideMutaion\tProteinMutation\tSIFT_pred\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_pred\tPROVEAN_pred" > $AnnovarOutput
awk -F"\t" 'ARGIND==1{Function[$1]=$2"\t"$3"\t"$4"\t"$5}ARGIND==2 && (FNR>1){if($1 in Function){out=Function[$1]}else{out=".\t.\t.\t."};print $1"\t"$2"\t"$6"\t"$12"\t"$13"\t"out}' $keyword".protein_function.all"  $LolliplotInput >> $AnnovarOutput
}

function ProteinAnnotation(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
ProteinAnnotionFold=$DatabaseFold"ProteinAnnotation/"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
TmpFold=$ProteinAnnotionFold"Tmp/"
mkdir $TmpFold
awk '(NR>1){split($1,a,"_");chromo=a[1];for(i=2;i<=length(a)-2;i++){chromo=chromo"_"a[i]};out=chromo"\t"a[length(a)-1]"\t"a[length(a)-1];if(a[length(a)]=="+"){out=out"\tA\tG"}else{out=out"\tT\tC"};print out}' $GeneAnnotationFold""$cancer".Table" > $TmpFold""$cancer".RNAeditinglist.avinput"
RNAediting2Protein $TmpFold""$cancer".RNAeditinglist.avinput" $ProteinAnnotionFold""$cancer".lolliplot" $TmpFold""$cancer
RNAediting2ProteinFunction $TmpFold""$cancer".RNAeditinglist.avinput" $ProteinAnnotionFold""$cancer".lolliplot" $ProteinAnnotionFold""$cancer".function" $TmpFold""$cancer
}

function MergeProteinAnnotation(){
ProteinAnnotionFold=$DatabaseFold"ProteinAnnotation/"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
TmpFold=$ProteinAnnotionFold"Tmp/"
#Lolliplot
rm -f  $TmpFold"Alllolliplot"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
tail -n +2 $ProteinAnnotionFold""$cancer".lolliplot" >> $TmpFold"Alllolliplot"
done
cat $TmpFold"Alllolliplot"  |sort | uniq > $TmpFold"Alllolliplot.1"
mv $TmpFold"Alllolliplot.1" $TmpFold"Alllolliplot"
awk -F"\t" '(NR==1){print "CAeditomeID\t"$0}' $ProteinAnnotionFold"BLCA.lolliplot"  > $ProteinAnnotionFold"ProteinLolliplot"
awk -F"\t" 'ARGIND==1 && (FNR>1){Editing[$2]=$1}\
            ARGIND==2{print Editing[$1]"\t"$0}'\
            $GeneAnnotationFold"GeneralInformation" $TmpFold"Alllolliplot" >> $ProteinAnnotionFold"ProteinLolliplot"
#Function
rm -f  $TmpFold"Allfunction"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
tail -n +2 $ProteinAnnotionFold""$cancer".function" >> $TmpFold"Allfunction"
done
cat $TmpFold"Allfunction"  |sort | uniq > $TmpFold"Allfunction.1"
mv $TmpFold"Allfunction.1" $TmpFold"Allfunction"
awk -F"\t" '(NR==1){out="CAeditomeID\t"$1"\t"$2"\tENSG\tGeneName";for(i=3;i<=NF;i++){out=out"\t"$i};print out}' $ProteinAnnotionFold"BLCA.function"  > $ProteinAnnotionFold"ProteinFunction"
awk -F"\t" 'ARGIND==1 && (FNR>1){Editing[$2]=$1;Gene[$6]=$3"\t"$5}\
            ARGIND==2{out=Editing[$1]"\t"$1"\t"$2"\t"Gene[$3];for(i=3;i<=NF;i++){out=out"\t"$i};print out}'\
            $GeneAnnotationFold"GeneralInformation" $TmpFold"Allfunction" >> $ProteinAnnotionFold"ProteinFunction"
}

function main(){
index=$1
ProteinAnnotation $index
}