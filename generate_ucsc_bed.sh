#!/bin/bash
cancertypes=("BLCA"  "LGG"  "BRCA"  "LIHC"  "CESC"  "LUAD"  "COAD"  "LUSC"  "GBM"  "PRAD"  "HNSC"  "STAD"  "KICH"  "THCA"  "KIRC"  "UCEC"  "KIRP"  "ACC"  "DLBC"  "ESCA"  "LAML"  "OV"  "PAAD"  "PCPG"  "SKCM"  "UCS"  "UVM"  "CHOL"  "MESO"  "READ"  "SARC"  "TGCT"  "THYM")
########################################################################
#parameter1: input file: CAeditomeID\tCAeditingInformation\tGene\tGenePosition\tCancer
#parameter2: output folder
#parameter3: Tmp filder
########################################################################
Filename=$1
OutputFold=$2
TmpFold=$3
colors=("0,191,255" "148,0,211" "255,0,0" "0,0,255" "100,149,237" "0,100,0" "139,0,0" "105,105,105" "255,127,80" "165,42,42" "0,128,0" "128,0,128" "72,209,204" "220,20,60" "124,252,0" "255,215,0" "218,112,214" "218,80,60" "255,255,0" "255,165,0" "107,142,35" "255,235,205" "221,160,221" "0,0,139" "245,222,179" "128,128,128" "240,230,140" "135,206,235" "176,196,222" "64,224,208" "205,92,92" "255,192,203" "0,0,0")
#ENSG
cut -f3 $Filename | sort | uniq > $TmpFold"ENSG.temp"
ENSG_uniq_number=`wc -l $TmpFold"ENSG.temp" |cut -d" " -f1`
#Iterate
for ((i = 1; i <= $ENSG_uniq_number; i++))
do
na1=`awk -v i="$i" '(NR==i){print $0}' $TmpFold"ENSG.temp"`
outFilename=$OutputFold""$na1".bed"
region=`grep $na1 $Filename | awk '(NR==1){split($2,b,"_");split($4,a,"-");print b[1]":"a[1]-1"-"a[2]}'`
echo -e "browser position "$region > $outFilename
echo -e "browser hide all" >> $outFilename
echo -e "browser pack ensGene" >> $outFilename
echo -e "browser pack wgEncodeGencodeV22" >> $outFilename
echo -e "browser pack wgEncodeGencodeBasicV22" >> $outFilename
echo -e "browser pack wgEncodeGencodePseudoGeneV22" >> $outFilename
echo -e "browser pack wgEncodeGencodeCompV22" >> $outFilename
echo -e "#chrom chromStart chromEnd name" >> $outFilename
for (( j=0; j<${#cancertypes[*]}; j++ ))
do
cancer=${cancertypes[$j]}
color=${colors[$j]}
echo -e "track name=CAeditome_"$cancer" description=\"Specific RNA editing in "$cancer"\" color="$color" visibility=3" >> $outFilename
grep $na1 $Filename | awk -F"\t" -v cancer="$cancer" '($5==cancer){split($2,a,"_");print a[1]"\t"a[2]-1"\t"a[2]"\t"$1}' >> $outFilename
done
done