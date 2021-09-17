#!/bin/bash

################################################
#update the folder
ReferenceFolder="/data9/swu13/CAeditome_hg38/RawData/ReferenceEditingDetection/"
BamFolder1="/data3/TCGA-RNA-BAM-files/"
BamFolder2="/data/myang9/TCGA-RNA-bam-files/"
TmpFolder="/data9/swu13/CAeditome_hg38/TCGARNAediting/Tmp/"
RNAeditingFolder="/data9/swu13/CAeditome_hg38/TCGARNAediting/"
###############################################


function RNAeditingDetection(){
InFile=$1
OutputFold=$2
python /data9/zfan2/RNAedit/REDItools-1.2.1/build/scripts-2.7/REDItoolKnown.py -i $InFile -o $OutputFold -f $ReferenceFolder"GRCh38.d1.vd1.fa" -l $ReferenceFolder"REDIportal_hg38_sorted.txt.gz" -t 10 -m 255 -G $ReferenceFolder"gencode.v22.annotation.sorted.gtf.gz" -u
}

#function RNAeditingDetection_rerun(){
#InFile=$1
#OutputFold=$2
#python /data9/zfan2/RNAedit/REDItools-1.2.1/build/scripts-2.7/REDItoolKnown.py -i $InFile -o $OutputFold -f $ReferenceFolder"GRCh38.d1.vd1.fa" -l $ReferenceFolder"REDIportal_hg38_sorted.txt.gz" -t 10 -m 255 -G $ReferenceFolder"gencode.v22.annotation.sorted.gtf.gz" -u
#}

function MergeEditing(){
InputCancer=$1
rm -f $RNAeditingFolder""$InputCancer"/matrix.Editing"
SubBamFolder1=$BamFolder1""$InputCancer
SubBamFolder2=$BamFolder2""$InputCancer
if [[ -d $SubBamFolder1 ]]
then
SubBamFolder=$SubBamFolder1
else
SubBamFolder=$SubBamFolder2
fi
for fold in `ls $RNAeditingFolder""$InputCancer`
do
if [[ -d $SubBamFolder"/"$fold ]] 
then
BamFile=`ls $SubBamFolder"/"$fold"/"*"gdc_realn_rehead.bam" | awk -F"/" '{print $NF}'`
File=`grep $BamFile $SubBamFolder"/"*"_gdc_sample_sheet."*".tsv" | cut -f7`
awk -F"\t" -v Sample="$File" '(NR>1) && ($4!=2) && ($8~"AG"){gsub(/\[/,"",$7);gsub(/]/,"",$7);gsub(/ /,"",$7);split($7,a,",");if($4==1){print Sample"\t"$1"_"$2"_+\t"a[3]/(a[1]+a[2]+a[3]+a[4])}else{print Sample"\t"$1"_"$2"_-\t"a[3]/(a[1]+a[2]+a[3]+a[4])}}' $RNAeditingFolder""$InputCancer"/"$fold"/known_"*"/outTable_"* >> $RNAeditingFolder""$InputCancer"/matrix.Editing"
fi
done
Rscript MergeEditing.R $RNAeditingFolder""$InputCancer"/matrix.Editing" $RNAeditingFolder""$InputCancer"/"$InputCancer"_Editing.txt"
}