cancertypes=("BLCA"  "LGG"  "BRCA"  "LIHC"  "CESC"  "LUAD"  "COAD"  "LUSC"  "GBM"  "PRAD"  "HNSC"  "STAD"  "KICH"  "THCA"  "KIRC"  "UCEC"  "KIRP"  "ACC"  "DLBC"  "ESCA"  "LAML"  "OV"  "PAAD"  "PCPG"  "SKCM"  "UCS"  "UVM"  "CHOL"  "MESO"  "READ"  "SARC"  "TGCT"  "THYM")
DatabaseFold="/data9/swu13/CAeditome_hg38/"

function MkdirPath(){
mkdir $DatabaseFold"Splicing/"
}

function LiftOver_start_end(){
#Attention: there are some problems with liftover, So we did the filtering as below. Otherwise, we cannot extract the splicing regions.
#1. we need to keep the smaller values which mapped to smaller positions.
#2. we need to keep the exons in plus/minus strand which mapped to the same strand. we need to keep the 5'SS/3'SS which mapped to the same splicing region.
InputFile=$1 # first line: chr start end name
ChainFile=$2
Num1=$3 # number of exons
OutputFile1=$4
OutputFile2=$5
Tmpfold=$6
mkdir $Tmpfold
ToolFold=$DatabaseFold"Tool/"
awk -F"\t" '(NR>1){split($4,a,":");for(i=1;i<=length(a);i++){print "chr"$3"\t"a[i]-1"\t"a[i]"\t"$1};\
                   split($5,b,":");for(i=1;i<=length(b);i++){print "chr"$3"\t"b[i]-1"\t"b[i]"\t"$1}}' $InputFile > $Tmpfold"Input"
$ToolFold"liftover/liftOver" $Tmpfold"Input" $ChainFile $Tmpfold"Output.mapped" $Tmpfold"Output.unmapped"
sed -i '1~2d' $Tmpfold"Output.unmapped"
awk -F"\t" 'ARGIND==1{a[$4]=$4} ARGIND==2 && (!($1 in a)){print $0}' $Tmpfold"Output.unmapped" $InputFile > $OutputFile1
head -n 1 $InputFile > $OutputFile2
awk -F"\t" -v num1="$Num1" 'function find_min(a,b){if(a<=b){return a}else{return b}}\
                            function find_max(a,b){if(a<b){return b}else{return a}}\
                            ARGIND==1{Po[$4]=Po[$4]":"$3}\
                            ARGIND==2 && (FNR>1) && ($1 in Po){split(substr(Po[$1],2),a,":");exons1="";exons2="";\
                                                 for(i=1;i<=num1;i++){exons1=exons1":"a[i]};\
                                                 for(i=num1+1;i<=length(a);i++){exons2=exons2":"a[i]};min_index=num1/2;max_index=min_index+1;\
                                                 split($4,a,":");split(substr(exons1,2),b,":");
                                                 split($5,c,":");split(substr(exons2,2),d,":");flag=1;\
                                                 if((a[min_index]<a[max_index] && b[min_index]>b[max_index]) || (a[min_index]>a[max_index] && b[min_index]<b[max_index])){flag=0};
                                                 for(i=1;i<=length(c);i=i+2){if((c[i]<c[i+1] && d[i]>d[i+1]) || (c[i]>c[i+1] && d[i]<d[i+1])){flag=0}};\
                                                 if(flag==1){print $1"\t"$2"\t"$3"\t"substr(exons1,2)"\t"substr(exons2,2)"\t"$6"\t"$7}}' $Tmpfold"Output.mapped" $OutputFile1 >> $OutputFile2
}

function SplicingSeq(){
InputFile=$1
OutputFile=$2
TmpFold=$3
DownloadFold=$DatabaseFold"RawData/"
DownloadReferenceFold=$DownloadFold"Reference/"
ReferenceFile=$DownloadReferenceFold"GRCh38.d1.vd1.fa"
mkdir $Tmpfold
awk -F"\t" '(FNR>1){if($8!="."){split($8,a,":");print "chr"$3"\t"a[1]-1"\t"a[2]"\t"$1"|5SS\t.\t"$7};\
                    if($9!="."){split($9,a,":");print "chr"$3"\t"a[1]-1"\t"a[2]"\t"$1"|3SS\t.\t"$7}}' $InputFile > $TmpFold"Input.bed"
bedtools getfasta -fi $ReferenceFile -bed $TmpFold"Input.bed" -tab -s -fo $TmpFold"Input.fa" -name
awk -F"\t" 'ARGIND==1{split($1,a,":");Seq[a[1]]=$2}\
            ARGIND==2{if(FNR==1){print $0"\t5SSseq\t3SSseq"}else{out=$0;if($1"|5SS" in Seq){out=out"\t"Seq[$1"|5SS"]}else{out=out"\t."};
                                                                        if($1"|3SS" in Seq){out=out"\t"Seq[$1"|3SS"]}else{out=out"\t."};print out}}'\
            $TmpFold"Input.fa"\
            $InputFile > $OutputFile
}

function SpliceStrengthEstimate(){
#Install MaxEntSan
#ToolFold=$DatabaseFold"Tool/"
#wget -c http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz -O $ToolFold"MaxEntScan.tar.gz"
#tar -zxvf $ToolFold"MaxEntScan.tar.gz" -C $ToolFold
#mv $ToolFold"fordownload" $ToolFold"MaxEntScan"
InputFile=$1
OutputFile=$2
TmpFold=$3
col1=$4
col2=$5
mkdir $TmpFold
awk -F"\t" -v col="$col1" '(FNR>1){if(($col!=".") && ($col!="NNNNNNNNN")){print ">"FNR"\n"$col}}' $InputFile > $TmpFold"ForEstimate.5ss"
awk -F"\t" -v col="$col2" '(FNR>1){if(($col!=".") && ($col!="NNNNNNNNNNNNNNNNNNNNNNN")){print ">"FNR"\n"$col}}' $InputFile > $TmpFold"ForEstimate.3ss"
ScriptFold=`pwd`
cd $ToolFold"MaxEntScan/"
perl $ToolFold"MaxEntScan/score3.pl" $TmpFold"ForEstimate.3ss" > $TmpFold"ForEstimate.3ss.strength"
perl $ToolFold"MaxEntScan/score5.pl" $TmpFold"ForEstimate.5ss" > $TmpFold"ForEstimate.5ss.strength"
awk -F"\t" -v col1="$col1" -v col2="$col2" 'ARGIND==1{Str5[$1]=$2}\
                                            ARGIND==2{Str3[$1]=$2}\
                                            ARGIND==3{if(FNR==1){print $0"\t5SSstrength\t3SSstrength"}\
                                                      else{out=$0;if($col1!="NNNNNNNNN" && $col!="."){out=out"\t"Str5[$col1]}else{out=out"\t."};\
                                                                  if($col2!="NNNNNNNNNNNNNNNNNNNNNNN" && $co2!="."){out=out"\t"Str3[$col2]}else{out=out"\t."};print out}}'\
                                            $TmpFold"ForEstimate.5ss.strength"\
                                            $TmpFold"ForEstimate.3ss.strength"\
                                            $InputFile > $OutputFile
cd $ScriptFold
}

function RNAeditingMatchSplicing(){
RNAediting=$1
Splicing=$2
OutputFile=$3
TmpFold=$4
mkdir $TmpFold
awk -F"\t" '(NR>1){split($1,a,"_");chromo=a[1];for(i=2;i<=length(a)-2;i++){chromo=chromo"_"a[i]};\
                   print chromo"\t"a[length(a)-1]-1"\t"a[length(a)-1]"\t.\t.\t"a[length(a)]}' $RNAediting | sort | uniq > $TmpFold"RNAediting"
awk -F"\t" '(NR>1){if($8!="."){split($8,a,":");print "chr"$3"\t"a[1]-1"\t"a[2]"\t"$1"|"$2"|"$3"|"$4"|"$5"|"$6"|"$7"|"$8"|"$10"|"$12"\t.\t"$7}}' $Splicing > $TmpFold"5SS"
awk -F"\t" '(NR>1){if($9!="."){split($9,a,":");print "chr"$3"\t"a[1]-1"\t"a[2]"\t"$1"|"$2"|"$3"|"$4"|"$5"|"$6"|"$7"|"$9"|"$11"|"$13"\t.\t"$7}}' $Splicing > $TmpFold"3SS"
head -n 1 $Splicing | awk '{print "RNAediting\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\tSplicingRegion\tSplicingSequence\tWTSplicingStrength\tEditedType\tEditedPosition\tEditedSS"}' > $OutputFile
if [[ -f $TmpFold"5SS" ]]
then
bedtools intersect -a $TmpFold"RNAediting" -b $TmpFold"5SS" -s -wa -wb > $TmpFold"RNAediting.5SS"
awk -F"\t" '{out=$1"_"$3"_"$6;split($10,a,"|");for(i=1;i<=length(a);i++){out=out"\t"a[i]};out=out"\t5SS";\
             if($6=="+"){po=$3-$8}else{po=$9-$2};if(po>3){out=out"\t"po-3"i"}else{out=out"\t"4-po"e"};out=out"\t"substr(a[9],1,po-1)"G"substr(a[9],po+1);if(a[9]!="."){print out}}'\
            $TmpFold"RNAediting.5SS" >> $OutputFile
fi
if [[ -f $TmpFold"3SS" ]]
then
bedtools intersect -a $TmpFold"RNAediting" -b $TmpFold"3SS" -s -wa -wb > $TmpFold"RNAediting.3SS"
awk -F"\t" '{out=$1"_"$3"_"$6;split($10,a,"|");for(i=1;i<=length(a);i++){out=out"\t"a[i]};out=out"\t3SS";\
             if($6=="+"){po=$3-$8}else{po=$9-$2};if(po>20){out=out"\t"po-20"e"}else{out=out"\t"21-po"i"};out=out"\t"substr(a[9],1,po-1)"G"substr(a[9],po+1);if(a[9]!="."){print out}}'\
            $TmpFold"RNAediting.3SS" >> $OutputFile
fi
}

function RNAeditingSpliceStrengthEstimate(){
#Install MaxEntSan
#ToolFold=$DatabaseFold"Tool/"
#wget -c http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz -O $ToolFold"MaxEntScan.tar.gz"
#tar -zxvf $ToolFold"MaxEntScan.tar.gz" -C $ToolFold
#mv $ToolFold"fordownload" $ToolFold"MaxEntScan"
ToolFold=$DatabaseFold"Tool/"
InputFile=$1
OutputFile=$2
TmpFold=$3
col1=$4
col2=$5
mkdir $TmpFold
awk -F"\t" -v col1="$col1" -v col2="$col2" '(FNR>1){if($col1=="5SS"){print ">"FNR"\n"$col2}}' $InputFile > $TmpFold"ForEstimate.5ss"
awk -F"\t" -v col1="$col1" -v col2="$col2" '(FNR>1){if($col1=="3SS"){print ">"FNR"\n"$col2}}' $InputFile > $TmpFold"ForEstimate.3ss"
ScriptFold=`pwd`
cd $ToolFold"MaxEntScan/"
perl $ToolFold"MaxEntScan/score3.pl" $TmpFold"ForEstimate.3ss" > $TmpFold"ForEstimate.3ss.strength"
perl $ToolFold"MaxEntScan/score5.pl" $TmpFold"ForEstimate.5ss" > $TmpFold"ForEstimate.5ss.strength"
awk -F"\t" -v col1="$col1" -v col2="$col2" 'ARGIND==1{Str5[$1]=$2}\
                                            ARGIND==2{Str3[$1]=$2}\
                                            ARGIND==3{if(FNR==1){print $0"\tEditedSSstrength"}\
                                                      else{out=$0;if($col1=="5SS"){out=out"\t"Str5[$col2]};\
                                                                  if($col1=="3SS"){out=out"\t"Str3[$col2]};print out}}'\
                                            $TmpFold"ForEstimate.5ss.strength"\
                                            $TmpFold"ForEstimate.3ss.strength"\
                                            $InputFile > $OutputFile
cd $ScriptFold
}

function RNAeditingPSI(){
MatchFile=$1
RNAeditingFile=$2
SplicingFold=$3
OutputFile=$4
CorrelationFold=$5
BoxplotFold=$6
TmpFold=$7
mkdir $TmpFold
#RNA editing and PSI match file
head -n 1 $RNAeditingFile | tr '\t' '\n' | awk -F"-" '(NR>1) && (substr(a[4],1,1)==0){print $0"\t"NR}' > $TmpFold"RNAediting.sample"
SplicingFile=`ls $SplicingFold* | head -n 1`
zcat $SplicingFile | head -n 1 | tr '\t' '\n' | awk -F"-" '(NR>6) && (substr(a[4],1,1)==0){print $1"-"$2"-"$3"-"$4"\t"NR}' > $TmpFold"PSI.sample"
awk -F"\t" 'ARGIND==1{PSI[$1]=PSI[$1]";"$2}\
            ARGIND==2 && ($1 in PSI){print $0"\t"substr(PSI[$1],2)}'\
            $TmpFold"PSI.sample"\
            $TmpFold"RNAediting.sample" > $TmpFold"RNAediting.PSI.sample"
#PSI file
PSIcolumn=`cut -f3 $TmpFold"RNAediting.PSI.sample"`
rm -f $TmpFold"PSI.All"
for file in `ls $SplicingFold*`
do
zcat $file | awk -F"\t" -v PSIcolumn="${PSIcolumn[*]}" '(FNR>1){split(PSIcolumn,a," ");out=$1;\
                                                                for(i=1;i<=length(a);i++){split(a[i],b,";");m=0;\
                                                                for(j=1;j<=length(b);j++){m=m+$b[j]};out=out"\t"m/length(b)};print out}' >> $TmpFold"PSI.All"
done
awk -F"\t" 'ARGIND==1{PSI[$1]=$0}\
            ARGIND==2 && (FNR>1){if($2~"mutex_exons"){print PSI[substr($2,1,length($2)-2)]}else{print PSI[$2]}}'\
            $TmpFold"PSI.All"\
            $MatchFile > $TmpFold"PSI"
#RNA editing file
RNAeditingColumn=`cut -f2 $TmpFold"RNAediting.PSI.sample"`
awk -F"\t" -v RNAeditingColumn="${RNAeditingColumn[*]}" '\
              ARGIND==1 && (FNR>1){split(RNAeditingColumn,a," ");out=$1;for(i=1;i<=length(a);i++){out=out"\t"$a[i]};RNAediting[$1]=out}\
              ARGIND==2 && (FNR>1){print RNAediting[$1]}'\
              $RNAeditingFile $MatchFile > $TmpFold"RNAediting"
#Correlation
Rscript SplicingCorrelation.R $TmpFold"RNAediting" $TmpFold"PSI" $BoxplotFold $CorrelationFold $TmpFold"correlation"
cut -f1-2 --complement $TmpFold"correlation" > $TmpFold"correlation.1"
paste $MatchFile $TmpFold"correlation.1" > $OutputFile
rm -f $TmpFold"correlation.1"
}

function WTSplicingAnnotation(){
Keyword=$1
DownloadFold=$DatabaseFold"RawData/"
DownloadReferenceFold=$DownloadFold"Reference/"
ToolFold=$DatabaseFold"Tool/"
RawSplicingFold=$DownloadFold"Splicing/"
SplicingFold=$DatabaseFold"Splicing/"
SplicingFileFold=$SplicingFold"Files/"
SplicingTmpFold=$SplicingFold"Tmp/"
SplicingBoxplot=$SplicingFold"Boxplot/"
SplicingCorrelation=$SplicingFold"Correlation/"
mkdir $SplicingFileFold
mkdir $SplicingTmpFold
mkdir $SplicingBoxplot
mkdir $SplicingCorrelation
PSIFile=$RawSplicingFold"merge_graphs_"$Keyword"_C2.confirmed.txt.gz"
Num=`zcat $PSIFile | head -n 2 | tail -n 1 | cut -f4 | awk -F":" '{print NF}'`
zcat $PSIFile | cut -f1-6 > $SplicingTmpFold""$Keyword".v19"
awk -F"\t" 'ARGIND==1{gsub(/"| /,"",$9);split($9,a,";");split(a[1],b,".");ENSG=substr(b[1],length("gene_id")+1);Strand[ENSG]=$7}\
            ARGIND==2{split($6,a,".");if(FNR==1){print $0"\tStrand"}else{print $0"\t"Strand[a[1]]}}'\
            $DownloadReferenceFold"gencode.v19.annotation.gtf"\
            $SplicingTmpFold""$Keyword".v19" > $SplicingTmpFold""$Keyword".v19.strand"
LiftOver_start_end $SplicingTmpFold""$Keyword".v19.strand" $ToolFold"liftover/hg19ToHg38.over.chain.gz" $Num $SplicingTmpFold""$Keyword".v19.strand.matched" $SplicingTmpFold""$Keyword".v38.strand.matched" $SplicingTmpFold"Tmp/"
awk -F"\t" 'ARGIND==1{gsub(/"| /,"",$9);split($9,a,";");split(a[1],b,".");ENSG=substr(b[1],length("gene_id")+1);Strand[ENSG]=$7}\
            ARGIND==2{split($6,a,".");if(FNR==1){print $0}else{if($7==Strand[a[1]]){print $0}}}'\
            $DownloadReferenceFold"gencode.v22.annotation.gtf"\
            $SplicingTmpFold""$Keyword".v38.strand.matched" > $SplicingTmpFold""$Keyword".v38.strand.matched.select"
if [[ $Keyword == "exon_skip" ]]
then
awk -F"\t" 'function find_min(a,b){if(a<=b){return a}else{return b}}\
            function find_max(a,b){if(a<b){return b}else{return a}}\
            {split($4,a,":");min=find_min(a[3],a[4]);max=find_max(a[3],a[4]);\
             if($7=="+"){print $0"\t"max-2":"max+6"\t"min-20":"min+2}\
             else if($7=="-"){print $0"\t"min-6":"min+2"\t"max-2":"max+20}\
             else if(NR==1){print $0"\t5SS\t3SS"}}' $SplicingTmpFold""$Keyword".v38.strand.matched.select" > $SplicingTmpFold""$Keyword".v38.matched.35ss"
elif [[ $Keyword == "intron_retention" ]]
then
awk -F"\t" 'function find_min(a,b){if(a<=b){return a}else{return b}}\
            function find_max(a,b){if(a<b){return b}else{return a}}\
            {split($4,a,":");min=find_min(a[2],a[3]);max=find_max(a[2],a[3]);\
            if($7=="+"){print $0"\t"min-2":"min+6"\t"max-20":"max+2}\
            else if($7=="-"){print $0"\t"max-6":"max+2"\t"min-2":"min+20}\
            else if(NR==1){print $0"\t5SS\t3SS"}}' $SplicingTmpFold""$Keyword".v38.strand.matched.select" > $SplicingTmpFold""$Keyword".v38.matched.35ss"
elif [[ $Keyword == "alt_3prime" ]]
then
awk -F"\t" 'function find_min(a,b){if(a<=b){return a}else{return b}}\
            function find_max(a,b){if(a<b){return b}else{return a}}\
            {split($4,a,":");min=find_min(a[3],a[4]);max=find_max(a[3],a[4]);\
            if(($7=="+")){print $0"\t.\t.\t"min-20":"min+2}\
            else if(($7=="-")){print $0"\t.\t.\t"max-2":"max+20}\
            else if(NR==1){print $0"\t5SS\t3SS"}}' $SplicingTmpFold""$Keyword".v38.strand.matched.select" > $SplicingTmpFold""$Keyword".v38.matched.35ss"
elif [[ $Keyword == "alt_5prime" ]] #The annotation of strand for this type is opposite
then
awk -F"\t" 'function find_min(a,b){if(a<=b){return a}else{return b}}\
            function find_max(a,b){if(a<b){return b}else{return a}}\
            {split($4,a,":");min=find_min(a[3],a[4]);max=find_max(a[3],a[4]);\
            if($7=="+"){print $0"\t"max-2":"max+6"\t.\t."}\
            else if($7=="-"){print $0"\t"min-6":"min+2"\t.\t."}\
            else if(NR==1){print $0"\t5SS\t3SS"}}' $SplicingTmpFold""$Keyword".v38.strand.matched.select" > $SplicingTmpFold""$Keyword".v38.matched.35ss"
elif [[ $Keyword == "mutex_exons" ]]
then
awk -F"\t" 'function find_min(a,b){if(a<=b){return a}else{return b}}\
            function find_max(a,b){if(a<b){return b}else{return a}}\
            {split($4,a,":");min1=find_min(a[3],a[4]);max1=find_max(a[3],a[4]);min2=find_min(a[5],a[6]);max2=find_max(a[5],a[6]);\
            if($7=="+"){print $1"_1\t"$2"_1\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"max1-2":"max1+6"\t"min1-20":"min1+2"\n"\
                              $1"_2\t"$2"_2\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"max2-2":"max2+6"\t"min2-20":"min2+2}\
            else if($7=="-"){print $1"_1\t"$2"_1\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"min1-6":"min1+2"\t"max1-2":"max1+20"\n"\
                                   $1"_2\t"$2"_2\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"min2-6":"min2+2"\t"max2-2":"max2+20}\
            else if(NR==1){print $0"\t5SS\t3SS"}}' $SplicingTmpFold""$Keyword".v38.strand.matched" > $SplicingTmpFold""$Keyword".v38.matched.35ss"
fi
SplicingSeq $SplicingTmpFold""$Keyword".v38.matched.35ss" $SplicingTmpFold""$Keyword".v38.matched.35ss.seq" $SplicingTmpFold"Tmp/"
SpliceStrengthEstimate $SplicingTmpFold""$Keyword".v38.matched.35ss.seq" $SplicingFileFold""$Keyword".WT.strength" $SplicingTmpFold"Tmp/" 10 11
}

function WTAnnotation(){
SplicingFold=$DatabaseFold"Splicing/"
WTSplicingAnnotation "exon_skip"
WTSplicingAnnotation "intron_retention"
WTSplicingAnnotation "alt_3prime"
WTSplicingAnnotation "alt_5prime"
WTSplicingAnnotation "mutex_exons"
cp $SplicingFileFold"exon_skip.WT.strength" $SplicingFileFold"WT.strength"
tail -n +2 $SplicingFileFold"intron_retention.WT.strength" >> $SplicingFileFold"WT.strength"
tail -n +2 $SplicingFileFold"alt_3prime.WT.strength" >> $SplicingFileFold"WT.strength"
tail -n +2 $SplicingFileFold"alt_5prime.WT.strength" >> $SplicingFileFold"WT.strength"
tail -n +2 $SplicingFileFold"mutex_exons.WT.strength" >> $SplicingFileFold"WT.strength"
#rm -f $SplicingFileFold*".WT.strength"
}

function RNAeditingStrengthAnnotation(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
RawSplicingFold=$DownloadFold"Splicing/"
DownloadRNAeditingFold=$DownloadFold"RNAediting/"
RNAeditingFile=$DownloadRNAeditingFold""$cancer"_Editing.txt"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
SplicingFold=$DatabaseFold"Splicing/"
SplicingFileFold=$SplicingFold"Files/"
SplicingTmpFold=$SplicingFold"Tmp/"
#Editing list
mkdir $SplicingTmpFold""$cancer"/"
RunningTableFile=$SplicingTmpFold""$cancer"/RunningEditing.Table"
cut -f1 $GeneAnnotationFold""$cancer".Table" | tail -n +2 > $RunningTableFile
#map editing sites and splicing regions
RNAeditingMatchSplicing $RunningTableFile $SplicingFileFold"WT.strength" $SplicingTmpFold""$cancer".RNAediting.match" $SplicingTmpFold""$cancer"/"
RNAeditingSpliceStrengthEstimate $SplicingTmpFold""$cancer".RNAediting.match" $SplicingFileFold""$cancer".WT.RNAediting.strength" $SplicingTmpFold""$cancer"/" 12 14
#PSI computation
mkdir $SplicingFold"Boxplot/"$cancer"/"
mkdir $SplicingFold"Correlation/"$cancer"/"
RNAeditingPSI $SplicingFileFold""$cancer".WT.RNAediting.strength" $RNAeditingFile $RawSplicingFold $SplicingFileFold""$cancer".WT.RNAediting.PSI" $SplicingFold"Correlation/"$cancer"/" $SplicingFold"Boxplot/"$cancer"/" $SplicingTmpFold""$cancer"/"
}

function MergeSplicing(){
SplicingFold=$DatabaseFold"Splicing/"
SplicingFileFold=$SplicingFold"Files/"
rm -f $SplicingFileFold"SplicingInformation.1"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
SplicingFile=$SplicingFileFold""$cancer".WT.RNAediting.strength"
awk -v cancer="$cancer" 'ARGIND==1{Editing[$2]=$1}\
                         ARGIND==2 && (FNR>1){print Editing[$1]"\t"$0}'\
                         $GeneAnnotationFold"GeneralInformation"\
                         $SplicingFile >> $SplicingFileFold"SplicingInformation.1"
done
head -n 1 /data9/swu13/CAeditome_hg38/Splicing/Files/BLCA.WT.RNAediting.strength | awk -F"\t" '{print "CAeditomeID\t"$0}' > $SplicingFileFold"SplicingInformation"
cat $SplicingFileFold"SplicingInformation.1" | sort | uniq >> $SplicingFileFold"SplicingInformation"
rm -f $SplicingFileFold"SplicingInformation.1"
#
head -n 1 $SplicingFileFold"BLCA.WT.RNAediting.PSI" |awk -F"\t" '{print "Cancer\tCAeditomeID\t"$0}' > $SplicingFileFold"SplicingPSIInformation"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
SplicingFile=$SplicingFileFold""$cancer".WT.RNAediting.PSI"
awk -v cancer="$cancer" 'ARGIND==1{Editing[$2]=$1}\
                         ARGIND==2 && (FNR>1) && ($16<0.05 || $20<0.05){print cancer"\t"Editing[$1]"\t"$0}'\
                         $GeneAnnotationFold"GeneralInformation"\
                         $SplicingFile >> $SplicingFileFold"SplicingPSIInformation"
done
}

function SplicingHeatmap(){
SplicingFold=$DatabaseFold"Splicing/"
SplicingFileFold=$SplicingFold"Files/"
SplicingHeatmapFold=$SplicingFold"Heatmap/"
SplicingTmpFold=$SplicingFold"Tmp/"
GeneAnnotationFold=$DatabaseFold"GeneralAnnotation/"
mkdir $SplicingHeatmapFold
rm -f $SplicingTmpFold"Dif.Heatmap.data"
rm -f $SplicingTmpFold"Cor.Heatmap.data"
for((i=0;i<${#cancertypes[*]};i++))
do
cancer=${cancertypes[$i]}
SplicingFile=$SplicingFileFold""$cancer".WT.RNAediting.PSI"
if [[ -f $SplicingFile ]]
then
awk -F"\t" -v cancer="$cancer" 'ARGIND==1{Editing[$2]=$1;split($3,a,".");if(!(Gene[$1] ~ a[1])){Gene[$1]=Gene[$1]"|"a[1]}}\
                                ARGIND==2 && (FNR>1){split(Gene[$1],a,"|");split($8,b,".");for(i=2;i<=length(a);i++){if(a[i]==b[1]){Known[$2"\t"$3]=$2"\t"$3}}}\
                                ARGIND==3 && (FNR>1) && ($16!="NA") && ($16<0.05) && ($1"\t"$2 in Known){split($7,a,".");\
                                                             print cancer"\t"Editing[$1]"|"$2"\t"$1"\t"a[1]"\t"$19}'\
                                $GeneAnnotationFold"GeneralInformation"\
                                $SplicingFileFold"SplicingInformation"\
                                $SplicingFile >> $SplicingTmpFold"Dif.Heatmap.data"
awk -F"\t" -v cancer="$cancer" 'ARGIND==1{Editing[$2]=$1;split($3,a,".");if(!(Gene[$1] ~ a[1])){Gene[$1]=Gene[$1]"|"a[1]}}\
                                ARGIND==2 && (FNR>1){split(Gene[$1],a,"|");split($8,b,".");for(i=2;i<=length(a);i++){if(a[i]==b[1]){Known[$2"\t"$3]=$2"\t"$3}}}\
                                ARGIND==3 && (FNR>1) && ($20!="NA") && ($20<0.05) && ($1"\t"$2 in Known){split($7,a,".");\
                                                             print cancer"\t"Editing[$1]"|"$2"\t"$1"\t"a[1]"\t"$21}'\
                                $GeneAnnotationFold"GeneralInformation"\
                                $SplicingFileFold"SplicingInformation"\
                                $SplicingFile >> $SplicingTmpFold"Cor.Heatmap.data"
fi
done
mkdir $SplicingHeatmapFold"Dif/"
Rscript Heatmap.R $SplicingTmpFold"Dif.Heatmap.data" 1 4 2 5 0 99 $SplicingHeatmapFold"Dif/" #cancer\tGene\tEditing and splicing\tvalue
mkdir $SplicingHeatmapFold"Cor/"
Rscript Heatmap.R $SplicingTmpFold"Cor.Heatmap.data" 1 4 2 5 0 99 $SplicingHeatmapFold"Cor/" #cancer\tGene\tEditing and splicing\tvalue
}

function DifferentialPSI(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
SplicingFold=$DownloadFold"Splicing/"
ClinicalFile=$DownloadFold"Clinical/"$cancer".txt"
ProcessSplicingFold=$DatabaseFold"DifferentialSplicing/"
mkdir $ProcessSplicingFold
TmpFold=$ProcessSplicingFold"Tmp/"
mkdir $TmpFold
#RNA editing and PSI match file
TmpCancerFold=$TmpFold""$cancer"/"
mkdir $TmpCancerFold
awk -F"\t" '(NR>1){print $1}' $ClinicalFile > $TmpCancerFold""$cancer".sample"
SplicingFile=`ls $SplicingFold* | head -n 1`
zcat $SplicingFile | head -n 1 | tr '\t' '\n' | awk -F"." '(NR>6){print $1"\t"NR}' > $TmpCancerFold"PSI.all.sample"
awk -F"\t" 'ARGIND==1{Sample[$1]=$1}\
            ARGIND==2{split($1,a,"-");if(a[1]"-"a[2]"-"a[3] in Sample){print $0}}'\
            $TmpCancerFold""$cancer".sample"\
            $TmpCancerFold"PSI.all.sample" > $TmpCancerFold"PSI."$cancer".sample"
#PSI file
PSIcolumn=`cut -f2 $TmpCancerFold"PSI."$cancer".sample"`
cut -f1 $TmpCancerFold"PSI."$cancer".sample" | awk -F"\t" 'BEGIN{out="Events"}{out=out"\t"$1}END{print out}' > $TmpCancerFold"PSI."$cancer
for file in `ls $SplicingFold*`
do
zcat $file | awk -F"\t" -v PSIcolumn="${PSIcolumn[*]}" '(FNR>1){split(PSIcolumn,a," ");out=$1;\
                                                                for(i=1;i<=length(a);i++){out=out"\t"$a[i]};print out}' >> $TmpCancerFold"PSI."$cancer
done
#Clinical
awk -F"-" '{if(substr($4,1,1)==0){print "Tumor"}else{print "Normal"}}' $TmpCancerFold"PSI."$cancer".sample" > $TmpCancerFold"PSI.clinical"
#Correlation
ProcessSplicingFileFold=$ProcessSplicingFold"Files/"
mkdir $ProcessSplicingFileFold
ProcessSplicingBoxplotFold=$ProcessSplicingFold"Boxplot/"
mkdir $ProcessSplicingBoxplotFold
mkdir $ProcessSplicingBoxplotFold""$cancer"/"
Rscript DifferentialSplicing.R $TmpCancerFold"PSI."$cancer $TmpCancerFold"PSI.clinical" $ProcessSplicingFileFold""$cancer".DifferentialResults" $ProcessSplicingBoxplotFold""$cancer"/"
#rm -r $TmpCancerFold
}

function StageAssociatedPSI(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
SplicingFold=$DownloadFold"Splicing/"
ClinicalFile=$DownloadFold"Clinical/"$cancer".txt"
ProcessSplicingFold=$DatabaseFold"StageAssociationSplicing/"
mkdir $ProcessSplicingFold
TmpFold=$ProcessSplicingFold"Tmp/"
mkdir $TmpFold
#RNA editing and PSI match file
TmpCancerFold=$TmpFold""$cancer"/"
mkdir $TmpCancerFold
PathologicalStageColumn=`head -n 1 $ClinicalFile | tr '\t' '\n' | grep -n "pathologic_stage" | cut -d":" -f1`
if [[ -n $PathologicalStageColumn ]]
then
awk -F"\t" -v stageindex="$PathologicalStageColumn" '(NR>1){gsub(/stage |a|b|c|1|2|3/,"",$stageindex);if($stageindex=="x"){print $1"\tnone"}\
                                                                                                      else{print $1"\t"$stageindex}}'\
           $ClinicalFile > $TmpCancerFold""$cancer".sample"
Rscript Roman2Num.R $TmpCancerFold""$cancer".sample" 2 $TmpCancerFold""$cancer".sample.roman2num"
SplicingFile=`ls $SplicingFold* | head -n 1`
zcat $SplicingFile | head -n 1 | tr '\t' '\n' | awk -F"." '(NR>6){print $1"\t"NR}' > $TmpCancerFold"PSI.all.sample"
awk -F"\t" 'ARGIND==1{Sample[$1]=$3}\
            ARGIND==2{split($1,a,"-");if((a[1]"-"a[2]"-"a[3] in Sample) && substr(a[4],1,1)==0){print $0"\t"Sample[a[1]"-"a[2]"-"a[3]]}}'\
            $TmpCancerFold""$cancer".sample.roman2num"\
            $TmpCancerFold"PSI.all.sample" > $TmpCancerFold"PSI."$cancer".tumor.sample"
#PSI file
PSIcolumn=`cut -f2 $TmpCancerFold"PSI."$cancer".tumor.sample"`
cut -f1 $TmpCancerFold"PSI."$cancer".tumor.sample" | awk -F"\t" 'BEGIN{out="Events"}{out=out"\t"$1}END{print out}' > $TmpCancerFold"PSI."$cancer
for file in `ls $SplicingFold*`
do
zcat $file | awk -F"\t" -v PSIcolumn="${PSIcolumn[*]}" '(FNR>1){split(PSIcolumn,a," ");out=$1;\
                                                                for(i=1;i<=length(a);i++){out=out"\t"$a[i]};print out}' >> $TmpCancerFold"PSI."$cancer
done
#Correlation
#ProcessSplicingFileFold=$ProcessSplicingFold"Files/"
#mkdir $ProcessSplicingFileFold
#ProcessSplicingCorrelationFold=$ProcessSplicingFold"Pathological/"
#mkdir $ProcessSplicingCorrelationFold
#mkdir $ProcessSplicingCorrelationFold""$cancer"/"
#Rscript StageExpressionCorrelation.R $TmpCancerFold"PSI."$cancer $TmpCancerFold"PSI."$cancer".tumor.sample" 3 $ProcessSplicingCorrelationFold""$cancer"/" $ProcessSplicingFileFold"/"$cancer".PathStage.CorResults"
#rm -r $TmpCancerFold
fi
}

function main(){
Index=$1
RNAeditingStrengthAnnotation $Index 
}
