#!/bin/bash
cancertypes=("BLCA"  "LGG"  "BRCA"  "LIHC"  "CESC"  "LUAD"  "COAD"  "LUSC"  "GBM"  "PRAD"  "HNSC"  "STAD"  "KICH"  "THCA"  "KIRC"  "UCEC"  "KIRP"  "ACC"  "DLBC"  "ESCA"  "LAML"  "OV"  "PAAD"  "PCPG"  "SKCM"  "UCS"  "UVM"  "CHOL"  "MESO"  "READ"  "SARC"  "TGCT"  "THYM")
DatabaseFold="/data9/swu13/CAeditome_hg38/"
RNAdetectionFold="/data9/swu13/CAeditome_hg38/TCGARNAediting/"
#conda create --name CAeditome python=2.7

function MkdirPath(){
DownloadFold=$DatabaseFold"RawData/"
ToolFold=$DatabaseFold"Tool/"
mkdir $DownloadFold
mkdir $ToolFold
}

function LiftOver_pos(){
#Install
#mkdir $ToolFold"liftover" 
#wget -P $ToolFold"liftover" http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
#chmod +x $ToolFold"liftover/liftOver"
#wget -P $ToolFold"liftover" http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
InputFile=$1 # first line: chr; second line: position
ChainFile=$2
OutputFile1=$3
OutputFile2=$4
Tmpfold=$5
mkdir $Tmpfold
awk '{print $1"\t"$2-1"\t"$2}' $InputFile > $Tmpfold"Input"
$ToolFold"liftover/liftOver" $Tmpfold"Input" $ChainFile $Tmpfold"Output.mapped" $Tmpfold"Output.unmapped"
sed -i '1~2d' $Tmpfold"Output.unmapped"
awk -F"\t" 'ARGIND==1{a[$1"\t"$3]=$1"\t"$3} ARGIND==2 && (!($1"\t"$2 in a)){print $0}' $Tmpfold"Output.unmapped" $InputFile > $OutputFile1
awk 'ARGIND==1{SNP[FNR]=$1"\t"$3}ARGIND==2{out=SNP[FNR];for(i=3;i<=NF;i++){out=out"\t"$i};print out}' $Tmpfold"Output.mapped" $OutputFile1 > $OutputFile2
rm -r $Tmpfold
}

function DownloadSNPreference(){
#Install vcftools: conda install -c bioconda vcftools
#Install bedtools: conda install -c bioconda bedtools
DownloadFold=$DatabaseFold"RawData/"
ToolFold=$DatabaseFold"Tool/"
SNPfold=$DownloadFold"SNP/"
Reference=$DownloadFold"Reference/"
mkdir $SNPfold
#SNP array: "Chr\tPosition\tStrand\tRSID\tRef\tAlt"
wget -P $SNPfold http://www.affymetrix.com/Auth/analysis/downloads/na35/genotyping/GenomeWideSNP_6.na35.annot.csv.zip
unzip $SNPfold"GenomeWideSNP_6.na35.annot.csv.zip" -d $SNPfold
awk -F"," '(substr($1,1,1)!="#"){gsub(/"/,"",$0);if($3!="---" && $4!="---" && ($5=="+" || $5=="-"))print "chr"$3"\t"$4"\t"$5"\t"$2"\t"$9"\t"$10}' $SNPfold"GenomeWideSNP_6.na35.annot.csv" | tail -n +2 > $SNPfold"GenomeWideSNP_6.na35.annot.AlleleAB"
awk -F"\t" '{if($1=="chrMT"){chr="chrM"}else{chr=$1};\
             print chr"\t"$2-1"\t"$2"\t.\t.\t"$3}' $SNPfold"GenomeWideSNP_6.na35.annot.AlleleAB" > $SNPfold"GenomeWideSNP_6.na35.annot.AlleleAB.bed"
bedtools getfasta -fi $Reference"GRCh37.p13.genome.fa" -bed $SNPfold"GenomeWideSNP_6.na35.annot.AlleleAB.bed" -tab -s -fo $SNPfold"GenomeWideSNP_6.na35.annot.AlleleAB.bed.fa"
awk -F"\t" 'ARGIND==1{split($1,chr,":");split(chr[2],str,"(");split(str[1],pos,"-");SNP=chr[1]"\t"pos[2]"\t"substr(str[2],1,1);Ref[SNP]=$2}\
            ARGIND==2 {if($1=="chrMT"){Editing="chrM\t"$2"\t"$3}else{Editing=$1"\t"$2"\t"$3};\
                       if(Ref[Editing]==$5){Alt=$6}else if(Ref[Editing]==$6){Alt=$5};print $1"\t"$2"\t"$3"\t"$4"\t"Ref[Editing]"\t"Alt}'\
            $SNPfold"GenomeWideSNP_6.na35.annot.AlleleAB.bed.fa"\
            $SNPfold"GenomeWideSNP_6.na35.annot.AlleleAB" > $SNPfold"GenomeWideSNP_6.na35.annot.RefAlt"
LiftOver_pos $SNPfold"GenomeWideSNP_6.na35.annot.RefAlt" $ToolFold"liftover/hg19ToHg38.over.chain.gz" $SNPfold"GenomeWideSNP_6.na35.annot.RefAlt.hg19" $SNPfold"GenomeWideSNP_6.na35.annot.RefAlt.hg38" $SNPfold"Tmp/"
#dbSNP
wget -P $SNPfold https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
vcftools --gzvcf $SNPfold"All_20180418.vcf.gz" --out $SNPfold"All_20180418.SNPonly" --remove-indels --recode --recode-INFO-all
awk -F"\t" '(substr($1,1,1)!="#"){split($5,alt,",");for(i=1;i<=length(alt);i++){print "chr"$1"\t"$2"\tdbSNP\t"$3"\t"$4"\t"alt[i]}}' $SNPfold"All_20180418.SNPonly.recode.vcf" > $SNPfold"PossibleSNP.dbSNP"
}

function DownloadClinical(){
#Install firebrowse api: pip install firebrowse
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
DownloadClinicalFold=$DownloadFold"Clinical/"
mkdir $DownloadClinicalFold
fbget clinical cohort=$cancer page_size=10000 --outfile=$DownloadClinicalFold""$cancer".txt" 
}

function DownloadGeneExpression(){
#Install R packages: TCGAbiolinks, plyr, limma, biomaRt, SummarizedExperiment
#conda install -c conda-forge r-base #version>4.0
#BiocManager::install("TCGAbiolinks");install.packages("plyr");BiocManager::install("limma");BiocManager::install("biomaRt");BiocManager::install("SummarizedExperiment")
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
DownloadExpressionFold=$DownloadFold"GeneExpression/"
mkdir $DownloadExpressionFold
Rscript Download_Gene_expression_Data.R $cancer $DownloadExpressionFold "Gene"
}

function DownloadTranscriptExpression(){
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
DownloadIsorormFold=$DownloadFold"TranscriptExpression/"
mkdir $DownloadIsorormFold
wget "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/"$cancer"/20160128/gdac.broadinstitute.org_"$cancer".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.Level_3.2016012800.0.0.tar.gz" -O $DownloadIsorormFold""$cancer".tar.gz"
tar -zxvf $DownloadIsorormFold""$cancer".tar.gz" -C $DownloadIsorormFold
mv $DownloadIsorormFold"gdac.broadinstitute.org_"$cancer".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.Level_3.2016012800.0.0/"$cancer".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt" $DownloadIsorormFold"PCPG_expression.txt"
rm -r $DownloadIsorormFold"gdac.broadinstitute.org_"$cancer".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.Level_3.2016012800.0.0/"
#conversion
wget -O $DownloadIsorormFold"Ensemble.ucsc" 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "ucsc" /></Dataset></Query>' 
}

function DownloadSplicing(){
DownloadFold=$DatabaseFold"RawData/"
DownloadSplicingFold=$DownloadFold"Splicing/"
mkdir $DownloadSplicingFold
wget --no-check-certificate https://api.gdc.cancer.gov/data/3094de2b-1da1-40eb-96ca-27b432040ebf -O $DownloadSplicingFold"merge_graphs_alt_3prime_C2.confirmed.txt.gz"
wget --no-check-certificate https://api.gdc.cancer.gov/data/70237bc6-0144-4cb6-a328-f9a41113e91b -O $DownloadSplicingFold"merge_graphs_alt_5prime_C2.confirmed.txt.gz"
wget --no-check-certificate https://api.gdc.cancer.gov/data/281fe94f-2e49-4687-8087-552f3ab09776 -O $DownloadSplicingFold"merge_graphs_exon_skip_C2.confirmed.txt.gz"
wget --no-check-certificate https://api.gdc.cancer.gov/data/6e5dc1db-58c3-4d4f-aaf3-aa49ff25b6a1 -O $DownloadSplicingFold"merge_graphs_intron_retention_C2.confirmed.txt.gz"
wget --no-check-certificate https://api.gdc.cancer.gov/data/108c6e98-b851-48db-b708-258521617c37 -O $DownloadSplicingFold"merge_graphs_mutex_exons_C2.confirmed.txt.gz"
}

function TransferRNAediting(){
#Install
#install.packages("reshape2")
#BLCA, COAD and GBM has 3, 9, 1 repeated samples for RNA editing sites, we use the mean value of the repeated samples
#BRCA and OV each has one sample with problem of samtools index: /data3/TCGA-RNA-BAM-files/BRCA/6d135cc4-4eca-48be-8f87-bf94daf83e7a/108dc6d1-5612-4276-a682-0ba4d324cb00_gdc_realn_rehead.bam;/data3/TCGA-RNA-BAM-files/OV/397df160-9ae2-4311-a58f-f2e7327d9856/bcdf0550-9e6f-4500-9473-706df7beeae4_gdc_realn_rehead.bam
CancerIndex=$1
cancer=${cancertypes[$CancerIndex]}
DownloadFold=$DatabaseFold"RawData/"
DownloadRNAeditingFold=$DownloadFold"RNAediting/"
mkdir $DownloadRNAeditingFold
rm -f $RNAdetectionFold""$cancer"/"*Editing*
source RNAdetection.sh
MergeEditing $cancer
cp $RNAdetectionFold""$cancer"/"$cancer"_Editing.txt" $DownloadRNAeditingFold
}

function DownloadGenomeReference(){
#Install seqkit
#conda install -c bioconda seqkit
#
DownloadFold=$DatabaseFold"RawData/"
DownloadReferenceFold=$DownloadFold"Reference/"
mkdir $DownloadReferenceFold
#GTF
wget -P $DownloadReferenceFold ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz
gunzip $DownloadReferenceFold"gencode.v22.annotation.gtf.gz"
wget -P $DownloadReferenceFold ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip $DownloadReferenceFold"gencode.v19.annotation.gtf.gz"
wget -P $DownloadReferenceFold ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
gunzip $DownloadReferenceFold"gencode.v37.annotation.gtf.gz"
#GFF
wget -P $DownloadReferenceFold ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gff3.gz
gunzip $DownloadReferenceFold"gencode.v22.annotation.gff3.gz"
wget -P $DownloadReferenceFold ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gff3.gz
gunzip $DownloadReferenceFold"gencode.v37.annotation.gff3.gz"
#Fasta
wget -P $DownloadReferenceFold ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/GRCh38.p2.genome.fa.gz
gunzip $DownloadReferenceFold"GRCh38.p2.genome.fa.gz"
wget -P $DownloadReferenceFold ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz
gunzip $DownloadReferenceFold"GRCh37.p13.genome.fa.gz"
wget --no-check-certificate https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834 -O $DownloadReferenceFold"GRCh38.d1.vd1.fa.tar.gz"
tar -zxvf $DownloadReferenceFold"GRCh38.d1.vd1.fa.tar.gz" -C $DownloadReferenceFold
#Protein coding transcript sequences from gencode
wget -P $DownloadReferenceFold ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.pc_transcripts.fa.gz
gunzip $DownloadReferenceFold"gencode.v22.pc_transcripts.fa.gz"
#Protein sequences from gencode
wget -P $DownloadReferenceFold ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.pc_translations.fa.gz
gunzip $DownloadReferenceFold"gencode.v22.pc_translations.fa.gz"
#Exon transcripts sequences
awk -F"\t" '($3=="exon"){gsub(/"| /,"",$9);split($9,Gene,";");transID=substr(Gene[2],length("transcript_id")+1);Transcript[transID]=Transcript[transID]";"$4"-"$5;Strand[transID]=$7;Chr[transID]=$1}END{for(i in Transcript){print i"\t"substr(Transcript[i],2)"\t"Strand[i]"\t"Chr[i]}}' $DownloadReferenceFold"gencode.v22.annotation.gtf" > $DownloadReferenceFold"gencode.v22.annotation.exon"
awk -F"\t" '{split($2,a,";");for(i=1;i<=length(a);i++){split(a[i],b,"-");print $4"\t"b[1]-1"\t"b[2]"\t"$1"|"a[i]"|"$3"|"$4"\t.\t"$3}}' $DownloadReferenceFold"gencode.v22.annotation.exon" > $DownloadReferenceFold"gencode.v22.annotation.exon.bed"
#bedtools getfasta -fi $DownloadReferenceFold"GRCh38.p2.genome.fa" -bed $DownloadReferenceFold"gencode.v22.annotation.exon.bed" -tab -s -name -fo $DownloadReferenceFold"gencode.v22.annotation.exon.fa" 
#awk -F"\t" '{split($1,a,"|");Seq[a[1]]=Seq[a[1]]""$2;Exon[a[1]]=Exon[a[1]]";"a[2];Strand[a[1]]=a[3];split(a[4],b,":");Chr[a[1]]=b[1]}END{for(i in Exon){print ">"i"|"substr(Exon[i],2)"|"Strand[i]"|"Chr[i]"\n"Seq[i]}}' $DownloadReferenceFold"gencode.v22.annotation.exon.fa" > $DownloadReferenceFold"gencode.v22.annotation.transcripts.fa"
bedtools getfasta -fi $DownloadReferenceFold"GRCh38.d1.vd1.fa" -bed $DownloadReferenceFold"gencode.v22.annotation.exon.bed" -tab -s -name -fo $DownloadReferenceFold"gencode.v22.annotation.exon.d1.vd1.fa" 
awk -F"\t" '{split($1,a,"|");Seq[a[1]]=Seq[a[1]]""$2;Exon[a[1]]=Exon[a[1]]";"a[2];Strand[a[1]]=a[3];split(a[4],b,":");Chr[a[1]]=b[1]}END{for(i in Exon){print ">"i"|"substr(Exon[i],2)"|"Strand[i]"|"Chr[i]"\n"Seq[i]}}' $DownloadReferenceFold"gencode.v22.annotation.exon.d1.vd1.fa"  > $DownloadReferenceFold"gencode.v22.annotation.transcripts.d1.vd1.fa" #ENST/EXONS/STRAND/CHR
#protein sequences
awk -F"\t" 'ARGIND==1 && (FNR %2==1){split($1,a,"|");ENST=substr(a[1],2);getline;V1D1seq[ENST]=$0}\
            ARGIND==2 && (FNR %2==1){split($1,a,"|");ENST=substr(a[1],2);getline;if($0==V1D1seq[ENST]){Confi[ENST]=ENST}}\
            ARGIND==3 && (FNR %2==1){split($1,a,"|");ENST=a[2];header=$0;getline;if(ENST in Confi){print header"\n"$0}}'\
            $DownloadReferenceFold"gencode.v22.annotation.transcripts.d1.vd1.fa"\
            $DownloadReferenceFold"gencode.v22.pc_transcripts.fa"\
            $DownloadReferenceFold"gencode.v22.pc_translations.fa" >  $DownloadReferenceFold"gencode.v22.annotation.translation.d1.vd1.fa"
#UTR3 sequences
Rscript ExtractUTR.R $DownloadReferenceFold"gencode.v22.annotation.gtf" $DownloadReferenceFold"gencode.v22.UTR3"
Rscript ExtractUTR.R $DownloadReferenceFold"gencode.v37.annotation.gtf" $DownloadReferenceFold"gencode.v37.UTR3"
awk -F"\t" 'ARGIND==1{split($16,a,".");UTR37[a[1]]=UTR37[a[1]]";"$1"_"$2"_"$3"_"$5}\
            ARGIND==2{split($16,a,".");UTR22[a[1]]=UTR22[a[1]]";"$1"_"$2"_"$3"_"$5}\
            ARGIND==3 && ($3=="three_prime_UTR"){split($9,a,";");split(a[1],b,".");ENST=substr(b[1],length("ID=UTR3:")+1);ExonNum=substr(a[9],length("exon_number=")+1);\
                      if(UTR37[ENST]==UTR22[ENST]){print ENST"\t"ExonNum"\t"$1"\t"$4"\t"$5"\t"$7}}'\
            $DownloadReferenceFold"gencode.v37.UTR3"\
            $DownloadReferenceFold"gencode.v22.UTR3"\
            $DownloadReferenceFold"gencode.v37.annotation.gff3" > $DownloadReferenceFold"gencode.v22.UTR3.Really"
awk -F"\t" '{print $3"\t"$4-1"\t"$5"\t"$1"_"$2"\t.\t"$6}' $DownloadReferenceFold"gencode.v22.UTR3.Really" > $DownloadReferenceFold"gencode.v22.UTR3.Really.bed"
bedtools getfasta -fi $DownloadReferenceFold"GRCh38.d1.vd1.fa" -bed $DownloadReferenceFold"gencode.v22.UTR3.Really.bed" -s -name -fo $DownloadReferenceFold"gencode.v22.UTR3.Really.fa"
#lncRNA sequences
awk -F"\t" 'ARGIND==1 && ($3=="transcript"){gsub(/"| /,"",$9);split($9,a,";");ENSTtype=substr(a[5],length("transcript_type")+1);\
                                            split(a[2],b,".");ENST=substr(b[1],length("transcript_id")+1);\
                                            if(ENSTtype=="lncRNA"){LncRNA[ENST]=$1"\t"$4"\t"$5"\t"$7}}\
            ARGIND==2 && ($3=="transcript"){gsub(/"| /,"",$9);split($9,a,";");split(a[2],b,".");ENST=substr(b[1],length("transcript_id")+1);\
                                            if((ENST in LncRNA) && ($1"\t"$4"\t"$5"\t"$7==LncRNA[ENST])){print ENST"\t"$1"\t"$4"\t"$5"\t"$7}}'\
            $DownloadReferenceFold"gencode.v37.annotation.gtf"\
            $DownloadReferenceFold"gencode.v22.annotation.gtf" > $DownloadReferenceFold"gencode.v22.lncRNA.Really"
awk -F"\t" 'ARGIND==1{ENST[$1]=$1}\
            ARGIND==2 && (FNR % 2==1){split($0,a,"|");split(a[1],b,".");transcript=substr(b[1],2);if(transcript in ENST){out=$0;getline;print out"\n"$0}}'\
            $DownloadReferenceFold"gencode.v22.lncRNA.Really"\
            $DownloadReferenceFold"gencode.v22.annotation.transcripts.d1.vd1.fa" > $DownloadReferenceFold"gencode.v22.lncRNA.Really.fa"
#Protein sequences from uniprotKB
wget -P $DownloadReferenceFold ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip $DownloadReferenceFold"uniprot_sprot.fasta.gz"
seqkit fx2tab $DownloadReferenceFold"uniprot_sprot.fasta" > $DownloadReferenceFold"uniprot_sprot.txt"
#match files: Human genes (GRCh38.p13), Gene stable ID, Gene stable ID version, Transcript stable ID, Transcript stable ID version, UniProtKB Gene Name symbol, UniProtKB/Swiss-Prot ID
#UniProtKB Gene Name IDENSG-ENST-ENSP-UniprotKB  Biomart doenload http://asia.ensembl.org/biomart/martview/f5e08a3f26e2cc07d4bd52b26c975e92, martquery_0309041721_147.txt.gz
wget http://asia.ensembl.org/biomart/martresults/94?file=martquery_0312025852_27.txt.gz -O $DownloadReferenceFold"martquery_0312025852_27.txt.gz"
gunzip $DownloadReferenceFold"martquery_0312025852_27.txt.gz" 
mv $DownloadReferenceFold"martquery_0312025852_27.txt" $DownloadReferenceFold"ENSG_ENST_Uniprot.txt"
}


function DownloadANNOVARReference(){
#Install
#conda install -c bioconda ucsc-gtftogenepred
#DownloadFold=$DatabaseFold"RawData/"
#ToolFold=$DatabaseFold"Tool/"
#wget -P $ToolFold http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
#tar -zxvf $ToolFold"annovar.latest.tar.gz" -C $ToolFold
#
DownloadFold=$DatabaseFold"RawData/"
ToolFold=$DatabaseFold"Tool/"
DownloadReferenceFold=$DownloadFold"Reference/"
mkdir $DownloadReferenceFold"AnnovarReference"
DownloadAnnovarReferenceFold=$DownloadReferenceFold"AnnovarReference/GRCh38_v22/"
mkdir $DownloadAnnovarReferenceFold
ANNOVARToolFold=$ToolFold"annovar/"
#ln -s /home/swu13/anaconda2/envs/py3/lib/libssl.so.1.1 /home/swu13/anaconda2/envs/py3/lib/libssl.so.1.0.0
#ln -s /home/swu13/anaconda2/envs/py3/lib/libcrypto.so.1.1 /home/swu13/anaconda2/envs/py3/lib/libcrypto.so.1.0.0
gtfToGenePred $DownloadReferenceFold"gencode.v22.annotation.gtf" $DownloadAnnovarReferenceFold"hg38_refGencode2.txt"
#$ANNOVARToolFold"retrieve_seq_from_fasta.pl" $DownloadAnnovarReferenceFold"hg38_refGencode2.txt" -seqfile $DownloadReferenceFold"GRCh38.p2.genome.fa" -format ensGene -outfile $DownloadAnnovarReferenceFold"hg38_refGencode2Mrna.fa"
$ANNOVARToolFold"retrieve_seq_from_fasta.pl" $DownloadAnnovarReferenceFold"hg38_refGencode2.txt" -seqfile $DownloadReferenceFold"GRCh38.d1.vd1.fa" -format ensGene -outfile $DownloadAnnovarReferenceFold"hg38_refGencode2Mrna.fa"
nl $DownloadAnnovarReferenceFold"hg38_refGencode2.txt" > $DownloadAnnovarReferenceFold"hg38_refGencode21.txt"
mv $DownloadAnnovarReferenceFold"hg38_refGencode21.txt" $DownloadAnnovarReferenceFold"hg38_refGencode2.txt"
$ANNOVARToolFold"annotate_variation.pl" -buildver hg38 -downdb rmsk $DownloadAnnovarReferenceFold
$ANNOVARToolFold"annotate_variation.pl" -buildver hg38 -downdb -webfrom annovar dbnsfp41a $DownloadAnnovarReferenceFold
}

function DownloadmiRNAdata(){
DownloadFold=$DatabaseFold"RawData/"
DownloadmiRNAFold=$DownloadFold"Reference/"
wget -P $DownloadmiRNAFold ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3
awk '($3=="miRNA"){split($9,a,";");split(a[3],b,"=");print $1"\t"$4-1"\t"$5"\t"b[2]"\t.\t"$7}' $DownloadmiRNAFold"hsa.gff3" > $DownloadmiRNAFold"hsa.miRNA.bed"
bedtools getfasta -fi $DownloadmiRNAFold"GRCh38.d1.vd1.fa" -bed $DownloadmiRNAFold"hsa.miRNA.bed" -tab -s -name -fo $DownloadmiRNAFold"hsa.miRNA.fa"
awk '{split($1,a,":");gsub(/T/,"U",$2);print a[1]"\t"substr($2,2,7)"\t9606"}' $DownloadmiRNAFold"hsa.miRNA.fa" > $DownloadmiRNAFold"hsa_miRNA_miRBase.txt"
awk '{gsub(/T/,"U",$2);print ">"$1"\n"$2}' $DownloadmiRNAFold"hsa.miRNA.fa" > $DownloadmiRNAFold"hsa_miRNA_miRBase.fa"
}

function DownoadTumorRelatedGene(){
DownloadFold=$DatabaseFold"RawData/"
DownloadTumorGeneFold=$DownloadFold"TumorGene/"
mkdir $DownloadTumorGeneFold
wget -P $DownloadTumorGeneFold https://oncovar.org/resource/download/Onco_genes_OncoVar_TCGA.tar.gz
tar -zxvf $DownloadTumorGeneFold"Onco_genes_OncoVar_TCGA.tar.gz" -C $DownloadTumorGeneFold 
gunzip $DownloadTumorGeneFold"Onco_genes_OncoVar_TCGA/"* 
wget -P $DownloadTumorGeneFold https://oncovar.org/resource/download/Onco_genes_OncoVar_ICGC.tar.gz
tar -zxvf $DownloadTumorGeneFold"Onco_genes_OncoVar_ICGC.tar.gz" -C $DownloadTumorGeneFold 
gunzip $DownloadTumorGeneFold"Onco_genes_OncoVar_ICGC/"*
#wget -P $DownloadTumorGeneFold https://oncovar.org/resource/download/Onco_mutations_OncoVar_TCGA.tar.gz
#tar -zxvf $DownloadTumorGeneFold"Onco_mutations_OncoVar_TCGA.tar.gz" -C $DownloadTumorGeneFold 
#gunzip $DownloadTumorGeneFold"Onco_mutations_OncoVar_TCGA/"*
#wget -P $DownloadTumorGeneFold https://oncovar.org/resource/download/Onco_mutations_OncoVar_TCGA.tar.gz
#tar -zxvf $DownloadTumorGeneFold"Onco_mutations_OncoVar_ICGC.tar.gz" -C $DownloadTumorGeneFold 
#gunzip $DownloadTumorGeneFold"Onco_mutations_OncoVar_ICGC/"*
wget -P $DownloadTumorGeneFold https://oncovar.org/resource/download/Onco_pathways_OncoVar_TCGA.tar.gz
tar -zxvf $DownloadTumorGeneFold"Onco_pathways_OncoVar_TCGA.tar.gz" -C $DownloadTumorGeneFold 
gunzip $DownloadTumorGeneFold"Onco_pathways_OncoVar_TCGA/"* 
wget -P $DownloadTumorGeneFold https://oncovar.org/resource/download/Onco_pathways_OncoVar_ICGC.tar.gz
tar -zxvf $DownloadTumorGeneFold"Onco_pathways_OncoVar_ICGC.tar.gz" -C $DownloadTumorGeneFold 
gunzip $DownloadTumorGeneFold"Onco_pathways_OncoVar_ICGC/"*
#
cp /data9/swu13/CAeditome_hg38/RawData/TumorGene/Onco_genes_OncoVar_ICGC/ICGC.BRCA.onco.genes.OncoVar.tsv /data9/swu13/CAeditome_hg38/RawData/TumorGene/Onco_genes_OncoVar_TCGA
cp /data9/swu13/CAeditome_hg38/RawData/TumorGene/Onco_genes_OncoVar_ICGC/ICGC.OV.onco.genes.OncoVar.tsv /data9/swu13/CAeditome_hg38/RawData/TumorGene/Onco_genes_OncoVar_TCGA
cp /data9/swu13/CAeditome_hg38/RawData/TumorGene/Onco_genes_OncoVar_ICGC/ICGC.PanCancer.onco.genes.OncoVar.tsv /data9/swu13/CAeditome_hg38/RawData/TumorGene/Onco_genes_OncoVar_TCGA
cp /data9/swu13/CAeditome_hg38/RawData/TumorGene/Onco_genes_OncoVar_ICGC/ICGC.PRAD.onco.genes.OncoVar.tsv /data9/swu13/CAeditome_hg38/RawData/TumorGene/Onco_genes_OncoVar_TCGA
head -n 1 $DownloadTumorGeneFold"Onco_genes_OncoVar_TCGA/"*"ACC"* > $DownloadTumorGeneFold"/GenesAnnotation"
for file in `ls $DownloadTumorGeneFold"Onco_genes_OncoVar_TCGA/"* `
do
awk -F"\t" '(NR>1){print $0}' $file >> $DownloadTumorGeneFold"/GenesAnnotation"
done

#TSGene2.0
TSGCancers=("BLCA" "BRCA" "COAD" "HNSC" "KICH" "KIRC" "LUAD" "LUSC" "PRAD" "THCA" "UCEC")
wget https://bioinfo.uth.edu/TSGene/sig_exp.txt?csrt=18162788718577482961 -O $DownloadTumorGeneFold"TSG2016.P.txt"
awk -F"\t" -v cancers="${TSGCancers[*]}" '(NR>1){split(cancers,a," ");for(i=3;i<=NF;i++){if($i<0.05){out=a[i-2]"\t"$2"\t.\t.\t"$1"\tY\t"$i"\t.\t.\t.\tTSG";for(i=12;i<=29;i++){out=out"\t."};print out}}}' $DownloadTumorGeneFold"TSG2016.P.txt" >> $DownloadTumorGeneFold"/GenesAnnotation"
#Possible tumor-related gene
head -n 1 $DownloadTumorGeneFold"/GenesAnnotation" | cut -f2 > $DownloadTumorGeneFold"/Genes"
cut -f2 $DownloadTumorGeneFold"/GenesAnnotation" | tail -n +2 | sort | uniq >>  $DownloadTumorGeneFold"/Genes"
} 

function DownloadGeneSummary(){
DownloadFold=$DatabaseFold"RawData/"
DownloadGeneSummary=$DownloadFold"GeneSummary/"
mkdir $DownloadGeneSummary
wget -O $DownloadGeneSummary"Ensemble.GeneSummary" 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "external_synonym" /><Attribute name = "description" /><Attribute name = "uniprotswissprot" /><Attribute name = "band" /></Dataset></Query>' 
wget -O $DownloadGeneSummary"goa_human_rna.gaf.gz" 'http://geneontology.org/gene-associations/goa_human_rna.gaf.gz'
gunzip $DownloadGeneSummary"goa_human_rna.gaf.gz"
wget -O $DownloadGeneSummary"goa_human.gaf.gz" 'http://geneontology.org/gene-associations/goa_human.gaf.gz'
gunzip $DownloadGeneSummary"goa_human.gaf.gz"
}

function DownloadmiRNAENSG(){
DownloadFold=$DatabaseFold"RawData/"
DownloadGeneSummary=$DownloadFold"Reference/"
wget -O $DownloadGeneSummary"Ensemble.miRNA" 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "mirbase_id" /></Dataset></Query>' 
}

function DownloadGTEX(){
DownloadFold=$DatabaseFold"RawData/"
DownloadGTEXFold=$DownloadFold"GTEX/"
mkdir $DownloadGTEXFold
wget -O $DownloadGTEXFold"GTEX_phenotype.gz" https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/GTEX_phenotype.gz
gunzip $DownloadGTEXFold"GTEX_phenotype.gz"
wget -O $DownloadGTEXFold"gtex_RSEM_gene_fpkm.gz" https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_RSEM_gene_fpkm.gz
gunzip $DownloadGTEXFold"gtex_RSEM_gene_fpkm.gz"
}