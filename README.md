# CAeditome
# CAeditome: functional annotation of A-to-I RNA editing in cancer

#This project contains the codes for diverse down-stream analyses of A-to-I RNA editing in human pan-cancers from TCGA

#1. doanload links related to the project. 
Tool: LiftOver
Download.sh
Download_Gene_expression_Data.R
ExtractUTR.R


#2. RNA editing detection. 
RNAdetection.sh

#3. RNA editing annotations in genes, regions, and repeats.
RNAeditingGeneralAnnotation.sh

#4. RNA editing frequencies analysis. 
SpecificRNAeditingAnnotation.sh
DifferentialEditing.R
StageCorrelation.R
SurvivalAnalysis.R
Roman2Num.R
Heatmap.R
generate_ucsc_bed.sh

#5. Correlations between RNA editng and gene expression (ADAR enzymes). 
ExpressionRNAeditingAnnotation.sh
ExpressionCorrelation.R
Heatmap.R
ADARCorrelation.R
DifferentialExpression.R
StageExpressionCorrelation.R
Roman2Num.R
SurvivalExpressionAnalysis.R

#6. annotation of RNA editing in protein-coding regions. 
RNAeditingProteinAnnotation.sh

#7. Correlation of RNA editing and Splicing. 
Tools: MaxEntScan
RNAeditingSplicing.sh
SplicingCorrelation.R
DifferentialSplicing.R
Heatmap.R

#8. Annotation of RNA editing effects on miRNA binding and regulation. 
Tools: TargetScan and miRanda
RNAeditingMiRNA.sh
CompeteComputeTwo.R
Network.R

#9. Annotation of RNA editing effects on RNA structure.
Tool: RNAfold
RNAeditingStructureAnnotation.sh

#10. Statistics including all criteria for each part. 
Statistic.sh
BarPlot.R
Boxplot.R
BubblePlot.R
Enrichment.R
Histogram.R
VennPlot.R
Violin.R
VolcanoPlot.R

#11. Analyses results. 
https://ccsm.uth.edu/CAeditome/index.html

#12. citation. 
Wu S, Fan Z, Kim P, Huang L, and Zhou X, CAeditome: functional annotation of A-to-I RNA editing in patients with cancer.
