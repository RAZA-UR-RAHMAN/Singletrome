library("ggpubr")
library(ggplot2)
library(plyr)
library(ggstatsplot)
source("/data/Mullen_1/Raza/scRNA/Analysis/Scripts/Pipeline_scripts/CustomTheme.R")
source("/data/Mullen_1/Raza/scRNA/Analysis/Scripts/Pipeline_scripts/CustomTheme.R")


################################################# FUNCTIONS ################################################# 
obtainGeneAndTranscriptMatrix<-function(geneLengthPath,transcriptLengthPath){
  Genes<-read.table(geneLengthPath,sep = "\t",header = T, row.names = 1)
  Transcripts<-read.table(transcriptLengthPath,sep = "\t",header = T, row.names = 1)
  colnames(Transcripts)[5]<-"Transcript_Length"
  colnames(Genes)<-"Gene_Length"
  Genes$gene_id=rownames(Genes)
  GeneAndTranscriptLengthsDF<-merge(Genes,Transcripts,by="gene_id")
  return(GeneAndTranscriptLengthsDF)
}
################################################# 
obtainGeneToPlots<-function(bedFilePath,GeneAndTranscriptLengthsDF){
  bedFile<-read.table(bedFilePath,header = F,sep = "\t")
  genesToConsider<-GeneAndTranscriptLengthsDF[which(GeneAndTranscriptLengthsDF$tx_name %in% bedFile$V4),"gene_id"]
  tempMatrix<-GeneAndTranscriptLengthsDF[which(GeneAndTranscriptLengthsDF$gene_id %in% genesToConsider),]
  return(tempMatrix)
}
################################################# PLOT
plotScatter<-function(df,outputDir,fileName){
  df$Transcript_Length=df$Transcript_Length/1000
  df$Gene_Length=df$Gene_Length/1000
  p <- ggscatter(df, x = "Gene_Length", y = "Transcript_Length", add = "reg.line",cor.method = "pearson",
                 title = fileName ,xlab = "Gene length (kb)", ylab = "Transcript length (kb)", color = "type",shape = ".",palette = c("#386cb0","#fdb462"),fullrange = TRUE) + 
    stat_cor(aes(color = type), label.x = 3) 
  pdf(paste0(outputDir,"/",fileName,"_scatter.pdf"),width = 3,height = 3)
  #png(paste0(outputDir,"/",fileName,"_scatter.png"))
  print(p)
  dev.off()
}
################################################# 
plotViolinPlot<-function(df,outputDir,fileName){
  p <- ggplot(df, aes(x=type, y=Transcript_Length,fill=type)) + geom_violin()
  p <- p  +theme_Publication() + scale_fill_Publication() + scale_x_discrete(limits=c("proteincoding", "lncRNA")) + scale_fill_manual(values=c("#fdb462", "#386cb0")) +
    geom_boxplot(width=0.1, fill="white")+
    labs(title="Transcript lengths",x="", y = "Length (Kb)") + theme_classic()
  pdf(paste0(outputDir,"/",fileName,"_Violin.pdf"),width = 4,height = 4)
  print(p)
  dev.off()
}
################################################# VARIABLES ################################################# 
outputDir="/data/Mullen_1/Raza/scRNA/Genomes/GTFDetails/LengthAnalysis/"
bedFilesDir="/data/Mullen_1/Raza/scRNA/Genomes/GTFDetails/bedFiles/"
if(!dir.exists(outputDir)){dir.create(outputDir,mode = "0777")}
pcgGeneLengthPath="/data/Mullen_1/Raza/scRNA/Genomes/GTFDetails/gencode_v32_geneLength.txt"
pcgTranscriptLengthPath="/data/Mullen_1/Raza/scRNA/Genomes/GTFDetails/gencode_v32_transcriptLength.txt"
lncRNAGeneLengthPath="/data/Mullen_1/Raza/scRNA/Genomes/GTFDetails/LncBook_v2_geneLength.txt"
lncRNATranscriptLengthPath="/data/Mullen_1/Raza/scRNA/Genomes/GTFDetails/LncBook_v2_transcriptLength.txt"
pcgGeneAndTranscriptLengthsDF<-obtainGeneAndTranscriptMatrix(pcgGeneLengthPath,pcgTranscriptLengthPath)
lncRNAAndTranscriptLengthsDF<-obtainGeneAndTranscriptMatrix(lncRNAGeneLengthPath,lncRNATranscriptLengthPath)
################################################################################################## ANALYSIS ################################################################################################# 
################################################# All Genes to Transcripts ####################################### 
allPCGCorrlation<-round(cor(pcgGeneAndTranscriptLengthsDF$Gene_Length, pcgGeneAndTranscriptLengthsDF$Transcript_Length , method = c("pearson")),digits = 2)
allLncRNACorrelation<-round(cor(lncRNAAndTranscriptLengthsDF$Gene_Length, lncRNAAndTranscriptLengthsDF$Transcript_Length , method = c("pearson")),digits = 2)
print(paste0("allPCGCorrlation : ",allPCGCorrlation," allLncRNACorrelation: ",allLncRNACorrelation))
df<-rbind(pcgGeneAndTranscriptLengthsDF,lncRNAAndTranscriptLengthsDF)
#plotViolinPlot(df,outputDir,fileName)
plotScatter(df,outputDir,"All_Genes")
plotScatter(pcgGeneAndTranscriptLengthsDF,outputDir,"All_PCG")
plotScatter(lncRNAAndTranscriptLengthsDF,outputDir,"All_lncRNA")

################################################# GOING THROUGH TRANSCRIPT LENGTH FOR THE GENE TO MATCH RESEQC PLOTS FOR REASONING  ####################################### 
#TAKE ALL THE TRNACRIPTS OF DIFFERENT LENGTH FOR EACH GENE. HOWEVER THE GENES ARE CHOSEN BASED ON THE TRANSCRIPT LENTHS INITIALLY
TRANSCRIPT_TO_GENEOutputDir="/data/Mullen_1/Raza/scRNA/Genomes/GTFDetails/LengthAnalysis/TRANSCRIPT_TO_GENE/"
if(!dir.exists(TRANSCRIPT_TO_GENEOutputDir)){dir.create(TRANSCRIPT_TO_GENEOutputDir,mode = "0777")}
rangesDf<-read.table("/data/Mullen_1/Raza/scRNA/Genomes/GTFDetails/transcriptRanges.txt",sep="\t",header = T)
rangesTemp<-unique(rangesDf[,c("lowerLimit","upperLimit")])
correlationDf<-data.frame(lowerLimit=0,UpperLimit=0,PCG_Corrlation=0,lncRNA_Corrlation=0)
for(i in 1:(nrow(rangesTemp)-1)){
  lowerLimit=rangesTemp[i,"lowerLimit"]
  PCGUpperLimit=rangesTemp[i,"upperLimit"]
  lncRNAUpperLimit=rangesTemp[i,"upperLimit"]
  upperLimitLabel=rangesDf[i,"upperLimit"]
  if(lowerLimit==50000){PCGUpperLimit=109224;lncRNAUpperLimit=205012;upperLimitLabel="100000+"}
  
  lncRNABedFileName=paste0(bedFilesDir,rangesDf[which(rangesDf$lowerLimit==lowerLimit & rangesDf$upperLimit==lncRNAUpperLimit & rangesDf$type=="lncRNA"),"fileName"],".bed")
  pcgBedFileName=paste0(bedFilesDir,rangesDf[which(rangesDf$lowerLimit==lowerLimit & rangesDf$upperLimit==PCGUpperLimit & rangesDf$type=="proteincoding"),"fileName"],".bed")
  print(paste0("--- Running for lncRNABedFileName:",lncRNABedFileName, " pcgBedFileName:", pcgBedFileName))
  # Get matrices
  lncRNATempMatrix<-obtainGeneToPlots(lncRNABedFileName,lncRNAAndTranscriptLengthsDF)
  pcgTempMatrix<-obtainGeneToPlots(pcgBedFileName,pcgGeneAndTranscriptLengthsDF)
  
  fileName=paste0(lowerLimit,"-",upperLimitLabel)
  df<-rbind(pcgTempMatrix,lncRNATempMatrix)
  plotScatter(df,TRANSCRIPT_TO_GENEOutputDir,paste0(fileName,"_All"))
  plotScatter(pcgTempMatrix,TRANSCRIPT_TO_GENEOutputDir,paste0(fileName,"_PCG"))
  plotScatter(lncRNATempMatrix,TRANSCRIPT_TO_GENEOutputDir,paste0(fileName,"_lncRNA"))
  
  # Correlation
  PCGCorrlation<-round(cor(pcgTempMatrix$Gene_Length, pcgTempMatrix$Transcript_Length , method = c("pearson")),digits = 2)
  LncRNACorrelation<-round(cor(lncRNATempMatrix$Gene_Length, lncRNATempMatrix$Transcript_Length , method = c("pearson")),digits = 2)
  correlationDf<-rbind(correlationDf,c(lowerLimit,upperLimitLabel,PCGCorrlation,LncRNACorrelation))
}
correlationDf=correlationDf[-(1:3),]
write.table(file=paste0(TRANSCRIPT_TO_GENEOutputDir,"Transcript_To_Gene_To_All_Transcripts_Correlation.txt"),correlationDf,sep = "\t",row.names = F,quote = F)


################################################# GENE LENGTH ####################################### 
geneRangesOUTPUT="/data/Mullen_1/Raza/scRNA/Genomes/GTFDetails/LengthAnalysis/GENE_RANGES/"
if(!dir.exists(geneRangesOUTPUT)){dir.create(geneRangesOUTPUT,mode = "0777")}
geneSizes<-seq(1000,20000,1000)
for(i in 1:(length(geneSizes))){
  j=i+1
  lowerLimit=geneSizes[i]
  upperLimit=geneSizes[j]
  if(is.na(upperLimit)){
    upperLimit=max(c(pcgGeneAndTranscriptLengthsDF$Gene_Length,lncRNAAndTranscriptLengthsDF$Gene_Length))
  }
  # pcg
  pcgGenesPerLimits=pcgGeneAndTranscriptLengthsDF[which(pcgGeneAndTranscriptLengthsDF$Gene_Length>lowerLimit & pcgGeneAndTranscriptLengthsDF$Gene_Length<=upperLimit ),]
  pcgNGenesInRange<-length(unique(pcgGenesPerLimits$gene_id))
  pcgGenesFileName=paste0("PCG_",lowerLimit,"_",upperLimit,"_",pcgNGenesInRange)
  plotCorrelation(pcgGenesPerLimits,geneRangesOUTPUT,pcgGenesFileName)
  # lncRNA
  lncRNAGenesPerLimits=lncRNAAndTranscriptLengthsDF[which(lncRNAAndTranscriptLengthsDF$Gene_Length>lowerLimit & lncRNAAndTranscriptLengthsDF$Gene_Length<=upperLimit ),]
  lncRNANGenesInRange<-length(unique(lncRNAGenesPerLimits$gene_id))
  lncRNAGenesFileName=paste0("lncRNA_",lowerLimit,"_",upperLimit,"_",lncRNANGenesInRange)
  plotCorrelation(lncRNAGenesPerLimits,geneRangesOUTPUT,lncRNAGenesFileName)
}






# b <- width(transcriptsBy(TxDb,by='gene'))
# c <- transcriptsBy(TxDb,by='gene')
# txlens <- transcriptLengths(TxDb,with.cds_len=T)

txdb=TxDb
txlens <- transcriptLengths(txdb)
exons.list.per.gene <- transcriptsBy(txdb,by="gene")

getGeneLength<-function(outputDir,type,txdb){
  #Tutorial:Extract Total Non-Overlapping Exon Length Per Gene With Bioconductor (https://www.biostars.org/p/83901/)
  exons.list.per.gene <- exonsBy(txdb,by="gene")
  exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))
  colnames(exonic.gene.sizes)="Length"
  write.table(exonic.gene.sizes,paste0(outputDir,type,"_geneLength.txt"),quote = F,sep="\t",col.names=NA)
  return(exonic.gene.sizes)
}
##########################################################
getTranscriptLength<-function(outputDir,type,txdb){
  txlens <- transcriptLengths(txdb)
  colnames(txlens)[5]="Length"
  write.table(txlens,paste0(outputDir,type,"_transcriptLength.txt"),quote = F,sep="\t",col.names=NA)
  return(txlens)
}
