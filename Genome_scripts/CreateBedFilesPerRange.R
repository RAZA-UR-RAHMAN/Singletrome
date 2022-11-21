library(GenomicFeatures)
library(ggplot2)
require(gridExtra)

##########################################################Functions##########################################################
getGeneLength<-function(outputDir,type,txdb){
#Tutorial:Extract Total Non-Overlapping Exon Length Per Gene With Bioconductor (https://www.biostars.org/p/83901/)
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))
colnames(exonic.gene.sizes)="Gene_Length"
exonic.gene.sizes$gene_id=rownames(exonic.gene.sizes)
write.table(exonic.gene.sizes,paste0(outputDir,type,"_geneLength.txt"),quote = F,sep="\t",col.names=NA)
return(exonic.gene.sizes)
}
##########################################################
getTranscriptLength<-function(outputDir,dataSource,type,txdb){
  txlens <- transcriptLengths(txdb)
  colnames(txlens)[5]="Transcript_Length"
  txlens$type=type
  write.table(txlens,paste0(outputDir,dataSource,"_transcriptLength.txt"),quote = F,sep="\t",col.names=NA)
  return(txlens)
}
##########################################################
makeTxDbFromGTF<-function(file,format,dataSource,organism,taxonomyId){
  TxDb=makeTxDbFromGFF(file, format=format, dataSource=dataSource, organism=organism, taxonomyId=taxonomyId)
  return(TxDb)
}
##########################################################
writeBedFilesPerRange<-function(transcripts,bed12File,dataSource,type,featureSize,bedFilesDir){
  numberOfFeaturesDF <- data.frame(matrix(ncol = 5, nrow = 0))
  bed12Df<-read.table(bed12File,sep=" ",header = F)
  for(i in 1:(length(featureSize))){
    j=i+1
    lowerLimit=featureSize[i]
    upperLimit=featureSize[j]
    if(is.na(upperLimit)){
      upperLimit=max(transcripts$Transcript_Length)
    }
    print(paste0(lowerLimit,"   ",upperLimit))
    transcriptsPerLimits=transcripts[which(transcripts$Transcript_Length>lowerLimit & transcripts$Transcript_Length<=upperLimit ),]

    tempDf=bed12Df[which(bed12Df$V4 %in% transcriptsPerLimits$tx_name),]
    if(nrow(tempDf)!=nrow(transcriptsPerLimits)){stop(paste0("Not all transcripts per limit ",nrow(transcriptsPerLimits)," are found in the orignal bed file ",nrow(tempDf)))}
    fileName=paste0(dataSource,"_",lowerLimit,"_",upperLimit,"_",nrow(tempDf))
    outBedFileName=paste0(bedFilesDir,"/",fileName,".bed")
    write.table(tempDf,file=outBedFileName,sep="\t",quote = F,col.names = F,row.names = F)
    # just to keep track and may be plot letter
    row=c(lowerLimit=lowerLimit,upperLimit=upperLimit,NumberOfFeatures=nrow(transcriptsPerLimits),type=type,fileName=fileName)
    numberOfFeaturesDF=rbind(numberOfFeaturesDF,row)
  }
  colnames(numberOfFeaturesDF) <- c("lowerLimit", "upperLimit", "NumberOfFeatures","type","fileName")
  write.table(numberOfFeaturesDF,file=paste0(bedFilesDir,"/",dataSource,"_numbers.txt"),sep="\t",quote = F,row.names = F)
  return(numberOfFeaturesDF)
}
################################################### MAIN
configDF<-read.table("/data/Mullen_1/Raza/scRNA/Analysis/Configs/createBedFile_config.txt",sep="\t",header = T)
outputDir="/data/Mullen_1/Raza/scRNA/Genomes/GTFDetails/"
setwd(outputDir) # it downloads the database locally to the current working dir
bedFilesDir=paste0(outputDir,"/bedFiles")
if(!dir.exists(bedFilesDir)){dir.create(bedFilesDir,mode = "0777")}
featureSize<-c(0,100,200,300,400,500,1000,2000,3000,4000,5000,10000,15000,20000,30000,40000,50000)
transcriptRanges<-data.frame()
transcriptlengthTable<-data.frame()
for(i in 1:nrow(configDF)){
  dataSource=configDF$dataSources[i]
  bed12File=configDF$bed12Files[i]
  type=configDF$type[i]

  TxDb=makeTxDbFromGTF(configDF$files[i], format="gtf", dataSource=dataSource, organism="Homo sapiens", taxonomyId=9606)

  print('Running....')
  ### GENE LENGTH
  genes=getGeneLength(outputDir,dataSource,TxDb)
  ### TRANSCRIPT LENGTH
  transcripts=getTranscriptLength(outputDir,dataSource,type,TxDb)
  transcripts<-merge(genes,transcripts,by="gene_id")
  # MERGE FOR ANALYSIS LATER
  tempRangesDf=writeBedFilesPerRange(transcripts,bed12File,dataSource,type,featureSize,bedFilesDir)
  transcriptRanges=rbind(transcriptRanges,tempRangesDf)
  transcriptlengthTable=rbind(transcriptlengthTable,transcripts)
}
write.table(transcriptRanges,file=paste0(outputDir,"/transcriptRanges.txt"),sep="\t",quote = F,row.names = F)
## MERGE GENE NAMES TO TRANSCRIPTS FOR LATER ON
geneIdsToGeneNamesLncPedia<-read.table("/data/Mullen_1/Raza/scRNA/Genomes/reference_sources/GENES/lncPediaGeneIdsAndNames.txt",sep = ",")
geneIdsToGeneNamesLncPedia<-geneIdsToGeneNamesLncPedia[which(geneIdsToGeneNamesLncPedia$V2!=""),]
geneIdsToGeneNamesGENECODE<-read.table("/data/Mullen_1/Raza/scRNA/Genomes/reference_sources/GENES/defaultFrom10XGeneIdsAndNames.txt",sep = ",")
geneIdsToGeneNamesGENECODE<-geneIdsToGeneNamesGENECODE[which(geneIdsToGeneNamesGENECODE$V2!=""),]
geneNamesToGeneIds<-rbind(geneIdsToGeneNamesGENECODE,geneIdsToGeneNamesLncPedia)
colnames(geneNamesToGeneIds)<-c("geneID","geneName")
##transcript length Table with gene names
transcriptAndGeneDf=merge(transcriptlengthTable,geneNamesToGeneIds,by.x="gene_id",by.y="geneID",all.x =T)
rownames(transcriptAndGeneDf)<-transcriptAndGeneDf[,"tx_name"]
write.table(transcriptAndGeneDf,file=paste0(outputDir,"/transcriptAndGene.txt"),sep="\t",quote = F,col.names = NA)
print("-----------------DONE------------------")
