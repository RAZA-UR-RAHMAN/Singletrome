set.seed(349870519)
library(Seurat)
library(cowplot)
library(dplyr)
library(patchwork)
# ************************************************************************************ RSEQC  ****************************************************************************** 
############################################################ ObtainCovPerTranscriptBodyWithLengthDF ###########################################################
ObtainCovPerTranscriptBodyWithLengthDF<-function(rseqQcDir,rangesDf,transcriptlengthTable,filteredTranscriptsDIR){
  covPerTranscriptBodyDF<-data.frame()
  for(i in 1:nrow(rangesDf)){
    fileName=rangesDf[i,"fileName"]
    covPerTranscriptBodyFile=paste0(rseqQcDir,fileName,".perGeneBodyCoverage.txt")
    if(!file.exists(covPerTranscriptBodyFile)){
      stop(paste0("One of the files does not exists",covPerTranscriptBodyFile))
    }
    message(paste0("Reading file ",covPerTranscriptBodyFile))
    tempDf=read.table(covPerTranscriptBodyFile,sep=",",header = T,row.names = 1)
    tempDf[,"fileName"]=fileName
    #tempDf[,"type"]=rangesDf[i,"type"]
    tempDf[,"lowerLimit"]=rangesDf[i,"lowerLimit"]
    tempDf[,"upperLimit"]=rangesDf[i,"upperLimit"]
    tempDf[,"NumberOfTranscriptsInFile"]=rangesDf[i,"NumberOfFeatures"]
    if(nrow(covPerTranscriptBodyDF)==0){
      covPerTranscriptBodyDF=tempDf
    }
    else{
      covPerTranscriptBodyDF<-rbind(covPerTranscriptBodyDF,tempDf)
    }
  }
  rm('tempDf')
  #### Remove transcripts with no expression at all. all zero rows
  covPerTranscriptBodyDF=covPerTranscriptBodyDF[rowSums(covPerTranscriptBodyDF[,paste0("percent_",1:100)])>0,]
  
  ## merge with transcriptlengthTable for gene ids and exact transcript length
  covPerTranscriptBodyWithLengthDF=merge(covPerTranscriptBodyDF,transcriptlengthTable,by="row.names")
  rownames(covPerTranscriptBodyWithLengthDF)<-covPerTranscriptBodyWithLengthDF$Row.names
  covPerTranscriptBodyWithLengthDF$Row.names<-NULL
  if(nrow(covPerTranscriptBodyDF)!=nrow(covPerTranscriptBodyWithLengthDF)){stop("In merging can not find all the transcripts")}
  write.table(file=paste0(filteredTranscriptsDIR,"RSEQC_COVERAGE_MATRIX.txt"),covPerTranscriptBodyWithLengthDF,sep="\t",col.names = NA,quote = F)
  rm('covPerTranscriptBodyDF')
  return(covPerTranscriptBodyWithLengthDF)
}
############################################################ removeTranscriptsFromCoverageMatrix ###########################################################
removeTranscriptsFromCoverageMatrix<-function(filteredTranscriptsForReadDistDirPerDs,outputFileName,covPerTranscriptBodyWithLengthDF,deleteTranscriptsDf){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  nTranscriptsBeforeDeletion=length(unique(rownames(covPerTranscriptBodyWithLengthDF)))
  TranscriptsToDelete<-rownames(deleteTranscriptsDf)
  covPerTranscriptBodyWithLengthDF<-covPerTranscriptBodyWithLengthDF[which(rownames(covPerTranscriptBodyWithLengthDF) %!in% TranscriptsToDelete),]
  nTranscriptsAfterDeletion=length(unique(rownames(covPerTranscriptBodyWithLengthDF)))
  if(length(TranscriptsToDelete)+nTranscriptsAfterDeletion!=nTranscriptsBeforeDeletion){stop("All the marked transcripts were not deleted from the coverage matrix")}
  write.table(file=paste0(filteredTranscriptsForReadDistDirPerDs,outputFileName,".txt"),covPerTranscriptBodyWithLengthDF,sep="\t",col.names = NA,quote = F)
  return(covPerTranscriptBodyWithLengthDF)
}
############################################################ check If All Expressed Transcripts Are Marked To be deleted #################################################
checkIfAllExpressedTranscriptsAreMarked<-function(markedDf,transcriptlengthTable,covPerTranscriptBodyWithLengthDF,filteredTranscriptsForReadDistDirPerDs,outputFileName){
  transcriptsToGeneDf<-data.frame(gene_id=NA,geneName=NA,TotalTranscriptInGene=0,transcriptsFoundinRseq=0,transcriptsFoundToFilter=0,delete=F)
  for(r in 1:nrow(markedDf)){
    g=markedDf[r,"gene_id"]
    geneName=markedDf[r,"geneName"]
    delete=F
    allTranscriptsForGene<-rownames(transcriptlengthTable)[which(transcriptlengthTable$gene_id==g)]
    transcriptsFoundinRseq<-rownames(covPerTranscriptBodyWithLengthDF)[which(rownames(covPerTranscriptBodyWithLengthDF) %in% allTranscriptsForGene)]
    transcriptsFoundToFilter=rownames(markedDf)[which(rownames(markedDf) %in% transcriptsFoundinRseq)]
    if(length(transcriptsFoundinRseq)==length(transcriptsFoundToFilter)){
      delete=T
    }
    r=c(g,geneName,length(allTranscriptsForGene),length(transcriptsFoundinRseq),length(transcriptsFoundToFilter),delete)
    transcriptsToGeneDf<-rbind(transcriptsToGeneDf,r)
  }
  transcriptsToGeneDf=unique(transcriptsToGeneDf)
  transcriptsToGeneDf=transcriptsToGeneDf[order(transcriptsToGeneDf$delete),]
  transcriptsToGeneDf=transcriptsToGeneDf[-1,] # REMOVE THE NA ROW
  write.table(file=paste0(filteredTranscriptsForReadDistDirPerDs,outputFileName,".txt"),transcriptsToGeneDf,sep="\t",col.names = NA,quote = F)
  return(transcriptsToGeneDf)
}
############################################################ delete Genes From Singletrome Based On Qc #################################################
deleteGenesFromSingletromeBasedOnQc<-function(deleteFromSingletrome,seuObj){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  deleteFromSingletrome=deleteFromSingletrome[which(deleteFromSingletrome$DeletePeakBin=="PeakBin" | deleteFromSingletrome$DeleteFivePrime=="FivePrime"),]
  nGenesToDelete<-length(unique(deleteFromSingletrome$gene_id))
  nGenesInSingletromeBeforeDeletion<-nrow(seuObj)
  # delete
  seuObj=seuObj[(rownames(seuObj) %!in% deleteFromSingletrome$gene_id &  rownames(seuObj) %!in% deleteFromSingletrome$geneName
                 &  rownames(seuObj) %!in% paste0(deleteFromSingletrome$geneName,".1")
                 &  rownames(seuObj) %!in% paste0(deleteFromSingletrome$geneName,".2")
                 &  rownames(seuObj) %!in% paste0(deleteFromSingletrome$geneName,".3")
                 &  rownames(seuObj) %!in% paste0(deleteFromSingletrome$geneName,".4")
                 &  rownames(seuObj) %!in% paste0(deleteFromSingletrome$geneName,".5")
                 &  rownames(seuObj) %!in% paste0(deleteFromSingletrome$geneName,".6")),]
  nGenesInSingletromeAfterDeletion<-nrow(seuObj)
  if(nGenesToDelete+nGenesInSingletromeAfterDeletion!=nGenesInSingletromeBeforeDeletion){stop("All the genes were not sucessfully deleted from Singletrome")}
  return(seuObj)
}

#   ******************************************************************** SEURAT ANALYSIS  ******************************************************************** 
############################################################ validateSampleIDToBarcode ###########################################################
validateSampleIDToBarcode<-function(seuratObj){
  tempObject=seuratObj
  tempObject@meta.data$testID = lapply(rownames(tempObject@meta.data), function(e){unlist(strsplit(e,split = "_"))[[1]]})
  which(tempObject@meta.data$SampleID!=tempObject@meta.data$testID)
  print(length(which(tempObject@meta.data$SampleID==tempObject@meta.data$testID)))
  print(nrow(tempObject@meta.data))
  outcome=nrow(tempObject@meta.data)==length(which(tempObject@meta.data$SampleID==tempObject@meta.data$testID))
  rm(tempObject)
  return(outcome)
}

###############################################################loadSeuratObjectsFromCellRanger##########################################################################
loadSeuratObjectsFromCellRanger<-function(samples,datasetID,referenceName){
  if(length(samples)==1){
    s=samples
    samplePath=paste0(mappingDir,datasetID,"/",referenceName,"/",s,"/",s,'/outs/raw_feature_bc_matrix/')
    cur_data <- Read10X(samplePath)
    seuratObj <- CreateSeuratObject(cur_data,min.cells = 10)
    }else{
      seurat_list <- sapply(samples, function(s){
        samplePath=paste0(mappingDir,datasetID,"/",referenceName,"/",s,"/",s,'/outs/raw_feature_bc_matrix/')
        cur_data <- Read10X(samplePath)
        if(datasetID=="GSE136103"){s=str_replace(s,"_Plus", "+")} # to match the barcodes of +
        colnames(cur_data) <- paste(s,sapply(colnames(cur_data), '[[' , 1L ), sep="_")
        cur_data
        })
      experiment.data <- do.call("cbind", seurat_list)
      seuratObj <- CreateSeuratObject(experiment.data,min.cells = 10,names.field = 1)
      seuratObj$SampleID<-seuratObj$orig.ident
      if(validateSampleIDToBarcode(seuratObj)==FALSE){stop("Error: Samples and Barcodes not properly merged!")}
    }
  return(seuratObj)
}
#########################################################PerformDownStreamAnalysis################################################################################
PerformDownStreamAnalysis<-function(seuratObjTemp,currentOutDir,fileName){
  seuratObjTemp=PerformSeuratAnalysis(seuratObjTemp,currentOutDir,fileName)
  writeGenesExpressed(seuratObjTemp,currentOutDir,fileName)
  return(seuratObjTemp)
}
########################################################PerformSeuratAnalysis#################################################################################
PerformSeuratAnalysis<-function(seuratObjTemp,currentOutDir,fileName){
  seuratObjTemp<-Seurat::NormalizeData(seuratObjTemp,verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = seuratObjTemp@var.genes, npcs = 30, verbose = FALSE)%>%
    RunUMAP(reduction = "pca", dims = 1:30)
  saveRDS(seuratObjTemp,file = paste0(currentOutDir,fileName,".RDS"))
  pdf(paste0(currentOutDir,fileName,"_UMAP.pdf"),width = 12)
  p1<-DimPlot(seuratObjTemp, reduction = "umap", group.by='celltype',label = TRUE) +  ggtitle(fileName) +umap_theme
  print(p1)
  dev.off()
  return(seuratObjTemp)
}
###########################################################writeGenesExpressed##############################################################################
writeGenesExpressed<-function(seuratObjTemp,currentOutDir,fileName){
  counts <- GetAssayData(object = seuratObjTemp, slot = "counts")
  Keep_genes <- Matrix::rowSums(counts) >= 1
  Keep_genes_True=Keep_genes[which(Keep_genes==TRUE)]
  #length(unique(names(Keep_genes_True)))
  write.table(file=paste0(currentOutDir,fileName,"_genes_expressed.txt"),unique(names(Keep_genes_True)),quote = F,row.names = F,col.names = F)
}
########################################################################## Differential Expression ##################################################### 
plotGenes<-function(seuratPerCellType,markersForCelltype,currentDir,fileName,performDEOnColumn){
  topN=20
  g=rownames(head(markersForCelltype, n = topN))
  pdf(paste0(currentDir,"DEG_Top_",topN,"_",fileName,".pdf"), width=12,height = 14)
  plots <- VlnPlot(seuratPerCellType, features = g, split.by = performDEOnColumn, group.by = 'celltype', pt.size = 0.1, combine = FALSE) 
  print(wrap_plots(plots = plots, ncol = 4))
  dev.off()
}

########################################################################## Identify markers for each cell type ####################################### 
findAllMarkersForCellType<-function(seuratObj,N,outDir,fileName){
  #Note that the raw and normalized counts are stored in the counts and data slots of RNA assay. By default, the functions for finding markers will use normalized data.
  DefaultAssay(seuratObj) <- "RNA"
  seuratObj<- ScaleData(object = seuratObj, features = rownames(seuratObj))
  Idents(seuratObj)=seuratObj@meta.data$celltype
  seuratObj@misc$markers <- FindAllMarkers(object = seuratObj, only.pos = TRUE,logfc.threshold = 0.25,min.pct =  0.25) # ,min.diff.pct = 0.25
  markersForCelltype <- as.data.frame(seuratObj@misc$markers)
  write.table(markersForCelltype,paste0(outDir,"/",fileName,"_celltype_markers.txt"),  sep = '\t', col.names = NA, quote = F)
  topN <- markersForCelltype %>% group_by(cluster) %>% top_n(N, avg_log2FC)
  pdf(paste0(outDir,"/",fileName,"_celltype_markers_top_",N,".pdf"), width=7,height = 6) # main Figure
  print(DoHeatmap(object = seuratObj, features = topN$gene, label = TRUE,size = 4.3) + theme(axis.text.y = element_text(size = 12)))
  dev.off()
  return(seuratObj)
}
####
getNMarkersForCellTypes<-function(seuratObj,N){
  markersForCelltype <- as.data.frame(seuratObj@misc$markers)
  topN <- markersForCelltype %>% group_by(cluster) %>% top_n(N, avg_log2FC)
  return(topN)
}

########################################################################## Identify differential expressed genes across conditions ######################### 
### Identify differential expressed genes across conditions. we are confident in having identified common cell types across condition, we can ask what genes change in different conditions for cells of the same type
findDEGInCelltypesAcrossConditions<-function(seuratObj,performDEOnColumn,outDir,proteinCodingGenes,singletromeLncRNAs){
  destatsDf<-data.frame(celltype=NA,group1Name=NA, group2Name=NA,nCellsGroup1=0,nCellsGroup2=0,TotalDE=0,TotalUpRegulated=0,TotalDownRegulated=0,
                        nPcgUpRegulated=0,nPcgDownRegulated=0,nLncRNAUpRegulated=0,nLncRNADownRegulated=0)
  m=seuratObj@meta.data
  m$deColumn=paste0(m$celltype,sep="_",m[[performDEOnColumn]])
  seuratObj@meta.data=m
  Idents(seuratObj)=seuratObj@meta.data$deColumn
  deColumnUniqueValues=unique(m[[performDEOnColumn]])
  print(deColumnUniqueValues)
  for(ct in unique(seuratObj@meta.data$celltype)){
    currentDir=paste0(outDir,"/",ct,"/")
    dir.create(currentDir,recursive = T)
    group1=paste0(ct,"_",deColumnUniqueValues[2])
    group2=paste0(ct,"_",deColumnUniqueValues[1])
    fileName=paste0(group1,"_vs_",group2)
    print(fileName)
    seuratPerCellType <- subset(x = seuratObj, subset = (deColumn == group1 | deColumn == group2))
    
    possibleError <- tryCatch({
      dePerCellType<- FindMarkers(seuratPerCellType, ident.1 = group1, ident.2 = group2)
      write.table(dePerCellType,paste0(currentDir,fileName,"_unfiltered_de.txt"),  sep = '\t', col.names = NA, quote = F)
      # filtered for p_val_adj<=0.1 and abs(dePerCellType$avg_log2FC)>0.25
      dePerCellType=dePerCellType[which(dePerCellType$p_val_adj<=0.1 & abs(dePerCellType$avg_log2FC)>0.25),]   
      dePerCellType=dePerCellType %>% tibble::rownames_to_column() %>% arrange(desc(abs(avg_log2FC))) 
      rownames(dePerCellType)<-dePerCellType$rowname
      dePerCellType$rowname=NULL
      
      upRegulated=dePerCellType[dePerCellType$avg_log2FC>0,]
      downRegulated=dePerCellType[dePerCellType$avg_log2FC<0,]
      write.table(dePerCellType,paste0(currentDir,fileName,".txt"),  sep = '\t', col.names = NA, quote = F)
      
      # save number of DEGs for the cell type
      biotypeStats<-getStatsPerBiotype(upRegulated,downRegulated,proteinCodingGenes,singletromeLncRNAs)
      tempDf<-c(celltype=ct,group1Name=group1, group2Name=group2,nCellsGroup1=length(which(seuratObj$deColumn==group1)),nCellsGroup2=length(which(seuratObj$deColumn==group2)),
                TotalDE=nrow(dePerCellType),TotalUpRegulated=nrow(upRegulated),TotalDownRegulated=nrow(downRegulated),biotypeStats)
      destatsDf=rbind(destatsDf,tempDf)
      # plots
      plotGenes(seuratPerCellType,upRegulated,currentDir,paste0(fileName,"_upReg_"),performDEOnColumn)
      plotGenes(seuratPerCellType,downRegulated,currentDir,paste0(fileName,"_downReg_"),performDEOnColumn)
      # write the object , if need to plot differently
      #saveRDS(seuratPerCellType,file = paste0(currentDir,ct,"_.RDS"))
    },error=function(e) print(e))
  }
  destatsDf=destatsDf[-1,] # REMOVE THE FIRST EMPTY ROW
  write.table(destatsDf,paste0(outDir,"/DIFFERNETIAL_EXPRESSION_STATS.txt"),  sep = '\t', row.names = F, quote = F)
}
########################################################################## getStatsPerBiotype  ######################### 
getStatsPerBiotype<-function(upRegulated,downRegulated,proteinCodingGenes,singletromeLncRNAs){
  nPcgUpRegulated<-length(which(rownames(upRegulated) %in% proteinCodingGenes))
  nLncRNAUpRegulated<-length(which(rownames(upRegulated) %in% singletromeLncRNAs))
  
  nPcgDownRegulated<-length(which(rownames(downRegulated) %in% proteinCodingGenes))
  nLncRNADownRegulated<-length(which(rownames(downRegulated) %in% singletromeLncRNAs))
  
  if(nrow(upRegulated)!=(nPcgUpRegulated+nLncRNAUpRegulated)){stop("ERROR! ALL UP-REGULATED CAN NOT BE ASSIGNED A BIOTYPE")}
  if(nrow(downRegulated)!=(nPcgDownRegulated+nLncRNADownRegulated)){stop("ERROR! ALL UP-REGULATED CAN NOT BE ASSIGNED A BIOTYPE")}
  statsPerBioType<-c(nPcgUpRegulated=nPcgUpRegulated,nPcgDownRegulated=nPcgDownRegulated,nLncRNAUpRegulated=nLncRNAUpRegulated,nLncRNADownRegulated=nLncRNADownRegulated)
  return(statsPerBioType)
}
