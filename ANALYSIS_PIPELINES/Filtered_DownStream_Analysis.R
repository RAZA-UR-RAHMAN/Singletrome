######################################################## Reperform Clustering after filtering lncRNAs ######################################################## 
#TO SPLIT PCG and LNCRNAs
## Gene lists (protein coding and lncRNA genes)
referenceResources="Genomes/reference_sources/"
proteinCodingGenes=read.table(paste0(referenceResources,"protein_coding_genenames_and_ids_unique.txt"), quote = "",row.names=NULL,sep="\t")[[1]]
singletromeLncRNAs=read.table(paste0(referenceResources,"Singletrome_genenames_and_ids_unique.txt"), quote = "",row.names=NULL,sep="\t")[[1]]
GENCODELncRNAs=read.table(paste0(referenceResources,"GENECODE_lncRNA_genenames_and_ids_unique.txt"), quote = "",row.names=NULL,sep="\t")[[1]]

# "LINC00595', 'LINC00856,LINC00856,LINC00595" added this manually. It is from lncExpDB
# add .1 .2 and .3. Cellranger adds these if gene name is duplicated for different gene ids
proteinCodingGenes=c(proteinCodingGenes,paste0(proteinCodingGenes,".1"),paste0(proteinCodingGenes,".2"),paste0(proteinCodingGenes,".3"),paste0(proteinCodingGenes,".4"))
singletromeLncRNAs=c(singletromeLncRNAs,paste0(singletromeLncRNAs,".1"),paste0(singletromeLncRNAs,".2"),paste0(singletromeLncRNAs,".3"),paste0(singletromeLncRNAs,".4"),"LINC00595', 'LINC00856,LINC00856,LINC00595")
GENCODELncRNAs=c(GENCODELncRNAs,paste0(GENCODELncRNAs,".1"),paste0(GENCODELncRNAs,".2"),paste0(GENCODELncRNAs,".3"))

####################################################################### VARIABLES
library("ggpubr")
source("CustomTheme.R")
source("Utility/FUNCTIONS.R")
filteredTranscriptsDir="Results/DOWNSTREAM/FILTERED_TRANSCRIPTS/FILTERED_TRANSCRIPTS_READ_DISTRIBUTION/"
resultsDir="Results/"
clusteringDir=paste0(resultsDir,"DOWNSTREAM/FILTERED_CLUSTERING/")
datasets<-c("pbmc_10k_v3","GSE115469","GSE136103")

pbmcsCols = c("CD8 Naive"='darkorange',  'CD4 Memory' = '#6EC5E2', 'CD14+ Monocytes' = 'lightgreen', 'NK cell' = '#DBAD33','pre-B cell'='#7D93D5','CD4 Naive'='#D48CD1','pDC'='darkblue',
              'Double negative T cell'='#D95393','CD16+ Monocytes'='#DBD160',
              'Platelets'='#9FEA4C','CD8 effector'='#883AE4','B cell progenitor'='#E56557','Dendritic cell'='red')
GSE115469Cols<-c("LSECs_11"='darkorange',  'LSECs_12' = '#6EC5E2', 'LSECs_13' = 'lightgreen', 'Cholangiocytes' = '#DBAD33','ab T cells'='#7D93D5','gd T cells_9'='#D48CD1','gd T cells_18'='lightblue',
                 'Macrophages_4'='#D95393','Macrophages_10'='#DBD160','NK cells'='#9FEA4C','Mature B cells'='#883AE4','Hepatic Stellate Cells'='#E56557','Plasma cells'='red','Erythroid cells'='#ACE5B7','Hepatocytes_1'='#E46753','Hepatocytes_3'='#71B859','Hepatocytes_5'='#6356D9','Hepatocytes_6'='#9A7989',
                 'Hepatocytes_14'='#DEC04C','Hepatocytes_15'='#A534EC')
GSE136103Cols = c("MPs"='darkorange',  'Tcells' = '#6EC5E2', 'ILCs' = 'lightgreen', 'Bcells' = '#DBAD33','pDCs'='#7D93D5','Plasma Bcells'='#D48CD1',
                  'Mesenchyme'='lightblue','Hepatocytes'='#D95393','Mast cells'='#DBD160','Endothelia'='#9FEA4C','Cholangiocytes'='#883AE4','Mesothelia'='red')
for(ds in datasets){
  dsOutDir<-paste0(clusteringDir,ds,"/")
  if(!dir.exists(dsOutDir)){dir.create(dsOutDir,mode="0777",recursive = T)}
  cols = switch(  
    ds,  
    "pbmc_10k_v3"= pbmcsCols,
    "GSE115469"= GSE115469Cols,
    "GSE136103"= GSE136103Cols,
  )  
  # load gencode
  GENCODEObj<-readRDS(paste0(resultsDir,"DOWNSTREAM/CLUSTERING/",ds,"/GENCODE/COMBINED.RDS"))
  GENCODELncRNAObj<-readRDS(paste0(resultsDir,"DOWNSTREAM/CLUSTERING/",ds,"/GENCODE/LNCRNA.RDS"))
  
  # load singletrome
  filteredSingletromeObj<-readRDS(paste0(filteredTranscriptsDir,ds,"/",ds,"_QC_Filtered_Singletrome.RDS"))
  # choose lncRNA list to split
  proteinSeuratObj <- subset(filteredSingletromeObj, features = proteinCodingGenes)
  lncRNASeuratObj <- subset(filteredSingletromeObj, features = singletromeLncRNAs)
  if(nrow(filteredSingletromeObj)!=(nrow(proteinSeuratObj)+nrow(lncRNASeuratObj))){stop("nrow(filteredSingletromeObj)!=(nrow(proteinSeuratObj)+nrow(lncRNASeuratObj))")}
  
  filteredSingletromeObj=PerformDownStreamAnalysis(filteredSingletromeObj,dsOutDir,"COMBINED")
  proteinSeuratObj=PerformDownStreamAnalysis(proteinSeuratObj,dsOutDir,"PROTEINCODING")
  lncRNASeuratObj=PerformDownStreamAnalysis(lncRNASeuratObj,dsOutDir,"LNCRNA")
 
  
  
  GENCODEUmap<-DimPlot(GENCODEObj, reduction = "umap", group.by='celltype',label = TRUE, label.size = 4.3 , cols=cols) +  ggtitle("GENCODE") +umap_theme # no need to recalculate this umap, because we are not chaning anything
  singletromeUmap<-DimPlot(filteredSingletromeObj, reduction = "umap", group.by='celltype',label = TRUE, label.size = 4.3,cols=cols) +  ggtitle("Singletrome") +umap_theme
  pcgUmap<-DimPlot(proteinSeuratObj, reduction = "umap", group.by='celltype',label = TRUE,label.size = 4.3,cols=cols) +  ggtitle("Protein-coding") +umap_theme
  lncRNAUmap<-DimPlot(lncRNASeuratObj, reduction = "umap", group.by='celltype',label = TRUE, label.size = 4.3,cols=cols) +  ggtitle("lncRNA") +umap_theme 
  # plot
  AllUmaps <- ggarrange(GENCODEUmap,singletromeUmap,pcgUmap,lncRNAUmap, labels = c("A.", "B.","C.","D."), ncol = 2, nrow = 2,common.legend = TRUE, align="hv",font.label=list(size=20),legend = "bottom")
  ggexport(AllUmaps, filename = paste0(dsOutDir,ds,"_umap.pdf"),width = 8, height = 8)
  
  ## EXPRESSION
  # Merge by gene type (protein coding and lncRNA)
  proteinSeuratObj@meta.data$type="SingletromeProteinCoding"
  lncRNASeuratObj@meta.data$type="SingletromeLncRNAs"
  GENCODELncRNAObj@meta.data$type="GENCODELncRNAs"
  combined <- merge(proteinSeuratObj, y = c(lncRNASeuratObj,GENCODELncRNAObj), add.cell.ids = c("SingletromeProteinCoding", "SingletromeLncRNAs","GENCODELncRNAs"), project = "ProteinAndLncRNAsCombined")
  combined@meta.data$type <- factor(x = combined@meta.data$type, levels = c('SingletromeProteinCoding', 'SingletromeLncRNAs','GENCODELncRNAs')) # Reorder for plots
  # Number of reads
  nCount_RNA_Combined<-VlnPlot(combined, features = c("nCount_RNA"), group.by = "type",log=T,pt.size = 0)+geom_boxplot(width=0.1,fill="white") + ggtitle("Reads") +theme_Publication() + scale_fill_Publication() + ylab("log (Number of reads)") +xlab("") + theme(axis.text.x = element_blank(),axis.text.y = element_text(size = rel(1)),axis.title.y=element_text(size = 10), plot.title = element_text(size = 10, face = "bold"),legend.key.size = unit(0.5, 'cm'),legend.text = element_text(size=10))
  nFeature_RNA_Combined<-VlnPlot(combined, features = c("nFeature_RNA"), group.by = "type",log=T,pt.size = 0)+geom_boxplot(width=0.1,fill="white") + ggtitle("Genes") +theme_Publication() + scale_fill_Publication() + ylab("log (Number of genes)") +xlab("") + theme(axis.text.x = element_blank(),axis.text.y = element_text(size = rel(1)),axis.title.y=element_text(size = 10), plot.title = element_text(size = 10, face = "bold"),legend.key.size = unit(0.5, 'cm'),legend.text = element_text(size=10))
  
  # plot
  nCountRNAFigure <- ggarrange(nCount_RNA_Combined,nFeature_RNA_Combined, labels = c("E.", "F."), ncol = 2, nrow = 1, font.label=list(size=20),common.legend = TRUE)
  #ggexport(nCountRNAFigure,   filename = paste0(dsOutDir,ds,"_nCount_RNA-nFeature_RNA.pdf"),width = 5, height = 3) #. Main Figure Dimension
  ggexport(nCountRNAFigure,   filename = paste0(dsOutDir,ds,"_nCount_RNA-nFeature_RNA.pdf"),width = 5, height = 4) #. Supplementary Figure Dimension
  
  #### Cell type markers
  N=3
  proteinSeuratObj<-findAllMarkersForCellType(proteinSeuratObj,N,dsOutDir,"protein_coding")
  lncRNASeuratObj<-findAllMarkersForCellType(lncRNASeuratObj,N,dsOutDir,"lncRNA")
}
print("-- CLUSTERING AND UMAP DONE---")
#################################################################  Differential Expression ################################################################# 
library("ggpubr")
source("CustomTheme.R")
source("Analysis/Scripts/Utility/FUNCTIONS.R")
resultsDir="Results/"
clusteringDir=paste0(resultsDir,"DOWNSTREAM/FILTERED_CLUSTERING/")
differentialExpressionDir=paste0(resultsDir,"DOWNSTREAM/FILTERED_DIFFERENTIAL_EXPRESSION/")

analysisConfigPath="Analysis/Configs/analysisConfig.txt" # to load summary files
analysisConfigTable=read.table(analysisConfigPath,sep="\t",header = T)
analysisConfigTableSubset=analysisConfigTable[which(!is.na(analysisConfigTable$decolumn)),]
geneSets<-c("LNCRNA","COMBINED","PROTEINCODING")
for(deRow in 1:nrow(analysisConfigTableSubset)){
  performDEOnColumn=as.character(analysisConfigTableSubset[deRow,"decolumn"])
  datasetID=as.character(analysisConfigTableSubset[deRow,"baseDir"])
  seuratObjectsDir=paste0(clusteringDir,"/",datasetID,"/")
  for(gs in geneSets){
    print(paste0("Running for ", datasetID," and for gene set ",gs))
    seuratObj=readRDS(paste0(seuratObjectsDir,gs,".RDS"))
    # outdir creation
    differentialExpressionDIR=paste0(differentialExpressionDir,datasetID,"/",gs,"/")
    if(!dir.exists(differentialExpressionDIR)){dir.create(differentialExpressionDIR,mode="0777",recursive = T)}
    # perform DE
    #Note that the raw and normalized counts are stored in the counts and data slots of RNA assay. By default, the functions for finding markers will use normalized data.
    DefaultAssay(seuratObj) <- "RNA"
    seuratObj<- ScaleData(object = seuratObj, features = rownames(seuratObj)) # some genes are missing in the scaled version (when plotting the heatmap) so will rescale (https://github.com/satijalab/seurat/issues/1369)
    seuratObj<-findAllMarkersForCellType(seuratObj,3,differentialExpressionDIR,gs)
    findDEGInCelltypesAcrossConditions(seuratObj,performDEOnColumn,differentialExpressionDIR,proteinCodingGenes,singletromeLncRNAs)
    }
}

print("------------------DE COMPLETED------------------")
