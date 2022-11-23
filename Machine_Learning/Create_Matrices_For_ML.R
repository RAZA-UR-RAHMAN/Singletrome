###################################################### LOADS SEURAT OBJECTS AND CREATE MATRICES FOR CELL TYPE AND DISEASE PREDICTION ###################################################### 
library(Seurat)
library(rhdf5)
library(scrattch.io)
############ RUN FOR EACH DATASET ########################
outdir<-"MACHINE_LEARNING/"
######################################################### LNCRNA.RDS  ##############################################################################
seuratObjectPath="LNCRNA.RDS"
seuratObj<-readRDS(seuratObjectPath)
rawCounts=GetAssayData(object = seuratObj, slot = "counts")
rownames(rawCounts)[grep(x = rownames(rawCounts),pattern = ",")]<- gsub(',', '_', rownames(rawCounts)[grep(x = rownames(rawCounts),pattern = ",")] )
write_dgCMatrix_csv(rawCounts, paste0(outdir,"GSE136103-LncRNA-ScrattchMatrix.csv"), col1_name = "gene",chunk_size = 2000)
write.table(seuratObj@meta.data,paste0(outdir,'GSE136103-LncRNA-metadata.csv'), sep = ',', row.names = T, col.names = NA, quote = F)
######################################################### PROTEINCODING.RDS  ########################################################################
seuratObjectPath="PROTEINCODING.RDS"
seuratObj<-readRDS(seuratObjectPath)
rawCounts=GetAssayData(object = seuratObj, slot = "counts")
rownames(rawCounts)[grep(x = rownames(rawCounts),pattern = ",")]<- gsub(',', '_', rownames(rawCounts)[grep(x = rownames(rawCounts),pattern = ",")] )
write_dgCMatrix_csv(rawCounts, paste0(outdir,'GSE136103-Protein-ScrattchMatrix.csv'), col1_name = "gene",chunk_size = 2000)
write.table(seuratObj@meta.data,paste0(outdir,'GSE136103-Protein-metadata.csv'), sep = ',', row.names = T, col.names = NA, quote = F)
######################################################### SINGLETROME.RDS  ##########################################################################
seuratObjectPath="SINGLETROME.RDS"
seuratObj<-readRDS(seuratObjectPath)
rawCounts=GetAssayData(object = seuratObj, slot = "counts")
rownames(rawCounts)[grep(x = rownames(rawCounts),pattern = ",")]<- gsub(',', '_', rownames(rawCounts)[grep(x = rownames(rawCounts),pattern = ",")] )
write_dgCMatrix_csv(rawCounts, paste0(outdir,'GSE136103-SINGLETROME-ScrattchMatrix.csv'), col1_name = "gene",chunk_size = 2000)
write.table(seuratObj@meta.data,paste0(outdir,'GSE136103-SINGLETROME-metadata.csv'), sep = ',', row.names = T, col.names = NA, quote = F)
