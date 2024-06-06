################################################################################
###################CellChat Comparisons between all SEC stages##################
################################################################################

# Load the required libraries
library(CellChat)
library(patchwork)
library(Seurat)

#Set Working Directory
setwd("/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/scripts")

#### Data Preparation ####
#Load Data
Int=readRDS(file ="./data/DiestrusMice_mU7_mU30_Final_mergedCellID_02012024.rds",  # Filename
            refhook = NULL)
ncol(Int) #37,543 cells

Idents(Int) <- Int$seurat_clusters_cellType

Int <- subset(Int, idents='Meso', invert=TRUE)
ncol(Int) #36,913

#Removed Mesothelial Level
meta <- Int@meta.data
meta$seurat_clusters_cellType <- droplevels(meta$seurat_clusters_cellType, exclude='Meso')

#Split integrated object 
Int <- SplitObject(Int, split.by="RedSEC_stage")

#Prepare Data
Int.chat_cont <- createCellChat(as.matrix(Int$control), group.by='seurat_clusters_cellType', assay='RNA', do.sparse=T)  #normal
Int.chat_ear <- createCellChat(as.matrix(Int$earlySEC), group.by='seurat_clusters_cellType', assay='RNA', do.sparse=T)  #pre-dysplastic
Int.chat_late <- createCellChat(as.matrix(Int$lateSEC), group.by='seurat_clusters_cellType', assay='RNA', do.sparse=T)  #dysplastic

Int.chat_cont <- subsetCellChat(Int.chat_cont, idents.use='Meso', invert=TRUE)
Int.chat_ear <- subsetCellChat(Int.chat_ear, idents.use='Meso', invert=TRUE)
Int.chat_late <- subsetCellChat(Int.chat_late, idents.use='Meso', invert=TRUE)


#Set Ligand Receptor pair database - control
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB
Int.chat_cont@DB <- CellChatDB.use

Int.chat_cont <- subsetData(Int.chat_cont)

#Set Ligand Receptor pair database - early
CellChatDB.use <- CellChatDB
Int.chat_ear@DB <- CellChatDB.use

Int.chat_ear <- subsetData(Int.chat_ear)

#Set Ligand Receptor pair database - late
CellChatDB.use <- CellChatDB
Int.chat_late@DB <- CellChatDB.use

Int.chat_late <- subsetData(Int.chat_late)

#Compute Computational Probability - normal subgroup
Int.chat_cont <- identifyOverExpressedGenes(Int.chat_cont)
Int.chat_cont <- identifyOverExpressedInteractions(Int.chat_cont)
Int.chat_cont <- projectData(Int.chat_cont, PPI.mouse) #Projected interactions
Int.chat_cont <- computeCommunProb(Int.chat_cont, raw.use=TRUE, population.size=TRUE)
Int.chat_cont <- filterCommunication(Int.chat_cont, min.cells=10) #filter out communication if only few number of cells in group
Int.chat_cont <- computeCommunProbPathway(Int.chat_cont)
Int.chat_cont <- aggregateNet(Int.chat_cont)
groupsize_cont <- as.numeric(table(Int.chat_cont@idents))
par(mfrow=c(1,2), xpd=TRUE)

#Compute Computational Probability - pre-dysplastic subgroup
Int.chat_ear <- identifyOverExpressedGenes(Int.chat_ear)
Int.chat_ear <- identifyOverExpressedInteractions(Int.chat_ear)
Int.chat_ear <- projectData(Int.chat_ear, PPI.mouse) #Projected interactions
Int.chat_ear <- computeCommunProb(Int.chat_ear, raw.use=TRUE, population.size=TRUE)
Int.chat_ear <- filterCommunication(Int.chat_ear, min.cells=10) #filter out communication if only few number of cells in group
Int.chat_ear <- computeCommunProbPathway(Int.chat_ear)
Int.chat_ear <- aggregateNet(Int.chat_ear)
groupsize_ear <- as.numeric(table(Int.chat_ear@idents))
par(mfrow=c(1,2), xpd=TRUE)

#Compute Computational Probability - dysplastic subgroup
Int.chat_late <- identifyOverExpressedGenes(Int.chat_late)
Int.chat_late <- identifyOverExpressedInteractions(Int.chat_late)
Int.chat_late <- projectData(Int.chat_late, PPI.mouse) #Projected interactions
Int.chat_late <- computeCommunProb(Int.chat_late, raw.use=TRUE, population.size = TRUE)
Int.chat_late <- filterCommunication(Int.chat_late, min.cells=10) #filter out communication if only few number of cells in group
Int.chat_late <- computeCommunProbPathway(Int.chat_late)
Int.chat_late <- aggregateNet(Int.chat_late)
groupsize_late <- as.numeric(table(Int.chat_late@idents))
par(mfrow=c(1,2), xpd=TRUE)

#Save CellChat Objects
saveRDS(Int.chat_cont, './data/Control_CellChat_02202024_noMeso.rds')
saveRDS(Int.chat_ear, './data/Early_CellChat_02202024_noMeso.rds')
saveRDS(Int.chat_late, './data/Late_CellChat_02202024_noMeso.rds')