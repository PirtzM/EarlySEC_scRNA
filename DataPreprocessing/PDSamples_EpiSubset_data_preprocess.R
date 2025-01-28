################################################################################
#################Analysis of Pre-dysplastic Samples (7 total)###################
##############################Epi Subset########################################

#R v4.1.1
#Load Libraries
library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)
library(harmony)
library(tidyr)

#Set working directory
setwd("/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/scripts")

#Load Pre-dysplastic Subset####
EpiOnly <- readRDS(file = "./data/EarCan_recluster_mU7_mU30_01182025_simple.rds",  # Filename
                   refhook = NULL)
ncol(EpiOnly) #17,212 cells

Idents(EpiOnly) <- EpiOnly$seurat_clusters_CC #Set Active Identity
EpiOnly <- subset(EpiOnly, idents=c('LE','GE','COX/MAL','CE')) #epithelial clusters only
ncol(EpiOnly) #7,713 cells

#Rescale and Recluster
#EpiOnly <- NormalizeData(object = EpiOnly, assay = 'RNA')
EpiOnly <- FindVariableFeatures(object = EpiOnly, assay = 'RNA', 
                                selection.method = 'vst', nfeatures = 2000)
EpiOnly = ScaleData(EpiOnly, vars.to.regress = c("nFeature_RNA", "percent.mt"),
                    verbose = T)
EpiOnly = RunPCA(EpiOnly, verbose = F, npcs = 50)

# Calculate the number of PCs that contain some proportion (95%) of the variance
# This is a custom function.  You need to run this before using the function.  Once it is in your environment you do not need to run it again.
npcs <- function(seu, var.toal=0.95, reduction="pca_RNA"){
  if(is.null(seu@reductions[[reduction]])){
    cat("Reduction", reduction, "not found!")
    return(NULL)
  }
  tmp.var <- (seu@reductions[[reduction]]@stdev)^2
  var.cut <- var.toal*sum(tmp.var)
  n.pcs=0
  var.sum = 0
  while(var.sum < var.cut){
    n.pcs = n.pcs + 1
    var.sum <- var.sum + tmp.var[n.pcs]
  }
  return(n.pcs)
}

n.pcs=npcs(seu = EpiOnly, var.toal = 0.95, reduction = 'pca_RNA') # Use this value in the dimensions below
n.pcs
#PCs=39

# Run Harmony
gc()
EpiOnly <- RunHarmony(EpiOnly,                     # Object with PCA computed
                      group.by.vars = "Sample.ID",    # Which variables to remove
                      reduction = "pca_RNA",          # Name of reduction to use
                      assay = "RNA",                  # Name of assay to use.  Default is RNA assay
                      plot_convergence = TRUE)        # Whether to plot convergence

#Clustering
EpiOnly <- FindNeighbors(object = EpiOnly,           # Seurat object
                         reduction = "harmony",         # Name of reduction to use
                         k.param = 70,                  # Defines k for the k-nearest neighbor algorithm (smaller the number the more clusters)
                         dims = 1:n.pcs,                   # Dimensions of reduction to use (use the result from npcs())
                         graph.name = 'EarMiceEpi.DF.Clusters')    # Name for stored SNN graph

EpiOnly <- FindClusters(object = EpiOnly,            # Seurat object 
                        resolution = 0.7,
                        graph.name = 'EarMiceEpi.DF.Clusters') # Resolution parameter (smaller the number the fewer clusters) 

# Set Python UMAP via reticulate
umap.method = 'umap-learn'
metric = 'correlation'

EpiOnly <- RunUMAP(object = EpiOnly,      # Seurat object 
                   reduction = 'harmony',    # Reduction to use
                   dims = 1:n.pcs)              # Dimensions of recution to use (use the result from npcs())

EpiOnlyUMAP <- DimPlot(object = EpiOnly,                 # Seurat object 
                       reduction = 'umap',
                       group.by = "seurat_clusters",     # Labels to color the cells by ("seurat_clusters", "Age", "Time.Point)  
                       repel = TRUE,                       # Whether to repel the cluster labels
                       label = TRUE,                       # Whether to have cluster labels 
                       pt.size = 1,                      # Size of each dot is (0.1 is the smallest)
                       label.size = 5)                    # Font size for label
EpiOnlyUMAP

# Save as RDS files
saveRDS(EpiOnly, file = "./data/Ear_Epi_mU7_mU30_recluster_01182025.rds")

#Sample Analysis####
#Load Reclustered dataset ####
EpiOnly = readRDS(file ="./data/Ear_Epi_mU7_mU30_recluster_01182025.rds",  # Filename
                  refhook = NULL)
ncol(EpiOnly) #7713 cells

Idents(EpiOnly) <- EpiOnly$seurat_clusters2 #Set Active Identity

#Summary Epi Features ####
Epifeatures = c('Epcam','Pax8','TdTomato-UTR', #General Epithelium
                'Ly6a','Cd44','Ezh2','Fut4',
                'Itga6','Klf4','Klf5','Klf6','Lgr5','Lrig1','Neat1','Pou5f1',
                'Slc1a3','Slc38a2','Sox2','Sox9','Sox17','Zfp36', #Progenitor 
                'Foxa2','Foxj1','Aldh1a1','Axin2','Prss29','Lgr4','Wif1','Prom1', #Glandular Epithelium
                'Tacstd2','Krt13','Lpar3','Met','Sprr2f','Wnt7a','Wnt7b','Muc1','Ovgp1',#Luminal Epithelium
                'Cdkn2a','Ccne1','Cdk1','Cdk4','Cdkn1c','Pcna','Mki67' #Cycling genes
)
EpiDot = DotPlot(EpiOnly, features=Epifeatures, col.min=0, col.max=2.5)+
  scale_color_viridis(option = 'A', begin=0.9, end=0)+ #change to magma color scale
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic'))
EpiDot

#Stem Markers ####
SCfeatures = c('Aldh1a1','Ly6a','Cd44','Ezh2','Fut4','Hoxb5','Tert','Nt5e',
               'Itga6','Klf4','Klf5','Klf6','Lgr5','Lrig1','Neat1','Pax8','Prom1', 'Pou5f1',
               'Slc1a3','Slc38a2','Sox2','Sox9','Sox17','Tacstd2','Zfp36')
DotSC = DotPlot(EpiOnly, features=SCfeatures, col.min=0, col.max=2.5) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic'))
DotSC

#Cluster Renaming ####
#Fine Clusters
Idents(EarEpi) <- EarEpi$seurat_clusters

CellTypeNames <- RenameIdents(EarEpi, 
                              '0'= 'PD - LE 2', 
                              '1'= 'PD - LE 1', 
                              '2'= 'COX/MAL', 
                              '3'= 'Cycling',
                              '4'= 'DDP', 
                              '5'= 'PD - LE 3', 
                              '6'= 'GE',
                              '7'= 'PD - LE 4'
)

EarEpi$seurat_clusters2 <- Idents(CellTypeNames)

Idents(EarEpi) <- EarEpi$seurat_clusters2

#Fine Clusters - Simplified
Idents(EarEpi) <- EarEpi$seurat_clusters

CellTypeNames <- RenameIdents(EarEpi, 
                              '0'= 'LE 2',
                              '1'= 'LE 1',
                              '2'= 'COX/MAL', 
                              '3'= 'Cycling',
                              '4'= 'DDP',
                              '5'= 'LE 3',
                              '6'= 'GE',
                              '7'= 'LE 4'
)

EarEpi$seurat_clusters_SE <- Idents(CellTypeNames)

Idents(EarEpi) <- EarEpi$seurat_clusters_SE

#Coarse Clusters
Idents(EarEpi) <- EarEpi$seurat_clusters2

CellTypeNames <- RenameIdents(EarEpi, 
                              'PD - LE 1'= 'LE',
                              'PD - LE 2'= 'LE',
                              'COX/MAL'= 'CM', 
                              'Cycling'= 'CE',
                              'DDP'= 'DDP',
                              'PD - LE 3'= 'LE',
                              'GE'= 'GE',
                              'PD - LE 4'= 'LE'
)

EarEpi$seurat_clusters_CC <- Idents(CellTypeNames)

Idents(EarEpi) <- EarEpi$seurat_clusters_CC

saveRDS(EarEpi, "./data/Ear_Epi_mU7_mU30_recluster_01182025.rds")

#UMAP ####
UMAP = DimPlot(object=EpiOnly, reduction="umap", group.by = "seurat_clusters2",
               repel = FALSE,                     
               label =FALSE,
               pt.size = 0.75,                     
               label.size = 4
)
UMAP
