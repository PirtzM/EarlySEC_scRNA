################################################################################
################Analysis of Pre-dysplastic Samples (7 total)####################
################################################################################

#R v4.1.1
#Load Libraries
library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)

#Set working directory
setwd("/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/scripts")

#Load full object ####
IntData = readRDS(file ='./data/DiestrusMice_mU7_mU30_Final_01182025_simple.rds',  # Filename
                  refhook = NULL)
ncol(IntData) #37543 cells

#Subset by SEC stage and Recluster####
Idents(IntData) <- IntData$RedSEC_stage #Set Active Identity
EarOnly <- subset(IntData, idents=c('Pre-dysplastic')) #Pre-dysplastic samples
ncol(EarOnly) #17,212 cells

#Rescale and Recluster
EarOnly <- FindVariableFeatures(object = EarOnly, assay = 'RNA', 
                                selection.method = 'vst', nfeatures = 2000)
EarOnly = ScaleData(EarOnly,
                    vars.to.regress = c("nFeature_RNA", "percent.mt"),
                    verbose = F)
EarOnly = RunPCA(EarOnly, verbose = F, npcs = 50)

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

n.pcs = npcs(seu = EarOnly, var.toal = 0.95, reduction = 'pca_RNA') # Use this value in the dimensions below
n.pcs
#PCs=39

# Run Harmony
gc()
EarOnly <- RunHarmony(EarOnly,                     # Object with PCA computed
                      group.by.vars = "Sample.ID",    # Which variables to remove
                      reduction = "pca_RNA",          # Name of reduction to use
                      assay = "RNA",                  # Name of assay to use.  Default is RNA assay
                      plot_convergence = TRUE)        # Whether to plot convergence

#Clustering - mU7 with Doublet finder
EarOnly <- FindNeighbors(object = EarOnly,           # Seurat object
                         reduction = "harmony",         # Name of reduction to use
                         k.param = 50,                  # Defines k for the k-nearest neighbor algorithm (smaller the number the more clusters)
                         dims = 1:n.pcs,                   # Dimensions of reduction to use (use the result from npcs())
                         graph.name = 'EarMice.DF.Clusters')    # Name for stored SNN graph

EarOnly <- FindClusters(object = EarOnly,            # Seurat object 
                        resolution = 0.7,
                        graph.name = 'EarMice.DF.Clusters') # Resolution parameter (smaller the number the fewer clusters) 

# Set Python UMAP via reticulate
umap.method = 'umap-learn'
metric = 'correlation'

EarOnly<- RunUMAP(object = EarOnly,      # Seurat object 
                  reduction = 'harmony',    # Reduction to use
                  dims = 1:n.pcs)              # Dimensions of recution to use (use the result from npcs())

EarOnlyUMAP <- DimPlot(object = EarOnly,                 # Seurat object 
                       reduction = 'umap',
                       group.by = "seurat_clusters",     # Labels to color the cells by ("seurat_clusters", "Age", "Time.Point)  
                       repel = TRUE,                       # Whether to repel the cluster labels
                       label = TRUE,                       # Whether to have cluster labels 
                       pt.size = 1,                      # Size of each dot is (0.1 is the smallest)
                       label.size = 5)                    # Font size for label
EarOnlyUMAP
FeaturePlot(EarOnly, features='Epcam')

# Save as RDS files
saveRDS(EarOnly, file = "./data/EarCan_recluster_mU7_mU30_01182025_simple.rds")

#Sample Analysis####
#Load Status RDS Files ####
EarOnly = readRDS(file ="./data/EarCan_recluster_mU7_mU30_01182025_simple.rds",  # Filename
                  refhook = NULL)
ncol(EarOnly) #17,212 cells

Idents(EarOnly) <- EarOnly$seurat_clusters2 #Set Active Identity

#Summary of Markers ####
Sumfeatures = c('Epcam','Pax8','TdTomato-UTR', #General Epithelium
                'Tacstd2','Krt13', 'Klf5','Lpar3','Ly6a','Met','Sprr2f','Slc1a3','Wnt7a','Wnt7b','Ptgs1', #Luminal Epithelium
                'Foxa2', 'Aldh1a1','Axin2','Lgr5','Ltf','Prss29','Sox17','Lgr4','Muc1','Ovgp1','Wif1','Slc18a2', #Glandular Epithelium
                'Cdkn2a','Ccne1','Cdk1','Cdk4','Cdkn1c','Pcna','Mki67', #Cycling genes
                "Adamts1",'Col1a1','Col1a2','Col3a1','Creb5','Dcn','Dio2','Gfpt2','Igf1','Lrrn4','Medag','Mgp','Msln','Pdgfra','Rxfp1','Vcan','Vim', #Fibroblasts
                'Acta1','Des',"Myh11",'Acta2','Palld', # Muscle
                'Lyve1','Mmrn1', #lymphatic endothelium
                'Cavin2','Cldn5','Cdh5','Flt4','Icam1','Pecam1','Tek','Vwf' #Vascular Endothelium
                
)
SumDot = DotPlot(EarOnly, features=Sumfeatures, col.min=0, col.max=2.5) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic'))
SumDot

#Stem Cell Markers ####
SCfeatures = c('Aldh1a1','Ly6a','Cd44','Ezh2','Fut4',
               'Itga6','Klf4','Klf5','Klf6','Lgr5','Lrig1','Neat1','Pax8','Prom1', 'Pou5f1',
               'Slc1a3','Slc38a2','Sox2','Sox9','Sox17','Tacstd2','Zfp36')
DotSC = DotPlot(EarOnly, features=SCfeatures, col.min=0, col.max=2.5) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic'))
DotSC

#Diff Gene Expression - Progenitor Populations
DotPlot(EarOnly, features=c('Pax2','Hoxb5','Tert','Sox17','Klf5', #nucleus
                             'Nt5e','Lpar3','Lgr4', #membrane
                             'Wfdc2' #Secreted
), col.min=0, col.max=2.5
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic'))

#Immune Summaries ####
Imfeatures = c( 'Adgre1','Aif1','C1qa','Cd74','Csf1r','Fcgr1','Mrc1', #macrophage
                'Cd4','Cd3e','Cxcr6','Pdcd1','Ptprc', #T Cells
                'Ccl5','Itgam', 'Klrd1','Txk', 'Xcl1', #NK Cells
                'Cd24a','Cd38','Cd22', #B Cells
                'Cd14', 'Ccr2', 'Fcgr3' #monocytes
)
ImDot = DotPlot(EarOnly, features=Imfeatures, col.min=0, col.max=2.5) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic')) #angle x axis labels
ImDot

Idents(EarOnly) <- EarOnly$seurat_clusters

CellTypeNames <- RenameIdents(EarOnly, 
                              '0'= 'PD - LE',
                              '1'= 'PD - Fib 1',
                              '2'= 'Lymph',
                              '3'= 'COX/MAL Epi',
                              '4'= 'PD - Fib 3',
                              '5'= 'Cycling Epi',
                              '6'= 'Macrophage',
                              '7'= 'PD - Fib 2',
                              '8'= 'Vascular',
                              '9'= 'Lymphocyte',
                              '10'= 'Mesothelium',
                              '11'= 'GE',
                              '12'= 'Smooth Muscle',
                              '13'= 'PD - LE',
                              '14'= 'Vascular'
                              
)

EarOnly$seurat_clusters2 <- Idents(CellTypeNames)

Idents(EarOnly) <- EarOnly$seurat_clusters2

#Fine Grain - Simple Epithelium (SE)
Idents(EarOnly) <- EarOnly$seurat_clusters

CellTypeNames <- RenameIdents(EarOnly, 
                              '0'= 'Epithelium',
                              '1'= 'PD - Fib 1',
                              '2'= 'Lymph',
                              '3'= 'Epithelium',
                              '4'= 'PD - Fib 3',
                              '5'= 'Epithelium',
                              '6'= 'Macrophage',
                              '7'= 'PD - Fib 2',
                              '8'= 'Vascular',
                              '9'= 'Lymphocyte',
                              '10'= 'Mesothelium',
                              '11'= 'Epithelium',
                              '12'= 'Smooth Muscle',
                              '13'= 'Epithelium',
                              '14'= 'Vascular'
                              
)

EarOnly$seurat_clusters_SE <- Idents(CellTypeNames)

Idents(EarOnly) <- EarOnly$seurat_clusters_SE

#Fine Grain - for bar chart (BC)
Idents(EarOnly) <- EarOnly$seurat_clusters

CellTypeNames <- RenameIdents(EarOnly, 
                              '0'= 'Epithelium',
                              '1'= 'Fib 1',
                              '2'= 'Lymph',
                              '3'= 'Epithelium',
                              '4'= 'Fib 3',
                              '5'= 'Epithelium',
                              '6'= 'Macrophage',
                              '7'= 'Fib 2',
                              '8'= 'Vascular',
                              '9'= 'Lymphocyte',
                              '10'= 'Mesothelium',
                              '11'= 'Epithelium',
                              '12'= 'Smooth Muscle',
                              '13'= 'Epithelium',
                              '14'= 'Vascular'
                              
)

EarOnly$seurat_clusters_BC <- Idents(CellTypeNames)

Idents(EarOnly) <- EarOnly$seurat_clusters_BC

#Course Grain
Idents(EarOnly) <- EarOnly$seurat_clusters2

CellTypeNames <- RenameIdents(EarOnly, 
                              'COX/MAL Epi'='COX/MAL',
                              'PD - LE'='LE',
                              'Cycling Epi'='CE',
                              'GE'='GE',
                              'PD - Fib 1'='Fib',
                              'PD - Fib 2'='Fib',
                              'PD - Fib 3'='Fib',
                              'Mesothelium'='Meso',
                              'Lymph'='LV',
                              'Vascular'='BV',
                              'Macrophage'='Im',
                              'Lymphocyte'='Im',
                              'Smooth Muscle'='SM'
                              
)

EarOnly$seurat_clusters_CC <- Idents(CellTypeNames)


#Name by broad cell type
Idents(EarOnly) <- EarOnly$seurat_clusters2

CellTypeNames <- RenameIdents(EarOnly, 
                              'COX/MAL Epi'='Epi',
                              'PD - LE'='Epi',
                              'Cycling Epi'='Epi',
                              'GE'='Epi',
                              'PD - Fib 1'='Fib',
                              'PD - Fib 2'='Fib',
                              'PD - Fib 3'='Fib',
                              'Mesothelium'='Meso',
                              'Lymph'='LV',
                              'Vascular'='BV',
                              'Macrophage'='Im',
                              'Lymphocyte'='Im',
                              'Smooth Muscle'='SM'
                              
)

EarOnly$seurat_clusters_cellType <- Idents(CellTypeNames)
saveRDS(EarOnly, "./data/EarCan_recluster_mU7_mU30_01182025_simple.rds")

#UMAP ####
UMAP = DimPlot(object=EarOnly, reduction="umap", group.by='seurat_clusters2',
               repel = FALSE,                    
               label = FALSE,
               pt.size = 0.05,                    
               label.size = 4
)
UMAP
