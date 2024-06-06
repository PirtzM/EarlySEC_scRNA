################################################################################
##################Analysis of Final Normal Samples (5 total)####################
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
IntData = readRDS(file ='./data/DiestrusMice_mU7_mU30_Final_01252024.rds',  # Filename
                  refhook = NULL)
ncol(IntData) #37543 cells

#Subset by Cancer Status and Recluster####
Idents(IntData) <- IntData$Cancer_Status #Set Active Identity
ContOnly <- subset(IntData, idents=c('control')) #Normal samples
ncol(ContOnly) #7,614 cells

#Rescale and Recluster
ContOnly <- FindVariableFeatures(object = ContOnly, assay = 'RNA', 
                                 selection.method = 'vst', nfeatures = 2000)
ContOnly = ScaleData(ContOnly, vars.to.regress = c("nFeature_RNA", "percent.mt"),
                     verbose = FALSE)
ContOnly = RunPCA(ContOnly, verbose = FALSE, npcs = 50)

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

n.pcs=npcs(seu = ContOnly, var.toal = 0.95, reduction = 'pca_RNA') # Use this value in the dimensions below
n.pcs
#PCs=39

# Run Harmony
gc()
ContOnly <- RunHarmony(ContOnly,                     # Object with PCA computed
                       group.by.vars = "Sample.ID",    # Which variables to remove
                       reduction = "pca_RNA",          # Name of reduction to use
                       assay = "RNA",                  # Name of assay to use.  Default is RNA assay
                       plot_convergence = TRUE)        # Whether to plot convergence


#Clustering
ContOnly <- FindNeighbors(object = ContOnly,           # Seurat object
                          reduction = "harmony",         # Name of reduction to use
                          k.param =30,                  # Defines k for the k-nearest neighbor algorithm (smaller the number the more clusters)
                          dims = 1:n.pcs,                   # Dimensions of reduction to use (use the result from npcs())
                          graph.name = 'ContMice.DF.Clusters')    # Name for stored SNN graph

ContOnly <- FindClusters(object = ContOnly,            # Seurat object 
                         resolution = 0.7,
                         graph.name = 'ContMice.DF.Clusters') # Resolution parameter (smaller the number the fewer clusters) 

# Set Python UMAP via reticulate
umap.method = 'umap-learn'
metric = 'correlation'

ContOnly <- RunUMAP(object = ContOnly,      # Seurat object 
                    reduction = 'harmony',    # Reduction to use
                    dims = 1:n.pcs)              # Dimensions of recution to use (use the result from npcs())

ContOnlyUMAP <- DimPlot(object = ContOnly,                 # Seurat object 
                        reduction = 'umap',
                        group.by = "seurat_clusters",     # Labels to color the cells by ("seurat_clusters", "Age", "Time.Point)  
                        repel = TRUE,                       # Whether to repel the cluster labels
                        label = TRUE,                       # Whether to have cluster labels 
                        pt.size = 1,                      # Size of each dot is (0.1 is the smallest)
                        label.size = 5)                    # Font size for label
ContOnlyUMAP

# Save as RDS files
saveRDS(ContOnly, file = "./data/cont_mU7_mU30_Final_01262024.rds")

#Load Status RDS Files ####
ContOnly = readRDS(file ="./data/cont_mU7_mU30_Final_01262024.rds",  # Filename
                   refhook = NULL)
ncol(ContOnly) #7,614 cells

Idents(ContOnly) <- ContOnly$seurat_clusters2 #Set Active Identity

#Summary of Markers ####
Sumfeatures = c('Epcam','Pax8','TdTomato-UTR', #General Epithelium
                'Tacstd2','Prap1', 'Krt13', 'Klf5','Lpar3','Ly6a','Met','Sprr2f','Slc1a3','Wnt7a','Wnt7b','Ptgs1', #Luminal Epithelium
                'Foxa2', 'Aldh1a1','Axin2','Lgr5','Ltf','Sox17','Lgr4','Muc1','Ovgp1','Wif1','Slc18a2', #Glandular Epithelium
                'Cdkn2a','Ccne1','Cdk1','Cdk4','Cdkn1c','Pcna','Mki67', #Cycling genes
                "Adamts1",'Col1a1','Col1a2','Col3a1','Creb5','Dcn','Dio2','Gfpt2','Igf1','Lrrn4','Medag','Mgp','Msln','Pdgfra','Rxfp1','Vcan','Vim', #Fibroblasts
                'Acta1','Des',"Myh11",'Acta2','Palld', # Muscle
                'Lyve1','Mmrn1', #lymphatic endothelium
                'Cavin2','Cldn5','Cdh5','Flt4','Icam1','Pecam1','Tek','Vwf' #Vascular Endothelium
                
)
SumDot = DotPlot(ContOnly, features=Sumfeatures, col.min=0, col.max=2.5) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic'))
SumDot

#Stem Cell Markers ####
SCfeatures = c('Aldh1a1','Ly6a','Cd44','Ezh2','Fut4',
               'Itga6','Klf4','Klf5','Klf6','Lgr5','Lrig1','Neat1','Pax8','Prom1', 'Pou5f1',
               'Slc1a3','Slc38a2','Sox2','Sox9','Sox17','Tacstd2','Zfp36')
DotSC = DotPlot(ContOnly, features=SCfeatures, col.min=0, col.max=2.5) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic'))
DotSC

#Diff Gene Expression - Progenitor Populations
DotPlot(ContOnly, features=c('Pax2','Hoxb5','Tert','Sox17','Klf5', #nucleus
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
ImDot = DotPlot(ContOnly, features=Imfeatures, col.min=0, col.max=2.5) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic')) #angle x axis labels
ImDot

#Cluster Renaming ####
#Fine Grain
Idents(ContOnly) <- ContOnly$seurat_clusters

CellTypeNames <- RenameIdents(ContOnly, 
                              '0'= 'LE 1',
                              '1'= 'Lymph',
                              '2'= 'Progenitor-like',
                              '3'= 'Cycling',
                              '4'= 'Fib 1',
                              '5'= 'LE 2',
                              '6'= 'Vascular',
                              '7'= 'GE',
                              '8'= 'Fib 2',
                              '9'= 'Unidentified', #referred to as 'mt High' in manuscript
                              '10'= 'Macrophage',
                              '11'= 'NK cell',
                              '12'= 'Smooth Muscle',
                              '13'= 'T cell',
                              '14'= 'Mesothelium'
)

ContOnly$seurat_clusters2 <- Idents(CellTypeNames)

Idents(ContOnly) <- ContOnly$seurat_clusters2

#Course Grain
Idents(ContOnly) <- ContOnly$seurat_clusters2

CellTypeNames <- RenameIdents(ContOnly, 
                              'Progenitor-like'='EP',
                              'LE 1'='LE',
                              'LE 2'='LE',
                              'Cycling'='CE',
                              'GE'='GE',
                              'Unidentified'='UE', #referred to as 'mt High' in manuscript
                              'Fib 1'='Fib',
                              'Fib 2'='Fib',
                              'Mesothelium'='Meso',
                              'Smooth Muscle'='SM',
                              'Lymph'='LV',
                              'Vascular'='BV',
                              'Macrophage'='Im',
                              'NK cell'='Im',
                              'T cell'='Im'
                              
)

ContOnly$seurat_clusters_CC <- Idents(CellTypeNames)
saveRDS(ContOnly, "./data/cont_mU7_mU30_Final_01262024.rds")

#UMAP ####
UMAP = DimPlot(object=ContOnly, reduction="umap", group.by = "seurat_clusters2",
               repel = FALSE,                    
               label = FALSE,
               pt.size = 0.05,                     
               label.size = 4
)
UMAP