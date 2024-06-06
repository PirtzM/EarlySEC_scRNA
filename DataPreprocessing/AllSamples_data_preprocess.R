################################################################################
############Preprocessing of All Diestrus Samples at All Stages (20 total)######
################################################################################
##Workflow adapted from David McKellar (github: mckellardw) and Lauren Walter###

#R v4.1.1
#Load Libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(cowplot)
library(ggplot2)
library(SoupX)
library(tidyr)
library(DoubletFinder)

#Set working directory
setwd("/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/scripts")
#Load CSV with metadata
meta <- read.csv("/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/SC_diestrus_metadirectory_final.csv") #CSV with metadata and directory the "outs" are contained within
View(meta)


#Make SoupX Directory ####
soup.list <- lapply( 
  as.list(paste0(meta$data.dir,"/outs")),
  FUN = SoupX::load10X,
  keepDroplets=TRUE
)

soup.list.est <- lapply(
  soup.list,
  FUN = function(sc){
    return(tryCatch(autoEstCont(sc), error=function(e) NULL))
  }
)

#Estimated rho of 0.01 for all samples
adj.mat.list <- lapply(
  soup.list.est,
  FUN = function(sc){
    return(tryCatch(adjustCounts(sc), error=function(e) NULL))
  }
)

#Save SoupX adjust matrices to disk!
for(i in 1:length(adj.mat.list)){
  if(!is.null(adj.mat.list[[i]]) & !file.exists(paste0(meta$data.dir[i],"/outs/soupx/matrix.mtx"))){
    DropletUtils:::write10xCounts(
      path=paste0(meta$data.dir[i],"/outs/soupx"), #path to each sample's cellranger count output
      adj.mat.list[[i]]
    )
  }
}

# save Rho values
rhos <- list()
for(i in 1:length(soup.list.est)){
  rhos[[i]] <- mean(soup.list.est[[i]]$metaData$rho)
}
rhos <- do.call(rbind,rhos)
meta$soupx.rho <- rhos

# Save metadata file with the Rho values
write.table(x = meta, 
            file = "/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/SC_diestrus_metadirectory_final.csv",
            sep = ",",
            row.names = FALSE)

#Read in SoupX ####
seu.list <- list()
for(i in 1:length(meta$data.dir)){
  if(file.exists(paste0(meta$data.dir[i], '/outs/soupx'))){
    cat("Reading #",i, ": ", meta$data.dir[i], ' \n')
    seu.list[[i]] <- Seurat::Read10X(
      data.dir = paste0(meta$data.dir[i], '/outs/soupx')
    )
  }else{
    cat("Data not found for # ", i, " (", meta$data.dir[i], ")", "\n")
  }
}

# Seurat ####
seur.list <- list()

# Initialize Seurat objects
seur.list <- lapply( 
  seu.list,
  FUN = function(mat){
    return(CreateSeuratObject(
      counts = mat, 
      project = 'all_diestrus_mice'
    ))
  }
)  

#Add Metadata
for(i in 1:length(seur.list)){
  cat(' #####################################\n',
      '### Processing dataset number ', i, '###\n',
      '#####################################\n')
  # Add meta data
  for(md in colnames(meta)){
    seur.list[[i]][[md]] <- meta[[md]][i]
  }
  # add %MT (percent of features coming from mitochondrial genes)
  seur.list[[i]][["percent.mt"]]  <- PercentageFeatureSet(seur.list[[i]], pattern = "mt-") 
  
}

cat((sum(unlist(lapply(adj.mat.list, ncol)))-sum(unlist(lapply(seur.list, ncol)))),"cells (total) removed...\n")


# Preprocess seurat objects
seuPreProcess <- function(seu, assay='RNA', n.pcs=50, res=0.7){
  # NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  pca.name = paste0('pca_', assay)
  pca.key = paste0(pca.name,'_')
  umap.name = paste0('umap_', assay)
  
  seu = NormalizeData(
    seu
  ) %>% FindVariableFeatures(
    assay = assay,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = F
  ) %>% ScaleData(
    assay = assay
  ) %>% RunPCA(
    assay = assay,
    reduction.name = pca.name,
    reduction.key = pca.key,
    verbose = F,
    npcs = n.pcs
  )
  
  #find pcs to use
  tmp.var <- (seu@reductions[[pca.name]]@stdev)^2
  var.cut <- 0.95*sum(tmp.var)
  j=0
  var.sum = 0
  while(var.sum < 0.95*var.cut){
    j = j + 1
    var.sum <- var.sum + tmp.var[j]
  }
  n.pcs.use = j
  
  # FindNeighbors %>% RunUMAP, FindClusters
  seu <- FindNeighbors(
    seu,
    reduction = pca.name,
    dims = 1:n.pcs.use,
    force.recalc = TRUE,
    verbose = FALSE
  ) %>% RunUMAP(
    reduction = pca.name,
    dims = 1:n.pcs.use,
    reduction.name=umap.name
  )
  
  seu@reductions[[umap.name]]@misc$n.pcs.used <- n.pcs.use
  
  seu <- FindClusters(object = seu,resolution = res)
  seu[[paste0('RNA_res.',res)]] <- as.numeric(seu@active.ident)
  
  return(seu)
}

seur.list <- lapply(seur.list, seuPreProcess)

#DoubletFinder ####
bcmvn <- list()
pK <- list()
homotypic.prop <- list()
nExp_poi <- list()
nExp_poi.adj <- list()

# define the expected number of doublets - pK value calculation
for(i in 1:length(seur.list)){
  cat(' --------------------------------------------\n',
      '--- DoubletFinder for dataset number ', i, '---\n',
      '--------------------------------------------\n')
  
  ## pK Identification (no ground-truth)
  bcmvn[[i]]<- paramSweep_v3(             
    # pN-pK parameter sweeps on a 10,000-cell subset
    seu=seur.list[[i]], 
    PCs = 1:seur.list[[i]]@reductions$umap_RNA@misc$n.pcs.used, 
  ) %>% summarizeSweep(
    # computing the bimodality coefficient across pN and pK parameter space
    GT = FALSE
  ) %>% find.pK() 
  # Computes and visualizes the mean-variance normalized bimodality coefficient (BCmvn) score for each pK value tested during doubletFinder_ParamSweep
  
  ## Pull out max of bcmvn
  pK[[i]] <- as.numeric(as.character(bcmvn[[i]]$pK[bcmvn[[i]]$BCmetric==max(bcmvn[[i]]$BCmetric)])) # ugly, but functional...
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop[[i]] <- modelHomotypic(seur.list[[i]]$seurat_clusters) 
  
  nExp_poi[[i]] <- round(0.0242*length(colnames(seur.list[[i]])))  # Rate is 0.0242 for 5k "Target for Cell Recovery"
  nExp_poi.adj[[i]] <- round(nExp_poi[[i]]*(1-homotypic.prop[[i]]))
}

#Run Doublet finder
doublet <- function(seu) {
   seu = doubletFinder_v3(seu, 
                         pN = 0.25, #default value
                         pK = pK[[i]], 
                         nExp = nExp_poi.adj[[i]], 
                         PCs = 1:seur.list[[i]]@reductions$umap_RNA@misc$n.pcs.used)
  return(seu)
}

seur.list <- lapply(seur.list, FUN=doublet)

#Rename DF assigned columns for better merging
Col_rename <- function(seu) {
  colnames(seu@meta.data)[18] <- "DF.classifications"
  colnames(seu@meta.data)[17] <- "pANN"
  
  return(seu)
}

seur.list <- lapply(seur.list, FUN=Col_rename)

# Merging and Integration ####
# Merge datasets ####
# Add a prefix to the cell IDs so that cell IDs are NOT merged when the datasets are merged
Sample.ID <- c("mU07","mU13","mU14","mU15","mU16","mU17",
               'mU18','mU19','mU20','mU21','mU22','mU23','mU24',
               'mU26','mU27','mU28','mU29','mU30')

# Merge the datasets together
All.merged <- merge(seur.list[[1]], y = seur.list[2:length(Sample.ID)], add.cell.ids = Sample.ID[1:length(Sample.ID)])

# Save as RDS file
saveRDS(All.merged, file = "./data/DiestrusMice_mU7_mU30_Final_01252024.rds")

# Load RDS
All.merged <- readRDS(file = "./data/DiestrusMice_mU7_mU30_Final_01252024.rds", refhook = NULL)
ncol(All.merged) #42,694 cells

# Seurat workflow + Harmony integration ####
# Need to normalize and run PCA before running Harmony
All.merged <- NormalizeData(object = All.merged, assay = 'RNA')
All.merged <- FindVariableFeatures(object = All.merged, assay = 'RNA', 
                                   selection.method = 'vst', nfeatures = 2000)
All.merged <- ScaleData(object = All.merged, assay = 'RNA')
All.merged <- RunPCA(object = All.merged, assay = 'RNA', 
                     features = VariableFeatures(object = All.merged),
                     reduction.name = 'pca_RNA', reduction.key = 'pca_RNA_')

# Determine the 'dimensionality' of the dataset
ElbowPlot(All.merged,
          reduction = 'pca_RNA',
          ndims = 50) 

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

n.pcs=npcs(seu = All.merged, var.toal = 0.95, reduction = 'pca_RNA')
n.pcs
#39 Principle components

# Run Harmony
gc()
All.merged <- RunHarmony(All.merged, 
                         group.by.vars = "Sample.ID", 
                         reduction = "pca_RNA", 
                         assay = "RNA", 
                         plot_convergence = TRUE
                         )

# Normal Seurat clustering
# Set Python UMAP via reticulate
umap.method = 'umap-learn'
metric = 'correlation'

All.merged <- FindNeighbors(object = All.merged, 
                            reduction = "harmony",
                            k.param = 70, 
                            dims = 1:n.pcs, 
                            graph.name = 'harmony_snn')

All.merged <- FindClusters(object = All.merged, resolution = 0.7, graph.name = 'harmony_snn')
All.merged <- RunUMAP(object = All.merged, reduction = "harmony", dims = 1:n.pcs)

# UMAP 
DimPlot(object = All.merged, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5) + 
  labs(title = 'Colored by Cluster')

# Remove low_quality cells ####
Pre <- VlnPlot(All.merged, features=c('percent.mt',"nCount_RNA","nFeature_RNA"), group.by='Sample.ID')
Pre

HQ_Cells <- subset(All.merged, subset = percent.mt < 25 & nFeature_RNA > 200 & nCount_RNA > 750)
#38,298 cells remaining

Post <- VlnPlot(HQ_Cells, features=c('percent.mt', "nCount_RNA", "nFeature_RNA"), group.by='Sample.ID')
Post

saveRDS(All.merged,'./data/FinalSamples_preQC_02232024.rds') #save dataset before low-quality cells are removed (for later figures)
saveRDS(HQ_Cells,'./data/FinalSamples_preDF_02232024.rds') #save dataset with high-quality cells, but doublets still intact (for later figures)

#Check Doublets
#Violin Plot
VlnPlot(HQ_Cells, features = "nFeature_RNA", group.by = 'DF.classifications', pt.size = 0.1)

#UMAP
DFUMAP <- DimPlot(object = HQ_Cells,                 # Seurat object 
                  reduction = 'umap',
                  group.by = "seurat_clusters",     # Labels to color the cells by ("seurat_clusters", "Age", "Time.Point)  
                  repel = TRUE,                       # Whether to repel the cluster labels
                  label = TRUE,                       # Whether to have cluster labels 
                  pt.size = 0.5,                      # Size of each dot is (0.1 is the smallest)
                  label.size = 5) 
DFUMAP

DFUMAPFinder <- DimPlot(object = HQ_Cells,                 # Seurat object 
                        reduction = 'umap',
                        group.by = "DF.classifications",     # Labels to color the cells by ("seurat_clusters", "Age", "Time.Point)  
                        repel = TRUE,                       # Whether to repel the cluster labels
                        label = FALSE,                        # Whether to have cluster labels 
                        pt.size = 0.25,                      # Size of each dot is (0.1 is the smallest)
                        label.size = 5) 
DFUMAPFinder

#Remove Doublets ####
#Create Dataframe
Doublets = colnames(HQ_Cells@meta.data)[grepl("DF.classification", colnames(HQ_Cells@meta.data))]
HQ_Cells_DF = HQ_Cells[, HQ_Cells@meta.data[, Doublets] == "Singlet"]
ncol(HQ_Cells_DF) #37,543 cells remaining 

#Save RDS Objects with and without Doublets
saveRDS(HQ_Cells_DF, file = "./data/DiestrusMice_mU7_mU30_Final_01252024.rds")

#Run scaling etc. ####
HQ_Cells_DF = readRDS(file ="./data/DiestrusMice_mU7_mU30_Final_01252024.rds",  # Filename
                      refhook = NULL)
HQ_Cells_DF = ScaleData(HQ_Cells_DF, vars.to.regress = c("nFeature_RNA", "percent.mt"),   #RUN THESE AFTER FINDING AND REMOVING DOUBLETS
                        verbose = F)
HQ_Cells_DF = RunPCA(HQ_Cells_DF, verbose = F, npcs = n.pcs)
HQ_Cells_DF = RunUMAP(HQ_Cells_DF, dims = 1:10, verbose = F)

# Run Harmony
gc()
All.merged <- RunHarmony(All.merged, 
                         group.by.vars = "Sample.ID", 
                         reduction = "pca_RNA", 
                         assay = "RNA", 
                         plot_convergence = TRUE
)

#Clustering
HQ_Cells_DF <- FindNeighbors(object = HQ_Cells_DF,           # Seurat object
                             reduction = "harmony",         # Name of reduction to use
                             k.param = 90,                  # Defines k for the k-nearest neighbor algorithm (smaller the number the more clusters)
                             dims = 1:n.pcs,                   # Dimensions of reduction to use (use the result from npcs())
                             graph.name = 'AllMice.DF.Clusters')    # Name for stored SNN graph

HQ_Cells_DF <- FindClusters(object = HQ_Cells_DF,            # Seurat object 
                            resolution = 0.9,
                            graph.name = 'AllMice.DF.Clusters') # Resolution parameter (smaller the number the fewer clusters) 

# Set Python UMAP via reticulate
umap.method = 'umap-learn'
metric = 'correlation'

HQ_Cells_DF <- RunUMAP(object = HQ_Cells_DF,      # Seurat object 
                       reduction = 'harmony',    # Reduction to use
                       dims = 1:n.pcs)              # Dimensions of recution to use (use the result from npcs())

HQ_Cells_DFUMAP <- DimPlot(object = HQ_Cells_DF,                 # Seurat object 
                           reduction = 'umap',
                           group.by = "seurat_clusters",     # Labels to color the cells by ("seurat_clusters", "Age", "Time.Point)  
                           repel = TRUE,                       # Whether to repel the cluster labels
                           label = TRUE,                       # Whether to have cluster labels 
                           pt.size = 1,                      # Size of each dot is (0.1 is the smallest)
                           label.size = 5)                    # Font size for label
HQ_Cells_DFUMAP

#Overwrite for final filtered cells and clustering
saveRDS(HQ_Cells_DF, file = "./data/DiestrusMice_mU7_mU30_Final_01252024.rds")

# Add meta data - tdt postive/negative
#Find barcodes with Tdt+
TdT_expression = GetAssayData(object = HQ_Cells_DF, 
                              assay = "RNA", slot = "data")["TdTomato-UTR",]
pos_ids = names(which(TdT_expression>0))

## Get cell names
cellNames <- rownames(HQ_Cells_DF@meta.data)

## Mutate a column in original metadata
HQ_Cells_DF$barcode <- rownames(HQ_Cells_DF@meta.data)
HQ_Cells_DF@meta.data <- HQ_Cells_DF@meta.data %>% mutate(TdTomato = ifelse((HQ_Cells_DF$barcode %in% pos_ids), "Pos",  "Neg"))

#Check expression
VlnPlot(HQ_Cells_DF, features='TdTomato-UTR', group.by = 'TdTomato')

#Overwrite for final filtered cells and clustering
saveRDS(HQ_Cells_DF, file = "./data/DiestrusMice_mU7_mU30_Final_01252024.rds")

#Add 2 way coexpression data - Trp53 vs. Rb1 expression
# Add column of double neg/pos, etc. to object ####

# Dataframe of select genes
# Change the Seurat object and the gene IDs here
df <- data.frame("Trp53" = GetAssayData(object = HQ_Cells_DF, slot = "data")["Trp53",],
                      "Rb1" = GetAssayData(object = HQ_Cells_DF, slot = "data")["Rb1",])

#Add column with classification
for (i in 1:nrow(df)){
  if (df[i,1] > 0 & df[i,2] > 0){
    df[i,3] <- 'Double_Positive'
  }
  if (df[i,1] > 0 & df[i,2] == 0){
    df[i,3] <- 'Rb1_Negative'
  }
  if (df[i,1] == 0 & df[i,2] > 0){
    df[i,3] <- 'Trp53_Negative'
  }
  if (df[i,1] == 0 & df[i,2] == 0){
    df[i,3] <- 'Double_Negative'
  }
}

colnames(df)[3] <- 'Coexpression'

#Add new columns to the Seurat object
df_hash <- hash(keys=rownames(df), values=df$Coexpression)

HQ_Cells_DF$Coexpression <- NA

for(i in 1:nrow(HQ_Cells_DF@meta.data)){
  if((length(df_hash[[colnames(HQ_Cells_DF)[i]]]) == 1)==TRUE){
    HQ_Cells_DF@meta.data[i,23] <- df_hash[[colnames(HQ_Cells_DF)[i]]]
  }
}
View(HQ_Cells_DF$Coexpression)

saveRDS(HQ_Cells_DF, file = "./data/DiestrusMice_mU7_mU30_Final_01252024.rds")


#Add 3 way coexpression data - Tdt, Trp53, Rb1
# Add column of double neg/pos, etc. to object ####

# Dataframe of select genes
# Change the Seurat object and the gene IDs here
df <- data.frame("Trp53" = GetAssayData(object = HQ_Cells_DF, slot = "data")["Trp53",],
                 "Rb1" = GetAssayData(object = HQ_Cells_DF, slot = "data")["Rb1",],
                 "TdTomato" = GetAssayData(object = HQ_Cells_DF, slot = "data")["TdTomato-UTR",])

#Add column with classification
for (i in 1:nrow(df)){
  if (df[i,1] > 0 & df[i,2] > 0 & df[i,3] > 0){
    df[i,4] <- 'TdT+Trp53+Rb1+'
  }
  if (df[i,1] > 0 & df[i,2] > 0 & df[i,3] == 0){
    df[i,4] <- 'TdT-Trp53+Rb1+'
  }
  if (df[i,1] > 0 & df[i,2] == 0 & df[i,3] > 0){
    df[i,4] <- 'TdT+Trp53+Rb1-'
  }
  if (df[i,1] > 0 & df[i,2] == 0 & df[i,3] == 0){
    df[i,4] <- 'TdT-Trp53+Rb1-'
  }
  if (df[i,1] == 0 & df[i,2] > 0 & df[i,3] > 0){
    df[i,4] <- 'TdT+Trp53-Rb1+'
  }
  if (df[i,1] == 0 & df[i,2] > 0 & df[i,3] == 0){
    df[i,4] <- 'TdT-Trp53-Rb1+'
  }
  if (df[i,1] == 0 & df[i,2] == 0 & df[i,3] == 0){
    df[i,4] <- 'TdT-Trp53-Rb1-'
  }
  if (df[i,1] == 0 & df[i,2] == 0 & df[i,3] > 0){
    df[i,4] <- 'TdT+Trp53-Rb1-'
  }
}

colnames(df)[4] <- 'TdtExp_profile'

#Add new columns to the Seurat object
df_hash <- hash(keys=rownames(df), values=df$TdtExp_profile)

HQ_Cells_DF$TdtExp_profile <- NA

for(i in 1:nrow(HQ_Cells_DF@meta.data)){
  if((length(df_hash[[colnames(HQ_Cells_DF)[i]]]) == 1)==TRUE){
    HQ_Cells_DF@meta.data[i,24] <- df_hash[[colnames(HQ_Cells_DF)[i]]]
  }
}
View(HQ_Cells_DF$TdtExp_profile)

VlnPlot(HQ_Cells_DF, features='Rb1', group.by='Coexpression')

saveRDS(HQ_Cells_DF, file = "./data/DiestrusMice_mU7_mU30_Final_mergedCellID_02122024.rds")

#Sample Analysis####
#Load RDS File
IntData = readRDS(file ="./data/DiestrusMice_mU7_mU30_Final_mergedCellID_02122024.rds",  # Filename
                  refhook = NULL)
ncol(IntData) #37,543 Cells

#Summary of Markers ####
Sumfeatures = c('Epcam','Pax8','TdTomato-UTR', #General Epithelium
                'Tacstd2','Prap1','Krt13', 'Klf5','Lpar3','Ly6a','Met','Sprr2f','Slc1a3','Wnt7a','Wnt7b','Ptgs1', #Luminal Epithelium
                'Foxa2', 'Aldh1a1','Axin2','Lgr5','Ltf','Prss29','Sox17','Lgr4','Muc1','Ovgp1','Wif1','Slc18a2', #Glandular Epithelium
                'Cdkn2a','Ccne1','Cdk1','Cdk4','Cdkn1c','Pcna','Mki67', #Cycling genes
                "Adamts1",'Col1a1','Col1a2','Col3a1','Creb5','Dcn','Dio2','Gfpt2','Igf1','Lrrn4','Medag','Mgp','Msln','Pdgfra','Rxfp1','Vcan','Vim', #Fibroblasts
                'Acta1','Des',"Myh11",'Acta2','Palld', # Muscle
                'Lyve1','Mmrn1','Flt4', #lymphatic endothelium
                'Cavin2','Cldn5','Cdh5','Icam1','Pecam1','Tek','Vwf' #Vascular Endothelium
                
)
SumDot = DotPlot(IntData, features=Sumfeatures, col.min=0, col.max=2.5) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic'))
SumDot

#Stem Cell Markers ####
SCfeatures = c('Aldh1a1','Ly6a','Cd44','Ezh2','Fut4',
               'Itga6','Klf4','Klf5','Klf6','Lgr5','Lrig1','Neat1','Pax8','Foxa2','Wif1','Prom1', 'Pou5f1',
               'Slc1a3','Slc38a2','Sox2','Sox9','Sox17','Tacstd2','Zfp36')
DotSC = DotPlot(IntData, features=SCfeatures, col.min=0, col.max=2.5) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic'))
DotSC

#Diff Gene Expression - Progenitor Populations
DotPlot(IntData, features=c('Pax2','Hoxb5','Tert','Sox17','Klf5', #nucleus
                            'Nt5e','Lpar3','Lgr4', #membrane
                            'Wfdc2' #Secreted
), col.min=0, col.max=2.5) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic'))

#Immune Summaries ####
Imfeatures = c( 'Cd14','Csf3r','Itgal','Itgam','Lrp1','Ccr2', 'Fcgr3','Sirpa','Cd33', #monocytes
                'Adgre1','Aif1','C1qa','Cd74','Csf1r','Fcgr1','Mrc1', #macrophage
                'Ptprc', #Leukocytes
                'Cd4','Cd3e', #helper T cells
                'Il2ra','Foxp3', #regulatory T cells
                'Cd8a','Pdcd1',  #cytotoxic T Cells
                'Ncam1','Pecam1','Tnfrsf8','Cd38','Klrb1',
                'Ccl5', 'Klrd1','Txk', 'Xcl1', #NK Cells
                'Cd24a','Cd22','Cd19','Ms4a1' #B Cells
)
ImDot = DotPlot(IntData, features=Imfeatures, col.min=0, col.max=2.5) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face='italic')) #angle x axis labels
ImDot

#Cluster Renaming ####
#Fine Grain
Idents(IntData) <- IntData$seurat_clusters

CellTypeNames <- RenameIdents(IntData, 
                              '0'= 'Epi Progenitor',
                              '1'= 'Lymph',
                              '2'= 'Fib 1',
                              '3'= 'Cycling',
                              '4'= 'LE 1',
                              '5'= 'Fib 3',
                              '6'= 'Fib 2',
                              '7'= 'Macrophage',
                              '8'= 'Vascular',
                              '9'= 'GE',
                              '10'= 'Lymphocyte',
                              '11'= 'Mesothelium',
                              '12'= 'Smooth Muscle',
                              '13'= 'LE 2'
)

IntData$seurat_clusters2 <- Idents(CellTypeNames)

Idents(IntData) <- IntData$seurat_clusters2

#Course Grain
Idents(IntData) <- IntData$seurat_clusters2

CellTypeNames <- RenameIdents(IntData, 
                              'Epi Progenitor'='EP',
                              'LE 1'='LE',
                              'LE 2'='LE',
                              'Cycling'='CE',
                              'GE'='GE',
                              'Fib 1'='Fib',
                              'Fib 2'='Fib',
                              'Fib 3'='Fib',
                              'Mesothelium'='Meso',
                              'Smooth Muscle'='SM',
                              'Lymph'='LV',
                              'Vascular'='BV',
                              'Macrophage'='Im',
                              'Lymphocyte'='Im'
                              
)

IntData$seurat_clusters_CC <- Idents(CellTypeNames)

#Name by broad cell type - for later CellChat Analysis
Idents(IntData) <- IntData$seurat_clusters2

CellTypeNames <- RenameIdents(IntData, 
                              'Epi Progenitor'='Epi',
                              'LE 1'='Epi',
                              'LE 2'='Epi',
                              'Cycling'='Epi',
                              'GE'='Epi',
                              'Fib 1'='Fib',
                              'Fib 2'='Fib',
                              'Fib 3'='Fib',
                              'Mesothelium'='Meso',
                              'Smooth Muscle'='SM',
                              'Lymph'='LV',
                              'Vascular'='BV',
                              'Macrophage'='Im',
                              'Lymphocyte'='Im'
                              
)

IntData$seurat_clusters_cellType <- Idents(CellTypeNames)
saveRDS(IntData, file ="./data/DiestrusMice_mU7_mU30_Final_mergedCellID_02122024.rds")

#UMAP ####
UMAP = DimPlot(object=IntData, reduction="umap", group.by = "seurat_clusters2", 
               repel = TRUE,                     
               label = FALSE,
               pt.size = 0.05,    
               label.size = 4,
               shuffle=TRUE
)
UMAP
