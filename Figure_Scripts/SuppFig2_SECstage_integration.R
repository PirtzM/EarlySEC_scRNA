################################################################################
#############Supp Figure 2 - UMAP grouped by SEC stage (18 samples)#############
################################################################################

#R v4.1.1
#Load Libraries
library(Seurat)
library(ggplot2)
library(viridis)

#Set working directory
setwd("/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/scripts")

#Load color palette
SECcols=c('control'='#031273',
          'earlySEC'='#EDCD44',
          'lateSEC'='#DC3E26')

#Load All Sample Dataset
IntData = readRDS(file ="./data/DiestrusMice_mU7_mU30_Final_mergedCellID_02122024.rds",  # Filename
                  refhook = NULL)
ncol(IntData) #37,543 

#Create UMAP to visualize how well SEC stages integrated with each other
UMAP = DimPlot(object=IntData, reduction="umap", group.by = "RedSEC_stage",
               repel = TRUE,                   
               label = FALSE,
               pt.size = 0.05,                     
               label.size = 4,
               cols=SECcols,
               shuffle=TRUE
)
UMAP

ggsave('./plots/UMAP_SECintegration.pdf', UMAP, device='pdf', width=6, height=4, units='in',dpi=300)
