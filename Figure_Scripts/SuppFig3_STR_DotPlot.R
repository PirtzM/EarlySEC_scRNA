################################################################################
##############################Supp Figure 3#####################################
#################Dot Plots of cell type expression for all cells################
################################################################################


#R v4.1.1
#Load Libraries
library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)

#Load integrated RDS File ####
IntData = readRDS(file ="./data/DiestrusMice_mU7_mU30_Final_01182025_simple.rds",  # Filename
                  refhook = NULL)
ncol(IntData) #37,543 Cells
Idents(IntData) <- IntData$seurat_clusters_SE #Set Active Identity

#Plotting variable order ####
clusterord_int_SE= c('Epithelium',
                  'Mesothelium','Fib 1','Fib 2','Fib 3',
                  'Smooth Muscle','Lymph','Vascular', 'Macrophage','Lymphocyte')

#Gene list ####
AllFeatures_STR = c('TdTomato-UTR','Epcam','Pax8','Tacstd2','Foxa2',
                    'Msln','Lrrn4','Upk3b', #meso
                    'Pdgfra', #general Fib
                    'Adamts1', #fib1
                    'Pi16', #fib2
                    'Col6a4','Rxfp1',"Zyg11a", #Fibroblasts 
                    'Acta2', # Muscle
                    'Lyve1','Mmrn1','Flt4', #lymphatic endothelium
                    'Pecam1','Tek','Vwf', #Vascular Endothelium
                    'Adgre1','Fcgr1','Fcgr3','Cd3e','Cd4','Xcl1'
)

#DotPlot
STR <- DotPlot(IntData, features=c(AllFeatures_STR), col.min=0, col.max=2.5,
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=c(clusterord_int_SE
                            ))+
  coord_flip()
STR

ggsave('./plots/FinalSamples/DotPlot_ALL_IntData.pdf', last_plot(), device='pdf',width=9, height=15,units='in', dpi=300)
