################################################################################
#################################Supp Fig. 7####################################
###########Identification of epithelial clusters between SEC stages#############
################################################################################

#R v4.1.1
#Load Libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)

#Set working directory
setwd("/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/scripts")

#Load datasets ####
#Normal samples, n=5
ContEpi = readRDS(file ="./data/FinalSamples/cont_Epi_mU7_mU30_Final_01182025.rds",  
                   refhook = NULL)
ncol(ContEpi) #4,689 cells
Idents(ContEpi) <- ContEpi$seurat_clusters2 #Set Active Identity

#Pre-dysplastic samples, n=7
EarEpi = readRDS(file ="./data/Ear_Epi_mU7_mU30_recluster_01182025.rds",
                  refhook = NULL)
ncol(EarEpi) #7,713 cells
Idents(EarEpi) <- EarEpi$seurat_clusters2 #Set Active Identity

#Dysplastic samples, n=6
LEpi = readRDS(file ="./data/Late_Epi_mU7_mU30_recluster_01182025.rds",  # Filename
                refhook = NULL)
ncol(LEpi) #7,063 cells
Idents(LEpi) <- LEpi$seurat_clusters2 #Set Active Identity

#Plotting variable order ####
clusterord_cont= c('N - LE 1','N - LE 2',
                   'N - LE 3','GE',
                   'DDP','Cycling','COX/MAL') #Normal sample epithelial cell type order

clusterord_ear= c('PD - LE 1','PD - LE 2',
                  'PD - LE 3','PD - LE 4', 'GE', 
                  'DDP','Cycling','COX/MAL') #Pre-dysplastic sample epithelial cell type order

clusterord_L= c('D - LE 1','D - LE 2','D - LE 3',
                'D - LE 4','D - LE 5','GE',
                'Cycling', 'COX/MAL') #Dysplastic sample epithelial cell type order

#Supp. Fig 7 - Epithelial SEC stage-specific expression profile DotPlots (All final merging of plots were done in Adobe Illustrator####
#Gene List
AllFeatures_space = c('TdTomato-UTR','Epcam','Pax8', #General Epithelium
                      'Esr1','Pgr', #hormone receptors
                      'Ly6a','Klf4','Sox17','Tert','Nt5e','Hoxb5',
                      'Cd44','Itga6','Sox2','Lrig1','Klf6','Sox9',
                      'Tcf4','Dusp1','Btg2','Fos','Pax2','Wfdc2', #Putative Progenitor
                      'Tacstd2','Neat1','Met','Prap1','Krt13','Klf5','Sprr2f','Crabp2', #Luminal Epithelium
                      'Foxa2','Msx1','Gpx2','Wnt5a','Gstm7','Prom1','Aldh1a1','Axin2','Lgr5', 'Wif1', #Glandular Epithelium
                      'Foxj1','Pde4c','Dnah12','Dand5', #PD/DN
                      'Mki67','Cdkn2a','Top2a','Pcna', #Cycling genes
                      'COX1','COX2', 'Malat1' #CM                    
)

#Normal DotPlot
EPI <- DotPlot(ContEpi, features=AllFeatures_org, col.min=0, col.max=2.5
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0, limits=c(0,2.75))+ #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_cont)+
  coord_flip()
EPI

ggsave('./plots/FinalSamples/DotPlot_EpiOnly_Control_clusters_expEpi_06262025.pdf', last_plot(), device='pdf',width=6, height=12,units='in', dpi=300)

#Pre-dysplastic DotPlot
EPI <- DotPlot(EarEpi, features=AllFeatures_org, col.min=0, col.max=2.5
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0, limits=c(0,2.75))+ #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_ear)+
  coord_flip()
EPI

ggsave('./plots/FinalSamples/DotPlot_EpiOnly_Early_clusters_expEpi_06262025.pdf', last_plot(), device='pdf',width=6, height=12,units='in', dpi=300)

#Dysplastic DotPlot
EPI <- DotPlot(LEpi, features=AllFeatures_org, col.min=0, col.max=2.5
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0, limits=c(0,2.75))+ #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_L)+
  coord_flip()
EPI

ggsave('./plots/FinalSamples/DotPlot_EpiOnly_Late_clusters_expEpi_06262025.pdf', last_plot(), device='pdf',width=6, height=12,units='in', dpi=300)
