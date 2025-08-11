################################################################################
##############################Supp Figure 6#####################################
#Changes in epithelial gene expression across SEC stages with spatial location##
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

#Load All Sample epithelial subset ####
EpiOnly = readRDS(file ='./data/allDiestrus_Epi_mU7_mU30_recluster_final_01182025.rds',
                  refhook = NULL)
ncol(EpiOnly) #19,449 cells

#Set order of SEC stages
SECord= c('control', 'earlySEC','lateSEC') #(Normal, Pre-dysplastic, Dysplastic)

#Supp. Fig 5 - Epithelial gene expression compared across SEC stages####
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


DotPlot(EpiOnly, features=rev(AllFeatures_space), col.min=0, col.max=2.5, 
        group.by='RedSEC_stage'
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=SECord)+
  coord_flip()
ggsave('./plots/DotPlot_ALLEpi_IntDataEPI_SECstage.pdf', last_plot(), device='pdf',width=6, height=12,units='in', dpi=300)

