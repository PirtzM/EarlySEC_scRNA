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
EpiOnly = readRDS(file ='./data/allDiestrus_Epi_mU7_mU30_recluster_final_01302024_noLE3.rds',
                  refhook = NULL)
ncol(EpiOnly) #19,449 cells

#Set order of SEC stages
SECord= c('control', 'earlySEC','lateSEC') #(Normal, Pre-dysplastic, Dysplastic)

####Supp. Fig 6 - Epithelial gene expression compared across SEC stages with their spatial locations####
AllFeatures_space = c('Tacstd2','Prap1','Brca2','Cfb','Ckmt1','Clca1','Crabp2','Cxcl17',
                      'Irf7','Lcn2','Mogat1','Morrbid','Muc4','Oas2','Pla2g2e','Prim1', #LE
                      'Btg2','Dock11','Fat3','Met','Oasl2','Rb1','Trp53','Muc1','Padi1', #Both GE & LE
                      'Foxa2','Foxj1','Axin2','Lgr5','Prom1','Aoc1',#Glandular Epithelium
                      'Cdkn2a','Mki67',#'Ccne1','Cdk1','Cdk4','Pcna', #Cycling genes
                      'Pgr','Esr1' #hormonally regulated genes
)


DotPlot(EpiOnly, features=rev(AllFeatures_space), col.min=0, col.max=2.5, 
        group.by='RedSEC_stage'
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=SECord)+
  coord_flip()
ggsave('./plots/DotPlot_ALLEpi_Spatial_IntDataEPI_SECstage.pdf', last_plot(), device='pdf',width=6, height=12,units='in', dpi=300)

