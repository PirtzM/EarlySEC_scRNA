################################################################################
#################################Figure 2#######################################
#######Census of Epithelial Cell types in the mouse uterus-stages of SEC########
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

#Load each SEC stage Epithelial Subsets ####
#Integrated samples, n=18
EpiOnly = readRDS(file ="./data/allDiestrus_Epi_mU7_mU30_recluster_final_01182025.rds",
                  refhook = NULL)
ncol(EpiOnly) #19,449 cells
Idents(EpiOnly) <- EpiOnly$seurat_clusters_CC #Set Active Identity

#Normal samples, n=5
ContEpi = readRDS(file ="./data/FinalSamples/cont_Epi_mU7_mU30_Final_01182025.rds",  
                   refhook = NULL)
ncol(ContEpi) #4,689 cells
Idents(ContEpi) <- ContEpi$seurat_clusters_CC #Set Active Identity

#Pre-dysplastic samples, n=7
EarEpi = readRDS(file ="./data/Ear_Epi_mU7_mU30_recluster_01182025.rds",
                  refhook = NULL)
ncol(EarEpi) #7,713 cells
Idents(EarEpi) <- EarEpi$seurat_clusters_CC #Set Active Identity

#Dysplastic samples, n=6
LEpi = readRDS(file ="./data/Late_Epi_mU7_mU30_recluster_01182025.rds",  # Filename
                refhook = NULL)
ncol(LEpi) #7,063 cells
Idents(LEpi) <- LEpi$seurat_clusters_CC #Set Active Identity

#Plotting variable order
clusterord_CC= c('LE','GE',
                 'DDP',
                 'CE','CM')

#Supp. Fig 6 - comparison of the broad expression profiles of epithelial clusters between the integrated and SEC-specific dataset DotPlots (Final merging of DotPlots was done in Adobe Illustrator)
#Gene List
MidiFeatures_org = c('Ly6a','Dusp1','Wfdc2', 'Klf4', #putative progenitor
                     'Tacstd2','Neat1','Met','Prap1', #LE
                     'Foxa2','Msx1','Prom1','Lgr5', #GE
                     'Pde4c','Dnah12','Dand5', #DDP
                     'Mki67','Top2a','Cdkn2a', #Cycling
                     'COX1','COX2','Malat1' #CM
)
#Integrated epithelial Dataset DotPlot
midi_EPI <- DotPlot(EpiOnly, features=MidiFeatures_org, col.min=0, col.max=2.5, scale.min=0
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0, limits=c(0,2))+ #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_CC)+
  coord_flip()
midi_EPI

ggsave('./plots/FinalSamples/DotPlot_EpiOnly_int_midi-list_comparison_06262025.pdf', last_plot(), device='pdf',width=6, height=7,units='in', dpi=300)

#Normal epithelial dataset DotPlot
midi_EPI <- DotPlot(ContEpi, features=MidiFeatures_org, col.min=0, col.max=2.5, scale.min=0
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0, limits=c(0,2))+ #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_cont_CC)+
  coord_flip()
midi_EPI

ggsave('./plots/FinalSamples/DotPlot_EpiOnly_Control_midi-list_comparison_06262025.pdf', last_plot(), device='pdf',width=6, height=7,units='in', dpi=300)

#Pre-dysplastic epithelial dataset DotPlot
midi_EPI <- DotPlot(EarEpi, features=MidiFeatures_org, col.min=0, col.max=2.5, scale.min=0
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0, limits=c(0,2))+ #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_CC)+
  coord_flip()
midi_EPI

ggsave('./plots/FinalSamples/DotPlot_EpiOnly_Early_midi-list_comparison_06262025.pdf', last_plot(), device='pdf',width=6, height=7,units='in', dpi=300)

#Dysplastic epithelial dataset DotPlot
midi_EPI <- DotPlot(LEpi, features=MidiFeatures_org, col.min=0, col.max=2.5, scale.min=0
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0, limits=c(0,2))+ #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_CC)+
  coord_flip()
midi_EPI

ggsave('./plots/FinalSamples/DotPlot_EpiOnly_late_midi-list_comparison_06262025.pdf', last_plot(), device='pdf',width=6, height=7,units='in', dpi=300)
