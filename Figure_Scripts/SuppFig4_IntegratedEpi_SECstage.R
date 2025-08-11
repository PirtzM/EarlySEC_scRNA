################################################################################
##############################Supp Figure 4#####################################
#########Changes in Integrated epithelial cell types across SEC stages##########
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
EpiOnly = readRDS(file ="./data/allDiestrus_Epi_mU7_mU30_recluster_final_01182025.rds",
                  refhook = NULL)
ncol(EpiOnly) #19,449 cells
Idents(EpiOnly) <- EpiOnly$seurat_clusters_SE #Set Active Identity

# Color Palettes ####
#Integrated epithelial dataset cell types
mycols_SE=c('LE'='#066DD4',
            'GE'='#8F03C8',
            'Cycling'='#E86DE3',
            'COX/MAL'='#5BB40B',
            'DDP'='#9f2b68'
)

#SEC stage colors
SECcols=c('control'='#031273',
          'earlySEC'='#EDCD44',
          'lateSEC'='#DC3E26')

#Plotting variable order
clusterord_SE= c('LE','GE','DDP',
                 'Cycling','COX/MAL')


#Supp. Fig 4A - UMAP of the integrated epithelial dataset, coded by cell type####
UMAP = DimPlot(object=EpiOnly, reduction="umap", group.by = "seurat_clusters_SE",              
               repel = FALSE,                     
               label =FALSE,
               pt.size = 0.75,                
               label.size = 4,
               shuffle=TRUE,
               cols=mycols_SE
  ) + NoLegend()
UMAP 

ggsave('./plots/FinalSamples/IntEpiUMAP_clusters_06262025.pdf', UMAP, device='pdf', width=6, height=4, units='in', dpi=300)

#Supp. Fig 4B - UMAP of the integrated epithelial dataset, highlighting normal cells####
UMAP = DimPlot(object=EpiOnly, reduction="umap", group.by = "RedSEC_stage",              
               repel = FALSE,                     
               label =FALSE,
               pt.size = 0.75,                
               label.size = 4,
               shuffle=TRUE,
               cols=c('Normal'='#031273',
                      'Pre-dysplastic'='#bbbbbb'
                      'Dysplastic'='#bbbbbb')
  ) + NoLegend()
UMAP 

ggsave('./plots/FinalSamples/IntEpiUMAP_normal_06262025.pdf', last_plot(), device='pdf', width=6, height=4, units='in', dpi=300)

#Supp. Fig 4C - UMAP of the integrated epithelial dataset, highlighting pre-dysplastic cells####
UMAP = DimPlot(object=EpiOnly, reduction="umap", group.by = "RedSEC_stage",              
               repel = FALSE,                     
               label =FALSE,
               pt.size = 0.75,                
               label.size = 4,
               shuffle=TRUE,
               cols=c('Normal'='#bbbbbb',
                      'Pre-dysplastic'='#EDCD44'
                      'Dysplastic'='#bbbbbb')
  ) + NoLegend()
UMAP 

ggsave('./plots/FinalSamples/IntEpiUMAP_PD_06262025.pdf', last_plot(), device='pdf', width=6, height=4, units='in', dpi=300)

#Supp. Fig 4D - UMAP of the integrated epithelial dataset, highlighting dysplastic cells####
UMAP = DimPlot(object=EpiOnly, reduction="umap", group.by = "RedSEC_stage",              
               repel = FALSE,                     
               label =FALSE,
               pt.size = 0.75,                
               label.size = 4,
               shuffle=TRUE,
               cols=c('Normal'='#bbbbbb',
                      'Pre-dysplastic'='#bbbbbb'
                      'Dysplastic'='#DC3E26')
  ) + NoLegend()
UMAP 

ggsave('./plots/FinalSamples/IntEpiUMAP_D_06262025.pdf', last_plot(), device='pdf', width=6, height=4, units='in', dpi=300)

#Supp. Fig 4E - DotPlots of epithelial cell expression profiles, split by SEC stage (final plot merging done in Adobe Illustrator####
#Change identity for identity splitting
Idents(EpiOnly) <- EpiOnly$RedSEC_stage

#Gene list
MiniFeatures_org = c('TdTomato-UTR', #General Epithelium
                     'Tacstd2','Met',
                     'Foxa2','Msx1','Prom1',
                     'Dnah12',
                     'Mki67',
                     'COX1','Malat1'
)

#Normal DotPlot
mini_EPI <- DotPlot(EpiOnly, features=MiniFeatures_org, col.min=0, col.max=2.5,
                    group.by='seurat_clusters_SE', 
                    idents='Normal'
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0, limits=c(0,2))+# + #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_SE)+
  coord_flip()
mini_EPI

ggsave('./plots/FinalSamples/DotPlot_IntData_EpiOnly_Norm_clusters_MiniEpi_06262025.pdf', last_plot(), device='pdf',width=6, height=6,units='in', dpi=300)

#Pre-dysplastic DotPlot
mini_EPI <- DotPlot(EpiOnly, features=MiniFeatures_org, col.min=0, col.max=2.5,
                    group.by='seurat_clusters_SE', 
                    idents='Pre-dysplastic'
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0, limits=c(0,2))+# + #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_SE)+
  coord_flip()
mini_EPI

ggsave('./plots/FinalSamples/DotPlot_IntData_EpiOnly_PD_clusters_MiniEpi_06262025.pdf', last_plot(), device='pdf',width=6, height=6,units='in', dpi=300)

#Dysplastic DotPlot
mini_EPI <- DotPlot(EpiOnly, features=MiniFeatures_org, col.min=0, col.max=2.5,
                    group.by='seurat_clusters_SE', 
                    idents='Dysplastic'
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0, limits=c(0,2))+# + #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_SE)+
  coord_flip()
mini_EPI

ggsave('./plots/FinalSamples/DotPlot_IntData_EpiOnly_D_clusters_MiniEpi_06262025.pdf', last_plot(), device='pdf',width=6, height=6,units='in', dpi=300)
