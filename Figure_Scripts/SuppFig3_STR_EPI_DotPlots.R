################################################################################
##############################Supp Figure 3#####################################
#######Dot Plots of cell type expression for stromal and epithelial cells#######
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


####Load each SEC stage - All Cells####
#Normal samples, n=5
ContOnly = readRDS(file ="./data/cont_mU7_mU30_Final_01262024.rds",  
                   refhook = NULL)
ncol(ContOnly) #7614 cells
Idents(ContOnly) <- ContOnly$seurat_clusters2 #Set Active Identity

#Pre-dysplastic samples, n=7
EarOnly = readRDS(file ="./data/EarCan_recluster_mU7_mU30_10282023.rds",
                  refhook = NULL)
ncol(EarOnly) #17,212 cells
Idents(EarOnly) <- EarOnly$seurat_clusters2 #Set Active Identity

#Dysplastic samples, n=6
LOnly = readRDS(file ="./data/LateCan_recluster_mU7_mU30_10292023.rds",  # Filename
                refhook = NULL)
ncol(LOnly) #12,717 cells

Idents(LOnly) <- LOnly$seurat_clusters2 #Set Active Identity

####Load each SEC stage - Epithelial Subsets####
#Normal samples, n=5
ContEpi = readRDS(file ="./data/FinalSamples/cont_Epi_mU7_mU30_Final_01262024_noMeso.rds",  
                  refhook = NULL)
ncol(ContEpi) #4,689 cells
Idents(ContEpi) <- ContEpi$seurat_clusters2 #Set Active Identity

#Pre-dysplastic samples, n=7
EarEpi = readRDS(file ="./data/Ear_Epi_mU7_mU30_recluster_11212023_withUECE.rds",
                 refhook = NULL)
ncol(EarEpi) #7,713 cells
Idents(EarEpi) <- EarEpi$seurat_clusters2 #Set Active Identity

#Dysplastic samples, n=6
LEpi = readRDS(file ="./data/Late_Epi_mU7_mU30_recluster_01282024_CE.rds",  # Filename
               refhook = NULL)
ncol(LEpi) #7,063 cells
Idents(LEpi) <- LEpi$seurat_clusters2 #Set Active Identity

#Set Cell Cluster Order####
clusterord_conSTR= c('Mesothelium','Fib 1','Fib 2','Smooth Muscle', 'Lymph','Vascular',
                     'Macrophage','T cell','NK cell') #Normal sample cell type order

clusterord_earSTR= c('Mesothelium','Fib 1','Fib 2','Fib 3','Smooth Muscle',
                     'Lymph','Vascular','Macrophage','Lymphocyte') #Pre-dysplastic sample cell type order

clusterord_LSTR= c('Mesothelium','Fib 1','Fib 2','Fib 3', 
                   'Lymph','Vascular','Macrophage','Lymphocyte') #Dysplastic sample cell type order

clusterord_contepi= c('Progenitor','LE 1','LE 2','GE','Epithelial - Foxj1+',
                      'Cycling','Unidentified') #Normal sample epithelial cell type order

clusterord_earepi= c('Progenitor 1','Progenitor 2','LE 1','LE 2','GE','Epithelial - Foxj1+',
                     'Cycling','Unidentified') #Pre-dysplastic sample epithelial cell type order

clusterord_Lepi= c('Progenitor 1','Progenitor 2','Progenitor 3','LE 1','LE 2','GE',
                   'Epithelial - Foxj1+','Cycling')

#List genes of interest ####
#stromal genes
STRFeatures = c('Msln','Dcn','Col3a1','Pdgfra',"Adamts1",'Gfpt2','Medag',#'Lrrn4','Mgp','Vim',
                'Vcan','Dio2','Rxfp1','Igf1','Zyg11a', #Fibroblasts
                'Palld','Des',"Myh11",'Acta2', # Muscle
                'Lyve1','Mmrn1','Flt4', #lymphatic endothelium
                'Cavin2','Cldn5','Cdh5','Icam1','Pecam1','Tek','Vwf', #Vascular Endothelium
                'Adgre1','Fcgr1','Fcgr3', #macrophage
                'Cd3e','Cd4','Xcl1' #lymphocytes
)

EPIFeatures = c('TdTomato-UTR','Epcam','Pax8', #General Epithelium
                'Ly6a','Klf4','Sox17','Tert','Nt5e','Hoxb5',
                'Cd44','Itga6','Sox2','Lrig1','Klf6','Sox9',
                'Tcf4','Fosb','Dusp1','Btg2','Nedd9','Fos','Jun','Pax2','Wfdc2', #progenitor-like
                'Tacstd2','Neat1','Met','Prap1','Krt13','Klf5','Lpar3','Sprr2f', #Luminal Epithelium
                'Wnt7a','Wnt7b','Oas2','Oasl2', #LE
                'Foxa2','Aldh1a1','Prom1','Foxj1','Axin2','Lgr5','Ltf',
                'Lgr4','Muc1','Ovgp1','Wif1', #Glandular Epithelium
                'Cdkn2a','Ccne1','Cdk1','Cdk4','Pcna','Mki67', #Cycling genes
                'Malat1','COX1','COX2','percent.mt'
)

####Supp. Fig 3A - stromal expression in each SEC stage ####
#Normal samples
DotPlot(ContOnly, features=c(rev(STRFeatures)), col.min=0, col.max=2.5, 
        idents=c('Mesothelium','Fib 1','Fib 2','Smooth Muscle', 'Lymph','Vascular',
                 'Macrophage','T cell','NK cell')) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.y=element_text(face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_conSTR)+
  coord_flip()
ggsave('./plots/DotPlot_STRIM_control.pdf', last_plot(), device='pdf',width=6, height=8, units='in', dpi=300)

#Pre-dysplastic samples
DotPlot(EarOnly, features=c(rev(STRFeatures)), col.min=0, col.max=2.5, 
        idents=c('Mesothelium','Fib 1','Fib 2','Fib 3','Smooth Muscle',
                 'Lymph','Vascular','Macrophage','Lymphocyte')) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.y=element_text(face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_earSTR)+
  coord_flip()
ggsave('./plots/DotPlot_STRIM_early.pdf', last_plot(), device='pdf',width=6, height=8, units='in', dpi=300)

#Dysplastic samples
DotPlot(LOnly, features=c(rev(STRFeatures)), col.min=0, col.max=2.5, 
        idents=c('Mesothelium','Fib 1','Fib 2','Fib 3', 
                 'Lymph','Vascular','Macrophage','Lymphocyte')) +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.y=element_text(face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_LSTR)+
  coord_flip()
ggsave('./plots/DotPlot_STRIM_late.pdf', last_plot(), device='pdf',width=6, height=8, units='in', dpi=300)

####Supp. Fig 3B - epithelial expression in each SEC stage ####
#Normal Samples
DotPlot(ContEpi, features=c(rev(EPIFeatures)), col.min=0, col.max=2.5, 
        group.by='seurat_clusters2') +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.y=element_text(face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=c(clusterord_contepi))+
  coord_flip()
ggsave('./plots/DotPlot_EpiOnly_control.pdf', last_plot(), device='pdf',width=6, height=12, units='in', dpi=300)

#Pre-dysplastic Samples
DotPlot(EarEpi, features=c(rev(EPIFeatures)), col.min=0, col.max=2.5, 
        group.by='seurat_clusters2') +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.y=element_text(face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=c(clusterord_earepi))+
  coord_flip()
ggsave('./plots/DotPlot_EpiOnly_early.pdf', last_plot(), device='pdf',width=6, height=12, units='in', dpi=300)

#Dysplastic Samples
DotPlot(EarEpi, features=c(rev(EPIFeatures)), col.min=0, col.max=2.5, 
        group.by='seurat_clusters2') +
  scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
  theme(axis.text.y=element_text(face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=c(clusterord_Lepi))+
  coord_flip()
ggsave('./plots/DotPlot_EpiOnly_late.pdf', last_plot(), device='pdf',width=6, height=12, units='in', dpi=300)
