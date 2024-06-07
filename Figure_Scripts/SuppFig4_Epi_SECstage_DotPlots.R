################################################################################
##############################Supp Figure 4#####################################
############Changes in epithelial gene expression across SEC stages#############
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

####Supp. Fig 4 - Epithelial gene expression compared across SEC stages####
#All Epi
AllFeatures_epi = c('TdTomato-UTR','Epcam','Pax8', #General Epithelium
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

DotPlot(EpiOnly, features=c(rev(AllFeatures_epi)), col.min=0, col.max=2.5, 
        group.by='RedSEC_stage') +
        scale_color_viridis(option = 'A', begin=0.9, end=0) + #change to magma color scale
        theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
        scale_y_discrete(limits=c(SECord))+
        coord_flip()
ggsave('./plots/DotPlot_EpiOnly_IntData_SECstage.pdf', last_plot(), device='pdf',width=9, height=15, units='in', dpi=300)
