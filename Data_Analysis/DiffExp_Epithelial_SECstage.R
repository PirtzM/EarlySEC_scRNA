################################################################################
###################Analysis of Final Samples (18 total)#########################
###################Epithelial Differential Expression###########################
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

#Load full object ####
EpiOnly = readRDS(file ='./data/allDiestrus_Epi_mU7_mU30_recluster_final_01302024_noLE3.rds',
                  refhook = NULL)
ncol(EpiOnly) #19,449 cells

Idents(EpiOnly) <- Epi_Only$RedSEC_stage #Set Active Identity

#Differential Gene Expression####
#Epithelial - early####
StatusDiff <- FindMarkers(EpiOnly, ident.1="earlySEC", ident.2="control")

EarEpi <- StatusDiff %>%
  arrange(desc(avg_log2FC)) #sort in descending order by log2FC
View(EarEpi)

EarEpi100 <-top_n(EarEpi, 100, wt=avg_log2FC) #Take only top 100 - upregulated genes
View(EarEpi100)

#Epithelial - late
StatusDiff <- FindMarkers(EpiOnly, ident.1="lateSEC", ident.2="control")

LEpi <- StatusDiff %>%
  arrange(desc(avg_log2FC)) #sort in descending order by log2FC
View(LEpi)

LEpi100 <-top_n(LEpi, 100, wt=avg_log2FC) #Take only top 100 - upregulated genes
View(LEpi100)

