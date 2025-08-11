################################################################################
###################################Supp Fig. 8##################################
########################fGSEA Analysis of LE populations########################
################################################################################

#Load Libraries
library(data.table)
library(fgsea)
library(ggplot2)
library(Seurat)
library(tidyverse)


#Set working directory
setwd("/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/scripts")

#Load Reclustered datasets ####
#Integrated Epithelial, n=18
EpiOnly = readRDS(file ="./data/allDiestrus_Epi_mU7_mU30_recluster_final_01182025.rds",
                  refhook = NULL)
ncol(EpiOnly) #19,449 cells

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

#Subset datasets ####
#Normal epithelial dataset
Cont_LE <- subset(ContEpi, ident=c('N - LE1','N - LE2','N - LE3'))
Idents(Cont_LE) <- Cont_LE$barcode #set active identity for later subsetting

#Pre-dysplastic epithelial dataset
PD_LE <- subset(EarEpi, ident=c('PD - LE4'))
Idents(PD_LE) <- PD_LE$barcode #set active identity for later subsetting

#Pull barcodes from datasets for cross-referenceing ####
#Integrated epithelial for reference
Int_bar <- EpiOnly$barcode
write.csv(Int_bar, "./data/FinalSamples/GSEA/IntData_barcodes.csv") #save lists for later

#Control LE
cont_LE_bar <- Cont_LE$barcode
write.csv(cont_LE_bar, "./data/FinalSamples/GSEA/Normal_LE_barcodes.csv") #save lists for later

#Pre-dysplastic LE4
PD_LE4_bar <- PD_LE4_bar$barcode
write.csv(cont_LE_bar, "./data/FinalSamples/GSEA/PD_LE4_barcodes.csv") #save lists for later

#Read barcode CSVs
Int_bar <- read.csv("./data/FinalSamples/GSEA/IntData_barcodes.csv")
PD_LE4_bar <- read.csv("./data/FinalSamples/GSEA/PD_LE4_barcodes.csv")
Cont_LE_bar <- read.csv("./data/FinalSamples/GSEA/Normal_LE_barcodes.csv")

#Subset EpiOnly by barcodes to cross-match and make sure all cells are being properly compared in the Integrated dataset ####
Int_cells.use <- Int_bar$x

PD_LE4_cells.use <- PD_LE4_bar$x
PD_LE4_cells.use <- Int_cells.use[Int_cells.use %in% PD_LE4_cells.use]

Cont_LE_cells.use <- Cont_LE_bar$x
Cont_LE_cells.use <- Int_cells.use[Int_cells.use %in% Cont_LE_cells.use]


#Get DEG
Idents(EpiOnly) <- EpiOnly$barcode #Set Active Identity

LE4_PDvN <- FindMarkers(EpiOnly, ident.1=PD_LE4_cells.use, ident.2=Cont_LE_cells.use,
                        logfc.threshold=0)

#Save as CSV
write.csv(LE4_PDvN, './data/FinalSamples/GSEA/DiffExp_LE4_LE4PDvLEN.csv')

#Read csv ####
LE4_PDvN <- read.csv('./data/FinalSamples/GSEA/DiffExp_LE4_LE4PDvLEN.csv')

#Read DGE of single-cell data
LE4_PDvN$gene <- LE4_PDvN$X

#Order and filter the differentially expressed genes ####
genes_of_interest_LE4_PD <-LE4_PDvN[order(-LE4_PDvN$avg_log2FC),] #order differentially expressed genes from largest to smallest based on avg_log2FC
LE4_PD_gene_list <- genes_of_interest_LE4_PD$avg_log2FC  #extract gene names from list
names(LE4_PD_gene_list) <- genes_of_interest_LE4_PD$gene

#Download mSigDB databases: gsea-msigdb.org/gsea/msigdb/index.jsp ####
Reactome <- fgsea::gmtPathways("/workdir/mgp73/MSigDB_GeneSets/Mouse/m2.cp.reactome.v2025.1.Mm.symbols.gmt")
Regulatory <- fgsea::gmtPathways("/workdir/mgp73/MSigDB_GeneSets/Mouse/m3.all.v2025.1.Mm.symbols.gmt")
GO_BP <- fgsea::gmtPathways("/workdir/mgp73/MSigDB_GeneSets/Mouse/m5.go.bp.v2025.1.Mm.symbols.gmt")
GO_MF <- fgsea::gmtPathways("/workdir/mgp73/MSigDB_GeneSets/Mouse/m5.go.mf.v2025.1.Mm.symbols.gmt")
MPT <- fgsea::gmtPathways("/workdir/mgp73/MSigDB_GeneSets/Mouse/m5.mpt.v2025.1.Mm.symbols.gmt")
CT_sig <- fgsea::gmtPathways("/workdir/mgp73/MSigDB_GeneSets/Mouse/m8.all.v2025.1.Mm.symbols.gmt")
Hallmark <- fgsea::gmtPathways("/workdir/mgp73/MSigDB_GeneSets/Mouse/mh.all.v2025.1.Mm.symbols.gmt")

GO_gene_sets <- c(GO_BP, GO_MF, MPT)
All_gene_sets <- c(GO_gene_sets, Reactome, Regulatory, CT_sig, Hallmark)
base_gene_sets <- c(GO_gene_sets, Reactome, Regulatory, Hallmark)

# GSEA Analysis ####
GSEA_LE4_PDvN <- fgsea(pathways = base_gene_sets,
                       stats=LE4_PD_gene_list,
                       eps=0,
                       minSize=15,
                       maxSize=500)

#Make a table with top pathways ####
SigPath_LE4_PD <- GSEA_LE4_PDvN[GSEA_LE4_PDvN$padj<0.05] #filter only significant pathways

#Save as csv
df_LE4_PDvN <- apply(SigPath_LE4_PD, 2, as.character)
write.csv(df_LE4_PDvN, './data/FinalSamples/GSEA/TopPath_LE4_PDvN.csv')

#Prep for enrichment plot ####
#filter for only the top pathways enriched in PD-LE4
topPathwaysUp_LE4_PD <- GSEA_LE4_PDvN[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown_LE4_PD <- GSEA_LE4_PDvN[ES < 0][head(order(pval), n=20), pathway]
topPathways_LE4_PD <- c(topPathwaysUp_LE4_PD, rev(topPathwaysDown_LE4_PD))

#Supp. Fig 8A - GSEA Enrichment plot for PD-LE4 vs N-LE, Differentiation ####
PD_exmp <- plotEnrichment(base_gene_sets[['GOBP_EPITHELIAL_CELL_DIFFERENTIATION']], 
                          LE4_PD_gene_list) + labs(title='GOBP_EPITHELIAL_CELL_DIFFERENTIATION')

PD_exmp

ggsave('./plots/FinalSamples/GSEA/GOBP_EPITHELIAL_CELL_DIFFERENTIATION_LE4PDvLEN.pdf', last_plot(), dpi=300)

#Supp. Fig 8B - GSEA Enrichment plot for PD-LE4 vs N-LE, Development ####
PD_exmp_2 <- plotEnrichment(base_gene_sets[['GOBP_EPIDERMIS_DEVELOPMENT']], 
                            LE4_PD_gene_list) + labs(title='GOBP_EPIDERMIS_DEVELOPMENT')

PD_exmp_2

ggsave('./plots/FinalSamples/GSEA/GOBP_EPIDERMIS_DEVELOPMENT_LE4PDvLEN.pdf', last_plot(), dpi=300)
