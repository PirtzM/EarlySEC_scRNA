################################################################################
##############Supp Figure 1 - QC of all Final samples (18 samples)##############
################################################################################

#R v4.1.1
#Load Libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#Set working directory
setwd("/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/scripts")


#Set sample and cluster Order####
sampOrd <- c("mU15","mU16",'mU20','mU21','mU26', #normal
              'mU22','mU23','mU24','mU27','mU28','mU29','mU30', #PD
              "mU07","mU13","mU14","mU17", 'mU18','mU19' #D
              )

clusterordcon= c('Progenitor-like','LE 1','LE 2', "GE",'Cycling','Unidentified',
                 'Fib 1','Fib 2','Mesothelium','Smooth Muscle', 'Lymph','Vascular',
                 'Macrophage','T cell','NK cell') #Normal sample cell type order

clusterordear= c('Progenitor-like','LE','GE','Cycling','Unidentified',
                 'Fib 1','Fib 2','Fib 3', 'Mesothelium', 'Smooth Muscle',
                 'Lymph','Vascular','Macrophage','Lymphocyte') #Pre-dysplastic sample cell type order

clusterordL= c('Progenitor-like','Progenitor-like 2','LE','GE','Cycling','Unidentified',
               'Mesothelium','Fib 1','Fib 2','Fib 3', 
               'Lymph','Vascular','Macrophage','Lymphocyte') #Dysplastic sample cell type order

clusterord_contepi= c('Progenitor','LE 1','LE 2','GE','Epithelial - Foxj1+',
                      'Cycling','Unidentified') #Normal sample epithelial cell type order

clusterord_earepi= c('Progenitor 1','Progenitor 2','LE 1','LE 2','GE','Epithelial - Foxj1+',
                  'Cycling','Unidentified') #Pre-dysplastic sample epithelial cell type order

clusterord_Lepi= c('Progenitor 1','Progenitor 2','Progenitor 3','LE 1','LE 2','GE',
                'Epithelial - Foxj1+','Cycling')


#Load Color Palettes #####
#samples
sampCol <- c("mU15"='#48CAE4',"mU16"="#00B4D8",'mU20'='#0096C7','mU21'='#0077B6','mU26'='#023E8A', #normal/blues
             'mU22'='#f1b04c','mU23'='#ee9f27','mU24'='#ec9006','mU27'='#e88504','mU28'='#e27602','mU29'='#dc6601','mU30'='#d24e01', #PD
             "mU07"='#f19798',"mU13"='#e8585a',"mU14"='#e4383b',"mU17"='#bd191b', 'mU18'='#7e1012','mU19'='#5e0c0d' #D
)

#Total cell type
mycols=c('Progenitor-like'='#85b7b6',
         'Progenitor-like 1'='#85b7b6',
         'Epi Progenitor 2'='#2bb38f', #referred to as 'Progenitor-like 2' in manuscript
         'LE'='#2b7b7b',
         'LE 1'='#2b7b7b',
         'LE 2'='#005f7b',
         'GE'='#34535e',
         'Cycling'='#08314a',
         'Unidentified'='#8f9779', #referred to as 'mt High' in manuscript
         'mt High'='#8f9779',
         'Fib 1'='#d1b8e8',
         'Fib 2'='#ab8cb1',
         'Fib 3'='#722f80',
         'Mesothelium'='#d79fd7',
         'Smooth Muscle'='#f2d9f2',
         'Lymph'='#f6eca9',
         'Vascular'='#ffa07a',
         'Macrophage'='#CC8899',
         'Lymphocyte'='#FFc1CC',
         'T cell' = '#FFc1CC',
         'NK cell' = '#Fc8eac'
         
)

#Epithelial cell types
Epi_mycols=c('Progenitor'='#85b7b6', #referred to as 'Progenitor-like' in manuscript
         'Progenitor 1'='#85b7b6', #referred to as 'Progenitor-like 1' in manuscript
         'Progenitor 2'='#2bb38f', #referred to as 'Progenitor-like 2' in manuscript
         'Progenitor 3'='#5A9F68', #referred to as 'Progenitor-like 3' in manuscript
         'LE'='#2b7b7b',
         'LE 1'='#2b7b7b',
         'LE 2'='#005f7b',
         'GE'='#34535e',
         'Cycling'='#08314a',
         'Unidentified'='#8f9779', #referred to as 'mt high' in manuscript
         'Epithelial - Foxj1+'='#9f2b68'
)

#Load Non-QC with doublets data
All.merged <- readRDS('./data/FinalSamples_preQC_02232024.rds')
ncol(All.merged) #42,694 cells

####Supp. Fig 1A - percent mitochondrial genes per sample####
Pre_mt <- VlnPlot(All.merged, features=c('percent.mt'), group.by='Sample.ID', 
                  cols=sampCol, pt.size=0)+
  geom_hline(yintercept=25, linetype='dotted', color='red')+
  scale_x_discrete(limits=sampOrd)
Pre_mt

ggsave('./plots/FinalSamples/Vln_preprocess_mt_V2.pdf', Pre_mt, device='pdf', width=11, height=8.5, units='in', dpi=300)

####Supp. Fig 1B - RNA counts per sample####
Pre_Count <- VlnPlot(All.merged, features=c("nCount_RNA"), group.by='Sample.ID', 
                     cols=sampCol, pt.size=0)+
  geom_hline(yintercept=200, linetype='dotted', color='red')+
  scale_x_discrete(limits=sampOrd)
Pre_Count

ggsave('./plots/FinalSamples/Vln_preprocess_Counts_V2.pdf', Pre_Count, device='pdf', width=11, height=8.5, units='in', dpi=300)

####Supp. Fig 1C - RNA features per sample ####
Pre_feat <- VlnPlot(All.merged, features=c("nFeature_RNA"), group.by='Sample.ID',
                    cols=sampCol, pt.size=0)+
  geom_hline(yintercept=25, linetype='dotted', color='red')+
  scale_x_discrete(limits=sampOrd)
Pre_feat

ggsave('./plots/FinalSamples/Vln_preprocess_features_V2.pdf', Pre_feat, device='pdf', width=11, height=8.5, units='in', dpi=300)

#Supp. Fig 1D - UMAP visualizing doublets in each cluster####
#Load Non-QC with doublets data
HQ_Cells <- readRDS('./data/FinalSamples_preDF_02232024.rds')
ncol(HQ_Cells) #38,298 cells

#UMAP
DFUMAPFinder <- DimPlot(object = HQ_Cells,                 # Seurat object 
                        reduction = 'umap',
                        group.by = "DF.classifications",     # Labels to color the cells by ("seurat_clusters", "Age", "Time.Point)  
                        repel = TRUE,                       # Whether to repel the cluster labels
                        label = FALSE,                        # Whether to have cluster labels 
                        pt.size = 0.15,                      # Size of each dot is (0.1 is the smallest)
                        label.size = 5,
                        shuffle=TRUE,
                        cols=c("Singlet"='lightgrey','Doublet'='#FF007F')) 
DFUMAPFinder

DotPlot(HQ_Cells, features=c('Epcam','Acta2','Foxa2', 'COX2'))

ggsave('./plots/FinalSamples/UMAP_preDF_doublets_paper_V2.pdf', DFUMAPFinder, device='pdf', width=11, height=8.5, units='in', dpi=300)

####Supp. Fig 1E - barchart of doublets per cluster####
#Doublets per Cluster
tab1 <- table(HQ_Cells$seurat_clusters)
tab1
Freq<-table(HQ_Cells$seurat_clusters, HQ_Cells$DF.classifications)
Freq

write.csv(Freq, "./data/FinalSamples/DiestrusMice_mU7_mU30_DFClusterFreq_02232024.csv")

#Rename rownames to be rough cell type
rownames(Freq) <- c('Epithelium','Lymphatic Endothelium',
                    'Fibroblast','Epithelium','Immune','Vascular Endothelium',
                    'Fibroblast','Fibroblast','Epithelium','Epithelium','Immune',
                    'Epithelium','Epithelium','Mesothelium','Epithelium','Smooth Muscle')

Freq<-data.frame(Freq)
Freq.list<-split(Freq, Freq$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab

#Stacked Bar Plot
NewTab$Var2 <- factor(x=NewTab$Var2, levels=c('Singlet','Doublet'
))

plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(factor(Var1, levels=lvls),    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = Var2 # Variable to fill the bars   
               )) +
  theme_classic() +                  # There are various plot themes you can choose from (https://ggplot2.tidyverse.org/reference/ggtheme.html)
  # Bar plot
  geom_bar(position = 'fill',       # Position of bars.  Dodge means the bars are next to each other.
           stat = 'identity',        # Height of bars represent values in the data
           size = 2)  +                           # Name of plot you want to add customizations to
  labs(x = "Cell Type",                 # x-axis label
       y = "Fraction of Cells") +  # Size of bars
  theme(text = element_text(size = 15),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "right")+
  scale_x_discrete(limits=c('Lymphatic Endothelium',"Epithelium",'Immune','Fibroblast',
                            'Vascular Endothelium','Smooth Muscle','Mesothelium')
                  ) +
 scale_fill_manual(values=c("Singlet"='lightgrey','Doublet'='hotpink'))
plot # View plot

ggsave('./plots/Bar_DoubletCount_roughCluster_paper.pdf', plot, device='pdf', width=8, height=8, units='in', dpi=300)

####Supp. Fig 1F - barchart of cell type composition of Normal Samples####
#load dataset
ContOnly = readRDS(file ="./data/cont_mU7_mU30_Final_01262024.rds",  # Filename
                   refhook = NULL)
ncol(ContOnly) #7,614 cells

#create dataframe with Sample ID and clusters
tab1 <- table(ContOnly$Sample.ID)
tab1
Freq<-table(ContOnly$Sample.ID, ContOnly$seurat_clusters2)
Freq

Freq<-data.frame(Freq)
Freq.list<-split(Freq, Freq$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterordcon) #reorder cell types

plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(x = Var1,    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = CellType  # Variable to fill the bars   
               )) +
  theme_classic() +                  # There are various plot themes you can choose from (https://ggplot2.tidyverse.org/reference/ggtheme.html)
  # Bar plot
  geom_bar(position = 'fill',       # Position of bars.  Dodge means the bars are next to each other.
           stat = 'identity',        # Height of bars represent values in the data
           size = 2)  +                           # Name of plot you want to add customizations to
  labs(x = "Sample ID",                 # x-axis label
       y = "Fraction of Cells") +
  theme(text = element_text(size = 15),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "none")+
  scale_fill_manual(values=mycols)
plot # View plot

ggsave('./plots/Bar_ContAll_Clust_fillSamp_noLegend.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)

####Supp. Fig 1G - barchart of cell type composition of Pre-dysplastic Samples####
#load dataset
EarOnly = readRDS(file ="./data/EarCan_recluster_mU7_mU30_10282023.rds",  # Filename
                  refhook = NULL)
ncol(EarOnly) #17,212 cells

#create dataframe with Sample ID and clusters
tab1 <- table(EarOnly$Sample.ID)
tab1
Freq<-table(EarOnly$Sample.ID, EarOnly$seurat_clusters2)
Freq

Freq<-data.frame(Freq)
Freq.list<-split(Freq, Freq$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterordear) #reorder cell types

plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(x = Var1,    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = CellType  # Variable to fill the bars   
               )) +
  theme_classic() +                  # There are various plot themes you can choose from (https://ggplot2.tidyverse.org/reference/ggtheme.html)
  # Bar plot
  geom_bar(position = 'fill',       # Position of bars.  Dodge means the bars are next to each other.
           stat = 'identity',        # Height of bars represent values in the data
           size = 2)  +                           # Name of plot you want to add customizations to
  labs(x = "Sample ID",                 # x-axis label
       y = "Fraction of Cells") +
  theme(text = element_text(size = 15),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "none")+
  scale_fill_manual(values=mycols)
plot # View plot

ggsave('./plots/Bar_EarAll_Clust_fillSamp_noLegend.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)

####Supp. Fig 1H - barchart of cell type composition of Dysplastic Samples####
#Load samples
LOnly = readRDS(file ="./data/LateCan_recluster_mU7_mU30_10292023.rds",  # Filename
                refhook = NULL)
ncol(LOnly) #12,717 cells

#create dataframe with Sample ID and clusters
tab1 <- table(LOnly$Sample.ID)
tab1
Freq<-table(LOnly$Sample.ID, LOnly$seurat_clusters2)
Freq

Freq<-data.frame(Freq)
Freq.list<-split(Freq, Freq$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterordL)

plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(x = Var1,    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = CellType  # Variable to fill the bars   
               )) +
  theme_classic() +                  # There are various plot themes you can choose from (https://ggplot2.tidyverse.org/reference/ggtheme.html)
  # Bar plot
  geom_bar(position = 'fill',       # Position of bars.  Dodge means the bars are next to each other.
           stat = 'identity',        # Height of bars represent values in the data
           size = 2)  +                           # Name of plot you want to add customizations to
  labs(x = "Sample ID",                 # x-axis label
       y = "Fraction of Cells") +
  theme(text = element_text(size = 15),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "right")+
  scale_fill_manual(values=mycols)
plot # View plot

ggsave('./plots/Bar_LateAll_Clust_fillSamp_withLegend.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)

####Supp. Fig 1I - barchart of epithelial cell type composition of Normal epithelial subset samples####
#Load samples
ContEpi = readRDS(file ='./data/FinalSamples/cont_Epi_mU7_mU30_Final_01262024_noMeso.rds', # Filename
                  refhook = NULL)
ncol(ContEpi) #4,689 cells

#create dataframe with Sample ID and clusters
tab1 <- table(ContEpi$Sample.ID)
tab1
Freq<-table(ContEpi$Sample.ID, ContEpi$seurat_clusters2)
Freq

Freq<-data.frame(Freq)
Freq.list<-split(Freq, Freq$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterord_contepi) #reorder cell types

plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(x = Var1,    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = CellType  # Variable to fill the bars   
               )) +
  theme_classic() +                  # There are various plot themes you can choose from (https://ggplot2.tidyverse.org/reference/ggtheme.html)
  # Bar plot
  geom_bar(position = 'fill',       # Position of bars.  Dodge means the bars are next to each other.
           stat = 'identity',        # Height of bars represent values in the data
           size = 2)  +                           # Name of plot you want to add customizations to
  labs(x = "Sample ID",                 # x-axis label
       y = "Fraction of Cells") +
  theme(text = element_text(size = 15),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "none")+
  scale_fill_manual(values=Epi_mycols)
plot # View plot

ggsave('./plots/Bar_ContEpi_Clust_fillSamp.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)

####Supp. Fig 1J - barchart of epithelial cell type composition of Pre-dysplastic epithelial subset samples####
#Load samples
EarEpi = readRDS(file ="./data/Ear_Epi_mU7_mU30_recluster_11212023_withUECE.rds",  # Filename
                  refhook = NULL)
ncol(EarEpi) #7,713 cells

#create dataframe with Sample ID and clusters
tab1 <- table(EarEpi$Sample.ID)
tab1
Freq<-table(EarEpi$Sample.ID, EarEpi$seurat_clusters2)
Freq

Freq<-data.frame(Freq)
Freq.list<-split(Freq, Freq$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterord_earepi) #reorder cell types

plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(x = Var1,    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = CellType  # Variable to fill the bars   
               )) +
  theme_classic() +                  # There are various plot themes you can choose from (https://ggplot2.tidyverse.org/reference/ggtheme.html)
  # Bar plot
  geom_bar(position = 'fill',       # Position of bars.  Dodge means the bars are next to each other.
           stat = 'identity',        # Height of bars represent values in the data
           size = 2)  +                           # Name of plot you want to add customizations to
  labs(x = "Sample ID",                 # x-axis label
       y = "Fraction of Cells") +
  theme(text = element_text(size = 15),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "none")+
  scale_fill_manual(values=Epi_mycols)
plot # View plot

ggsave('./plots/Bar_EarEpi_Clust_fillSamp_noLegend.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)

####Supp. Fig 1K - barchart of epithelial cell type composition of Dysplastic epithelial subset samples####
#Load samples
LEpi = readRDS(file ="./data/Late_Epi_mU7_mU30_recluster_01282024_CE.rds",  # Filename
                  refhook = NULL)
ncol(LEpi) #7,063 cells

#create dataframe with Sample ID and clusters
tab1 <- table(LEpi$Sample.ID)
tab1
Freq<-table(LEpi$Sample.ID, LEpi$seurat_clusters2)
Freq

Freq<-data.frame(Freq)
Freq.list<-split(Freq, Freq$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterord_Lepi) #reorder cell types

plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(x = Var1,    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = CellType  # Variable to fill the bars   
               )) +
  theme_classic() +                  # There are various plot themes you can choose from (https://ggplot2.tidyverse.org/reference/ggtheme.html)
  # Bar plot
  geom_bar(position = 'fill',       # Position of bars.  Dodge means the bars are next to each other.
           stat = 'identity',        # Height of bars represent values in the data
           size = 2)  +                           # Name of plot you want to add customizations to
  labs(x = "Sample ID",                 # x-axis label
       y = "Fraction of Cells") +
  theme(text = element_text(size = 15),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "none")+
  scale_fill_manual(values=Epi_mycols)
plot # View plot

ggsave('./plots/Bar_LateEpi_Clust_fillSamp_noLegend.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)
