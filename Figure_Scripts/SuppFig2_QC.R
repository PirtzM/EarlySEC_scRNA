################################################################################
##############Supp Figure 2 - QC of all Final samples (18 samples)##############
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
sampOrd = c("mU15","mU16",'mU20','mU21','mU26', #normal
            'mU22','mU23','mU24','mU27','mU28','mU29','mU30', #PD
            "mU07","mU13","mU14","mU17", 'mU18','mU19' #D
              )

clusterord_int_SE= c('Epithelium',
                     'Mesothelium','Fib 1','Fib 2','Fib 3',
                     'Smooth Muscle','Lymph','Vascular', 'Macrophage','Lymphocyte') #Integrated dataset cell type order

clusterordcon_SE= c('Epithelium',
                    'Mesothelium','N - Fib 1','N - Fib 2','Smooth Muscle', 
                    'Lymph','Vascular','Macrophage','T cell','NK cell') #Normal sample cell type order

clusterordear_SE= c('Epithelium',
                    'Mesothelium','PD - Fib 1','PD - Fib 2','PD - Fib 3',  'Smooth Muscle',
                    'Lymph','Vascular','Macrophage','Lymphocyte') #Pre-dysplastic sample cell type order

clusterordL_SE= c('Epithelium',
                  'Mesothelium','D - Fib 1','D - Fib 2','D - Fib 3', 
                  'Lymph','Vascular','Macrophage','Lymphocyte') #Dysplastic sample cell type order

clusterord_SE= c('LE','GE','DDP',
                 'Cycling','COX/MAL') #Integrated epithelial cell type order

clusterord_cont= c('N - LE 1','N - LE 2',
                   'N - LE 3','GE',
                   'DDP','Cycling','COX/MAL') #Normal sample epithelial cell type order

clusterord_ear= c('PD - LE 1','PD - LE 2',
                  'PD - LE 3','PD - LE 4', 'GE', 
                  'DDP','Cycling','COX/MAL') #Pre-dysplastic sample epithelial cell type order

clusterord_L= c('D - LE 1','D - LE 2','D - LE 3',
                'D - LE 4','D - LE 5','GE',
                'Cycling', 'COX/MAL') #Dysplastic sample epithelial cell type order

#Load Color Palettes #####
#samples
sampCol = c("mU15"='#031273',"mU16"="#031273",'mU20'='#031273','mU21'='#031273','mU26'='#031273', #normal
            'mU22'='#EDCD44','mU23'='#EDCD44','mU24'='#EDCD44','mU27'='#EDCD44','mU28'='#EDCD44','mU29'='#EDCD44','mU30'='#EDCD44', #pre-dysplastic
            "mU07"='#DC3E26',"mU13"='#DC3E26',"mU14"='#DC3E26',"mU17"='#DC3E26', 'mU18'='#DC3E26','mU19'='#DC3E26' #dysplastic
)

#Integrated Dataset
mycols=c('Epithelium'='#507ECF',
         'Fib 1'='#5CAF23',
         'Fib 2'='#065B65',
         'Fib 3'='#5DCDAD',
         'Mesothelium'='#8F03C8',
         'Smooth Muscle'='#FA78FA',
         'Lymph'='#FAA000',
         'Vascular'='#940A1D',
         'Macrophage'='#B84B84',
         'Lymphocyte'='#709EB5'
)
#Normal sample cell types
mycols_cont=c('Epithelium'='#507ECF',
              'N - Fib 1'='#5CAF23',
              'N - Fib 2'='#065B65',
              'Mesothelium'='#8F03C8',
              'Smooth Muscle'='#FA78FA',
              'Lymph'='#FAA000',
              'Vascular'='#940A1D',
              'Macrophage'='#B84B84',
              'T cell'='#709EB5',
              'NK cell'='#BEBEBE'
)

#Pre-dysplastic sample cell types
mycols_ear=c('Epithelium'='#507ECF',
             'PD - Fib 1'='#5CAF23',
             'PD - Fib 2'='#065B65',
             'PD - Fib 3'='#5DCDAD',
             'Mesothelium'='#8F03C8',
             'Smooth Muscle'='#FA78FA',
             'Lymph'='#FAA000',
             'Vascular'='#940A1D',
             'Macrophage'='#B84B84',
             'Lymphocyte'='#709EB5'
)

#Dysplastic sample cell types
mycols_L=c('Epithelium'='#507ECF',
           'D - Fib 1'='#5CAF23',
           'D - Fib 2'='#065B65',
           'D - Fib 3'='#5DCDAD',
           'Mesothelium'='#8F03C8',
           'Lymph'='#FAA000',
           'Vascular'='#940A1D', 
           'Macrophage'='#B84B84',
           'Lymphocyte'='#709EB5'
)

#Integrated epithelial dataset cell types
mycols_SE=c('LE'='#066DD4',
            'GE'='#8F03C8',
            'Cycling'='#E86DE3',
            'COX/MAL'='#5BB40B',
            'DDP'='#9f2b68'
)

#Normal sample epithelial subset cell types
mycols_ContEpi=c('N - LE 1'='#066DD4',
                 'N - LE 2'='#54F2E3',
                 'N - LE 3'='#5dcdad',
                 'GE'='#8F03C8',
                 'Cycling'='#E86DE3',
                 'COX/MAL'='#5BB40B',
                 'DDP'='#9f2b68'
)

#Pre-dysplastic sample epithelial subset cell types
mycols_EarEpi=c('PD - LE 1'='#066DD4',
                'PD - LE 2'='#54F2E3',
                'PD - LE 3'='#5dcdad',
                'PD - LE 4'='#F4877F', 
                'GE'='#8F03C8',
                'Cycling'='#E86DE3',
                'COX/MAL'='#5BB40B',
                'DDP'='#9f2b68'
)

#Dysplastic sample epithelial subset cell types
mycols_Lepi=c('D - LE 1'='#066DD4',
              'D - LE 2'='#54F2E3',
              'D - LE 3'='#5dcdad',
              'D - LE 4'='#F4877F',
              'D - LE 5'='#E3242B',
              'GE'='#8F03C8',
              'Cycling'='#E86DE3',
              'COX/MAL'='#5BB40B',
              'DDP'='#9f2b68'
) 

#Supp. Figs 2A-2C ####
#Load Non-QC with doublets data
All.merged <- readRDS('./data/FinalSamples_preQC_02232024.rds')
ncol(All.merged) #42,694 cells

#Supp. Fig 2A - percent mitochondrial genes per sample####
Pre_mt <- VlnPlot(All.merged, features=c('percent.mt'), group.by='Sample.ID', 
                  cols=sampCol, pt.size=0)+
  geom_hline(yintercept=25, linetype='dotted', color='red')+
  scale_x_discrete(limits=sampOrd)
Pre_mt

ggsave('./plots/FinalSamples/Vln_preprocess_mt_V2.pdf', Pre_mt, device='pdf', width=11, height=8.5, units='in', dpi=300)

#Supp. Fig 2B - RNA counts per sample####
Pre_Count <- VlnPlot(All.merged, features=c("nCount_RNA"), group.by='Sample.ID', 
                     cols=sampCol, pt.size=0)+
  geom_hline(yintercept=200, linetype='dotted', color='red')+
  scale_x_discrete(limits=sampOrd)
Pre_Count

ggsave('./plots/FinalSamples/Vln_preprocess_Counts_V2.pdf', Pre_Count, device='pdf', width=11, height=8.5, units='in', dpi=300)

#Supp. Fig 2C - RNA features per sample ####
Pre_feat <- VlnPlot(All.merged, features=c("nFeature_RNA"), group.by='Sample.ID',
                    cols=sampCol, pt.size=0)+
  geom_hline(yintercept=25, linetype='dotted', color='red')+
  scale_x_discrete(limits=sampOrd)
Pre_feat

ggsave('./plots/FinalSamples/Vln_preprocess_features_V2.pdf', Pre_feat, device='pdf', width=11, height=8.5, units='in', dpi=300)

#Supp. Figs 2D-2E ####
#Load Non-QC with doublets data
HQ_Cells <- readRDS('./data/FinalSamples_preDF_02232024.rds')
ncol(HQ_Cells) #38,298 cells

#Supp. Fig 2D - UMAP visualizing doublets in each cluster####
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

ggsave('./plots/FinalSamples/UMAP_preDF_doublets_paper_V2.pdf', DFUMAPFinder, device='pdf', width=11, height=8.5, units='in', dpi=300)

#Supp. Fig 2E - bar chart of doublets per cluster####
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

#Supp. Fig 2F - bar chart of cell type composition of the integrated dataset####
#load dataset
IntData = readRDS(file ="./data/DiestrusMice_mU7_mU30_Final_01182025_simple.rds",  # Filename
                  refhook = NULL)
ncol(IntData) #37,543 Cells

#create dataframe with Sample ID and clusters
tab1 <- table(IntData$Sample.ID)
tab1
Freq<-table(IntData$Sample.ID, IntData$seurat_clusters_SE)
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

write.csv(NewTab,'./data/FinalSamples/SupportingData_paper1/IntDataAll_ClusterPerSampleID_SimpleEpi_quantification_01182025.csv')

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterord_int_SE)

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
  scale_x_discrete(limits = c(sampOrd
  )) +  # Order of x axis
  theme(text = element_text(size = 15),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "right")+
  scale_fill_manual(values=mycols)
plot # View plot

ggsave('./plots/FinalSamples/Bar_IntDataAll_Samples_cluster_01182025.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)

#Supp. Fig 2G - bar chart of cell type composition of the integrated dataset epithelial cells####
#load dataset
EpiOnly = readRDS(file ="./data/allDiestrus_Epi_mU7_mU30_recluster_final_01182025.rds",
                  refhook = NULL)
ncol(EpiOnly) #19,449 cells

#create dataframe with Sample ID and clusters
tab1 <- table(EpiOnly$Sample.ID)
tab1
Freq<-table(EpiOnly$Sample.ID, EpiOnly$seurat_clusters_SE)
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

NewTab <- filter(NewTab, Freq>0)
write.csv(NewTab,'./data/FinalSamples/SupportingData_paper1/IntEpi_ClusterPerSampleID_SE_quantification_01182025.csv')

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterord_SE)

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
  scale_x_discrete(limits = c(sampOrd
  )) +  # Order of x axis
  theme(text = element_text(size = 15),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "none")+
  scale_fill_manual(values=mycols_SE)
plot # View plot

ggsave('./plots/FinalSamples/Bar_IntEpi_Clust_SE_fillSamp_01182025_NL.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)

#Supp. Fig 2H - bar chart of cell type composition of Normal Samples####
#Load samples
ContOnly = readRDS(file ="./data/cont_mU7_mU30_Final_01182025.rds",  # Filename
                  refhook = NULL)
ncol(ContOnly) #7614 cells

#create dataframe with Sample ID and clusters
tab1 <- table(ContOnly$Sample.ID)
tab1
Freq<-table(ContOnly$Sample.ID, ContOnly$seurat_clusters_SE)
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

write.csv(NewTab,'./data/FinalSamples/SupportingData_paper1/ContAll_ClusterPerSampleID_SimpleEpi_quantification_01182025.csv')

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterordcon_SE)

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
  scale_fill_manual(values=mycols_cont)
plot # View plot

ggsave('./plots/FinalSamples/Bar_ContAll_Clust_fillSamp_SE_01182025.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)

#Supp. Fig 2I - bar chart of cell type composition of Pre-dysplastic Samples####
#Load samples
EarOnly = readRDS(file ="./data/EarCan_recluster_mU7_mU30_01182025_simple.rds",  # Filename
                  refhook = NULL)
ncol(EarOnly) #17,212 cells

#create dataframe with Sample ID and clusters
tab1 <- table(EarOnly$Sample.ID)
tab1
Freq<-table(EarOnly$Sample.ID, EarOnly$seurat_clusters_SE)
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
write.csv(NewTab,'./data/FinalSamples/SupportingData_paper1/EarAll_ClusterPerSampleID_SimpleEpi_quantification_01182025.csv')

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterordear_SE)

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
  scale_fill_manual(values=mycols_Ear)
plot # View plot

ggsave('./plots/FinalSamples/Bar_EarAll_Clust_fillSamp_SE_01182025_noL.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)

####Supp. Fig 2J - bar chart of cell type composition of Dysplastic Samples####
#Load samples
LOnly = readRDS(file ="./data/LateCan_recluster_mU7_mU30_01182025_simple.rds",  # Filename
                refhook = NULL)
ncol(LOnly) #12,717 cells

#create dataframe with Sample ID and clusters
tab1 <- table(LOnly$Sample.ID)
tab1
Freq<-table(LOnly$Sample.ID, LOnly$seurat_clusters_SE)
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
write.csv(NewTab,'./data/FinalSamples/SupportingData_paper1/LateAll_ClusterPerSampleID_SimpleEpi_quantification_01182025.csv')

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterordL_SE)

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
  scale_fill_manual(values=mycols_L)
plot # View plot

ggsave('./plots/FinalSamples/Bar_LateAll_Clust_fillSamp_SE_01182025_noL.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)

#Supp. Fig 2K - bar chart of epithelial cell type composition of Normal epithelial subset samples####
#Load samples
ContEpi = readRDS(file ='./data/FinalSamples/cont_Epi_mU7_mU30_Final_01182025.rds', # Filename
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

write.csv(NewTab,'./data/FinalSamples/SupportingData_paper1/ContEpi_ClusterPerSampleID_quantification_01182025.csv')

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterord_cont)

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
  scale_fill_manual(values=mycols_ContEpi)
plot # View plot

ggsave('./plots/FinalSamples/Bar_ContEpi_Clust_fillSamp_01102025.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)

#Supp. Fig 2L - bar chart of epithelial cell type composition of Pre-dysplastic epithelial subset samples####
#Load samples
EarEpi = readRDS(file ="./data/Ear_Epi_mU7_mU30_recluster_01182025.rds",  # Filename
                  refhook = NULL)
ncol(EarEpi) #7713 cells

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

write.csv(NewTab,'./data/FinalSamples/SupportData_paper1_EarEpi_ClusterPerSampleID_quantification_01182025.csv')

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterord_ear)

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
  scale_fill_manual(values=mycols_EarEpi)
plot # View plot

ggsave('./plots/FinalSamples/Bar_EarEpi_Clust_fillSamp_01182025_noL.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)

#Supp. Fig 2M - bar chart of epithelial cell type composition of Dysplastic epithelial subset samples####
#Load samples
LEpi = readRDS(file ="./data/Late_Epi_mU7_mU30_recluster_01182025.rds",  # Filename
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

write.csv(NewTab,'./data/FinalSamples/SupportingData_paper1/LateEpi_ClusterPerSampleID_quantification_01182025.csv')

#Stacked Bar Plot
colnames(NewTab) <- c("Var1",'CellType', 'Freq','FractionCells') #rename columns
NewTab$CellType <- factor(NewTab$CellType, levels=clusterord_L)

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
  scale_fill_manual(values=mycols_Lepi)
plot # View plot

ggsave('./plots/FinalSamples/Bar_LateEpi_Clust_fillSamp_01182025_NL.pdf', last_plot(), device='pdf', width=11, height=8.5, units='in', dpi=300)
