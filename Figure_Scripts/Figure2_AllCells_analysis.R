################################################################################
#################################Figure 2#######################################
######Census of Cell types in the mouse uterus at different stages of SEC#######
################################################################################


#R v4.1.1
#Load Libraries
library(Seurat)
library(ggplot2)
library(patchwork)

#Set working directory
setwd("/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/scripts")


####Load each SEC stage####
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

####Color Palettes####
#UMAP color palettes
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

#Bar chart
CCmycols=c('Fib'='#ab8cb1',
           'LV'='#f6eca9',
           'LE'='#2b7b7b',
           'EP'='#85b7b6',
           'GE'='#34535e',
           'Im'='#FFc1CC',
           'BV'='#ffa07a',
           'CE'='#08314a',
           'SM'='#f2d9f2',
           "Meso"='#d79fd7',
           'UE'='#8f9779',
           'Foxj1+'='#9f2b68')

####Figure 2A - Normal Sample Cell Census UMAP ####
UMAP_N = DimPlot(object=ContOnly, reduction="umap", group.by = "seurat_clusters2",
               repel = FALSE,                    
               label = FALSE,
               pt.size = 0.05,                     
               label.size = 4,
               cols=mycols
) + NoLegend()
UMAP_N

ggsave('./plots/N_UMAP_clusters.pdf', last_plot(), device='pdf', width=6, height=4, units='in', dpi=300)

####Figure 2B - Pre-dysplastic Sample Cell Census UMAP ####
UMAP_PD = DimPlot(object=EarOnly, reduction="umap", group.by='seurat_clusters2',             
               # cells=cellfind_ear,
               repel = FALSE,                    
               label = FALSE,
               pt.size = 0.05,                    
               label.size = 4,
               #cols=mycols
) + NoLegend()
UMAP_PD 

ggsave('./plots/PD_UMAP_clusters.pdf', last_plot(), device='pdf', width=6, height=4, units='in', dpi=300)

####Figure 2C - Dysplastic Cell Census UMAP ####
UMAP = DimPlot(object=LOnly, reduction="umap", group.by='seurat_clusters2', 
               repel = FALSE,                     
               label = FALSE,
               pt.size = 0.05,                     
               label.size = 4,
               cols=mycols
) + NoLegend()
UMAP

ggsave('./plots/D_UMAP_clusters.pdf', last_plot(), device='pdf', width=6, height=4, units='in', dpi=300)

####Figure 2D - Cell Type quantification across SEC progression####
#Create Frequency Tables per object
#Frequency of Stage by CC cluster
#Normal Samples
ContTab <- table(ContOnly$RedSEC_stage)
ContTab
ContTab<-table(ContOnly$RedSEC_stage, ContOnly$seurat_clusters_CC)
ContTab

#Pre-dysplastic samples
EarTab <- table(EarOnly$RedSEC_stage)
EarTab
EarTab<-table(EarOnly$RedSEC_stage, EarOnly$seurat_clusters_CC)
EarTab

#Dysplastic samples
LTab <- table(LOnly$RedSEC_stage)
LTab
LTab<-table(LOnly$RedSEC_stage, LOnly$seurat_clusters_CC)
LTab

#Merge all tables
AllTab <- merge(ContTab, EarTab, all=T)
AllTab <- merge(AllTab, LTab, all=T)

#Calculate cell type frequencies
Freq<-data.frame(AllTab)
Freq.list<-split(Freq, Freq$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab

#Stacked Bar Plot
NewTab$Var2 <- factor(x=NewTab$Var2, levels=c('EP','LE','GE','CE','UE','Meso','Fib',
                                              'SM','LV','BV','Im'
)) #reorder cell types

plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(x = Var1,    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = Var2  # Variable to fill the bars   
               )) +
  theme_classic() +                  # There are various plot themes you can choose from (https://ggplot2.tidyverse.org/reference/ggtheme.html)
  # Bar plot
  geom_bar(position = 'fill',       # Position of bars.  Dodge means the bars are next to each other.
           stat = 'identity',        # Height of bars represent values in the data
           size = 2)  +                           # Name of plot you want to add customizations to
  labs(x = "Cell Type",                 # x-axis label
       y = "Fraction of Cells") +
  # Size of bars
  theme(text = element_text(size = 15),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "right")+
  scale_fill_manual(values=CCmycols)+
  scale_x_discrete(limits = c('control','earlySEC','lateSEC' #(Normal, Pre-dysplastic, Dysplastic)
  ))
plot # View plot

ggsave('./plots/FinalSamples/AllStage_combBar_final.pdf', plot, device='pdf', width=6, height=8, units='in', dpi=300)
