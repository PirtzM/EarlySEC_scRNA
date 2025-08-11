################################################################################
#################################Figure 3#######################################
#######Census of Epithelial Cell types in the mouse uterus-stages of SEC########
################################################################################

#R v4.1.1
#Load Libraries
library(Seurat)
library(ggplot2)
library(patchwork)

#Set working directory
setwd("/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/scripts")

####Load each SEC stage Epithelial Subsets####
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

#Load Color Palette ####
#UMAPs
mycols=c('Progenitor'='#85b7b6', #referred to as 'Progenitor-like' in manuscript
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

#Bar Chart
CCmycols=c('LE'='#2b7b7b',
           'EP'='#85b7b6',
           'GE'='#34535e',
           'CE'='#08314a',
           'UE'='#8f9779',
           'Foxj1+'='#9f2b68')

####Figure 3A - Normal Sample Epithelial Subset UMAP ####
UMAP_N = DimPlot(object=ContEpi, reduction="umap", group.by = "seurat_clusters2", 
               repel = FALSE,                     
               label =FALSE,
               pt.size = 0.75,                    
               label.size = 4,
               cols=mycols
)+ NoLegend()
UMAP_N

ggsave('./plots/N_EPI_UMAP_clusters.pdf', last_plot(), device='pdf', dpi=300)

####Figure 3B - Pre-dysplastic Sample Epithelial Subset UMAP ####
UMAP_PD = DimPlot(object=EarEpi, reduction="umap", group.by = "seurat_clusters2", 
                 repel = FALSE,                     
                 label =FALSE,
                 pt.size = 0.75,                    
                 label.size = 4,
                 cols=mycols
)+ NoLegend()
UMAP_PD

ggsave('./plots/PD_EPI_UMAP_clusters.pdf', last_plot(), device='pdf', dpi=300)

####Figure 3C - Dysplastic Sample Epithelial Subset UMAP ####
UMAP_D = DimPlot(object=LEpi, reduction="umap", group.by = "seurat_clusters2", 
                 repel = FALSE,                     
                 label =FALSE,
                 pt.size = 0.75,                    
                 label.size = 4,
                 cols=mycols
)+ NoLegend()
UMAP_N

ggsave('./plots/D_EPI_UMAP_clusters.pdf', last_plot(), device='pdf', dpi=300)

####Figure 3D - Epithelial Cell Type quantification across SEC progression####
#Frequency of Stage by CC cluster - Epi Clusters
#Normal Samples
ContTab <- table(ContEpi$RedSEC_stage)
ContTab
ContTab<-table(ContEpi$RedSEC_stage, ContEpi$seurat_clusters_CC)
ContTab

#Pre-dysplastic Samples
EarTab <- table(EarEpi$RedSEC_stage)
EarTab
EarTab<-table(EarEpi$RedSEC_stage, EarEpi$seurat_clusters_CC)
EarTab

#Dysplastic Samples
LTab <- table(LEpi$RedSEC_stage)
LTab
LTab<-table(LEpi$RedSEC_stage, LEpi$seurat_clusters_CC)
LTab

#Merge all tables
AllTab <- merge(ContTab, EarTab, all=T)
AllTab <- merge(AllTab, LTab, all=T)

#Calculate cell type frequencies
Freq<-data.frame(AllTab)
Freq.list<-split(Freq, AllTab$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab
NewTab <- filter(NewTab, Freq>0)

#Stacked Bar Plot
NewTab$Var2 <- factor(x=NewTab$Var2, levels=c('EP','LE','GE','Foxj1+','CE','UE' #reorder cell types
))

plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(x = Var1,    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = Var2  # Variable to fill the bars   
               )) +
  theme_classic() +                  # There are various plot themes you can choose from (https://ggplot2.tidyverse.org/reference/ggtheme.html)
  # Bar plot
  geom_bar(position ='fill',       # Position of bars.  Dodge means the bars are next to each other.
           stat = 'identity'        # Height of bars represent values in the data
  )  +                           # Name of plot you want to add customizations to
  labs(x = "Cell Type",                 # x-axis label
       y = "Fraction of Epithelial Cells") +
  # Size of bars
  theme(text = element_text(size = 30),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "right")+
  scale_fill_manual(values=CCmycols)+
  scale_x_discrete(limits = c('control','earlySEC','lateSEC' #(Normal, Pre-dysplastic, Dysplastic)
  ))
plot # View plot

ggsave('./plots/AllStage_EPIcombBar_final.pdf', plot, device='pdf', width=6, height=8, units='in', dpi=300)

