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

#Dysplastic samples, n=6
LEpi = readRDS(file ="./data/Late_Epi_mU7_mU30_recluster_01182025.rds",  # Filename
                refhook = NULL)
ncol(LEpi) #7,063 cells
Idents(LEpi) <- LEpi$seurat_clusters2 #Set Active Identity

#Load Color Palette ####
#UMAPs
#Normal Samples
mycols_ContEpi=c('N - LE 1'='#066DD4',
                 'N - LE 2'='#54F2E3',
                 'N - LE 3'='#5dcdad',
                 'GE'='#8F03C8',
                 'Cycling'='#E86DE3',
                 'COX/MAL'='#5BB40B',
                 'DDP'='#9f2b68'
)

#Pre-dysplastic Samples
mycols_EarEpi=c('PD - LE 1'='#066DD4',
                'PD - LE 2'='#54F2E3',
                'PD - LE 3'='#5dcdad',
                'PD - LE 4'='#F4877F', 
                'GE'='#8F03C8',
                'Cycling'='#E86DE3',
                'COX/MAL'='#5BB40B',
                'DDP'='#9f2b68'
)

#Dysplastic Samples
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

#SEC stage coded bar chart
SECcols=c('control'='#031273',
          'earlySEC'='#EDCD44',
          'lateSEC'='#DC3E26')

#Plotting variables order ####
#Bar chart
clusterord_group= c('LE 1','LE 2','LE 3',
                'LE 4','LE 5','GE','DDP',
                'Cycling', 'COX/MAL')

#Normal DotPlot
clusterord_cont= c('N - LE 1','N - LE 2',
                   'N - LE 3','GE',
                   'DDP','Cycling','COX/MAL')

#Pre-dysplastic DotPlot
clusterord_ear= c('PD - LE 1','PD - LE 2',
                  'PD - LE 3','PD - LE 4','GE', 
                  'DDP','Cycling','COX/MAL')

#Dysplastic DotPlot
clusterord_L= c('D - LE 1','D - LE 2','D - LE 3',
                'D - LE 4','D - LE 5','GE',
                'Cycling', 'COX/MAL')

#Figure 2A - Normal Sample Epithelial Subset UMAP ####
UMAP_N = DimPlot(object=ContEpi, reduction="umap", group.by = "seurat_clusters2", 
               repel = FALSE,                     
               label =FALSE,
               pt.size = 0.75,                    
               label.size = 4,
               cols=mycols_ContEpi
)+ NoLegend()
UMAP_N

ggsave('./plots/N_EPI_UMAP_clusters.pdf', last_plot(), device='pdf', dpi=300)

#Figure 2B - Pre-dysplastic Sample Epithelial Subset UMAP ####
UMAP_PD = DimPlot(object=EarEpi, reduction="umap", group.by = "seurat_clusters2", 
                 repel = FALSE,                     
                 label =FALSE,
                 pt.size = 0.75,                    
                 label.size = 4,
                 cols=mycols_EarEpi
)+ NoLegend()
UMAP_PD

ggsave('./plots/PD_EPI_UMAP_clusters.pdf', last_plot(), device='pdf', dpi=300)

#Figure 2C - Dysplastic Sample Epithelial Subset UMAP ####
UMAP_D = DimPlot(object=LEpi, reduction="umap", group.by = "seurat_clusters2", 
                 repel = FALSE,                     
                 label =FALSE,
                 pt.size = 0.75,                    
                 label.size = 4,
                 cols=mycols_LEpi
)+ NoLegend()
UMAP_N

ggsave('./plots/D_EPI_UMAP_clusters.pdf', last_plot(), device='pdf', dpi=300)

#Figure 3D - Epithelial Cell Type quantification across SEC progression####
#Create Frequency Tables per object ####
#Frequency of Stage by CC cluster - all clusters
ContTab <- table(ContOnly$RedSEC_stage)
ContTab
ContTab<-table(ContOnly$RedSEC_stage, ContOnly$seurat_clusters_BC)
ContTab

EarTab <- table(EarOnly$RedSEC_stage)
EarTab
EarTab<-table(EarOnly$RedSEC_stage, EarOnly$seurat_clusters_BC)
EarTab

LTab <- table(LOnly$RedSEC_stage)
LTab
LTab<-table(LOnly$RedSEC_stage, LOnly$seurat_clusters_BC)
LTab

AllTab <- merge(ContTab, EarTab, all=TRUE)
AllTab <- merge(AllTab, LTab, all=TRUE)

Freq<-data.frame(AllTab)
Freq.list<-split(Freq, Freq$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab

write.csv(NewTab, './data/FinalSamples/SupportingData_paper1/AllStageCombo_frequency_SECstage_01182025.csv')

#Stacked Bar Plot
NewTab$Var2 <- factor(x=NewTab$Var2, levels=c(clusterord_group))

plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(x = Var2,    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = Var1  # Variable to fill the bars   
               )) +
  theme_classic() +               
  # Bar plot
  geom_bar(position = 'fill',       
           stat = 'identity',       
          size = 2)  +                           # Name of plot you want to add customizations to
  labs(x = "Cell Type",                 # x-axis label
       y = "Fraction of Cells") +
  # Size of bars
  theme(text = element_text(size = 15),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
       legend.position = "none")+
  scale_fill_manual(values=SECcols))
plot # View plot

ggsave('./plots/FinalSamples/AllStage_combBar_norm_SECstage_SE.pdf', last_plot(), device='pdf', width=6, height=8, units='in', dpi=300)

#Figure 3E - Epithelial Cell expression profile DotPlots (final combining of plots in Adobe Illustrator)####
#Gene list
MiniFeatures_org = c('TdTomato-UTR', #General Epithelium
                    'Ly6a',
                    'Tacstd2','Met',
                    'Lrig1','Hoxb5',
                    'Krt13',
                    'Itga6','Sox17',
                    'Foxa2','Msx1','Prom1',
                    'Dnah12',
                    'Mki67',
                    'COX1','Malat1'
)

#Normal DotPlot
mini_EPI <- DotPlot(ContEpi, features=MiniFeatures_org, col.min=0, col.max=2.5, scale.min=0
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0, limits=c(0,2.5))+ #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_cont)+
  coord_flip()
mini_EPI

ggsave('./plots/FinalSamples/DotPlot_EpiOnly_Control_clusters_06262025.pdf', last_plot(), device='pdf',width=6, height=6,units='in', dpi=300)

#Pre-dysplastic DotPlot
mini_EPI <- DotPlot(EarEpi, features=MiniFeatures_org, col.min=0, col.max=2.5,scale.min=0
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0, limits=c(0,2.5))+ #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_ear)+
  coord_flip()
mini_EPI

ggsave('./plots/FinalSamples/DotPlot_EpiOnly_Early_clusters_06262025.pdf', last_plot(), device='pdf',width=6, height=6,units='in', dpi=300)

#Dysplastic DotPlot
mini_EPI <- DotPlot(LEpi, features=MiniFeatures_org, col.min=0, col.max=2.5, scale.min=0scale.max=80 
) +
  scale_color_viridis(option = 'A', begin=0.9, end=0, limits=c(0,2.5))+ #change to magma color scale
  theme(axis.text.y=element_text(vjust=0.5, face='italic'))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_y_discrete(limits=clusterord_L)+
  coord_flip()
mini_EPI

ggsave('./plots/FinalSamples/DotPlot_EpiOnly_Late_clusters_06262025.pdf', last_plot(), device='pdf',width=6, height=6,units='in', dpi=300)

