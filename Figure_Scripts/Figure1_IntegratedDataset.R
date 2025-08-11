################################################################################
#################################Figure 1#######################################
#################Census of All Cell types in the mouse uterus###################
################################################################################

#R v4.1.1
#Load Libraries
library(Seurat)
library(ggplot2)

#Set working directory
setwd("/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/scripts")

# Load color palettes ####
#Cell type coded UMAP
mycols=c('Epithelium'='#507ECF',
         'Fib 1'='#5CAF23',
         'Fib 2'='#065B65',
         'Fib 3'='#5DCDAD',
         'Mesothelium'='#8F03C8',
         'Smooth Muscle'='#FA78FA',
         'Lymph'='#FAA000',
         'Vascular'='#940A1D',
         'Macrophage'='#B84B84',#F88379'
         'Lymphocyte'='#709EB5'#ffb6c1'
)

#SEC stage coded UMAP and bar chart
SECcols=c('control'='#031273',
          'earlySEC'='#EDCD44',
          'lateSEC'='#DC3E26')

# Cell type orders ####
#Bar chart
clusterord_int_SE= c('Epithelium',
                     'Mesothelium','Fib 1','Fib 2','Fib 3',
                     'Smooth Muscle','Lymph','Vascular', 'Macrophage','Lymphocyte')

#Load All Sample Dataset ####
IntData = readRDS(file ="./data/DiestrusMice_mU7_mU30_Final_mergedCellID_02122024.rds",  # Filename
                  refhook = NULL)
ncol(IntData) #37,543 

#Figure 1C - Integrated dataset, all cells coded by cell type UMAP ####
UMAP = DimPlot(object=IntData, reduction="umap", group.by = "seurat_clusters_SE",
               repel = TRUE,                    
               label = FALSE,
               pt.size = 0.05,   
               label.size = 4,
               cols=mycols,
               shuffle=TRUE
)
UMAP
ggsave('./plots/FinalSamples/UMAP_IntData_simpleEpi_clusters_06262025.pdf', last_plot(), device='pdf', width=6, height=4, units='in', dpi=300)

#Figure 1D - Integrated dataset, all cells coded by SEC stage UMAP ####
UMAP = DimPlot(object=IntData, reduction="umap", group.by = "RedSEC_stage",
               repel = TRUE,                   
               label = FALSE,
               pt.size = 0.05,                     
               label.size = 4,
               cols=SECcols,
               shuffle=TRUE
)
UMAP
ggsave('./plots/UMAP_SECintegration.pdf', last_plot(), device='pdf', width=6, height=4, units='in',dpi=300)

#Figure 1E - Integrated dataset, composition of each cell type by SEC stage bar chart ####
#Frequency of Stage within each SE cluster
IntTab <- table(IntData$RedSEC_stage)
IntTab
IntTab<-table(IntData$RedSEC_stage, IntData$seurat_clusters_SE)
IntTab

Freq<-data.frame(IntTab)
Freq.list<-split(Freq, Freq$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab

write.csv(NewTab, "./data/FinalSamples/SupportingData_paper1/AllCluster_AllCells_IntData_clustersize_SECstage_01182025.csv")

#Stacked Bar Plot
NewTab$Var2 <- factor(x=NewTab$Var2, levels=c(clusterord_int_SE
))

NewTab$Var1 <- factor(x=NewTab$Var1, levels=c('Normal','Pre-dysplastic','Dysplastic'
))


plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(x = Var2,    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = Var1  # Variable to fill the bars   
               )) +
  theme_classic() +                 
  # Bar plot
  geom_bar(position ='fill',      
           stat = 'identity'       
  )  +                       
  labs(x = "Cell Type",                 # x-axis label
       y = "Fraction of Cells") +
  # Size of bars
  theme(text = element_text(size = 30),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "none")+
  scale_fill_manual(values=SECcols)
plot # View plot

ggsave('./plots/FinalSamples/IntData_AllcombBar_allcells.pdf', last_plot(), device='pdf', width=6, height=8, units='in', dpi=300)
