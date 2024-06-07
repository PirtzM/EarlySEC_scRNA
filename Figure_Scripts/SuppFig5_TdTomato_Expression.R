################################################################################
##############################Supp Figure 5#####################################
###################TdTomato positive cells per SEC stage########################
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

#Color Palette
tdtCols <- c("Pos"="darkred",
             'Neg'='#909090')


#Load epithelial subsets
ContEpi = readRDS(file ="./data/FinalSamples/cont_Epi_mU7_mU30_Final_01262024_noMeso.rds",  # Filename
                  refhook = NULL)
EarEpi = readRDS(file ="./data/Ear_Epi_mU7_mU30_recluster_11212023_withUECE.rds",  # Filename
                 refhook = NULL)
LEpi = readRDS(file ="./data/Late_Epi_mU7_mU30_recluster_01282024_CE.rds",  # Filename
               refhook = NULL)

####Supp. Fig 5A - tdTomato Expression and positive cells in Normal Samples####
#Feature Plot
Feat_Tdt <- FeaturePlot(ContEpi, features='TdTomato-UTR', pt.size=0.25)+
            scale_color_viridis(option = 'A', begin=0.9, end=0,na.value=gray(0.8),
                                limits=c(10^-30,2.5))
Feat_Tdt
ggsave('./plots/Feature_ContEpi_TdTomato.pdf', last_plot(), device='pdf', width=6, height=4, units='in',dpi=300)

#Bar chart
ContTab <- table(ContEpi$seurat_clusters2)
ContTab
ContTab<-table(ContEpi$seurat_clusters2, ContEpi$TdTomato)
ContTab

Freq<-data.frame(ContTab)
Freq.list<-split(Freq, Freq$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab

#Stacked Bar Plot
plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(x = Var1,    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = Var2  # Variable to fill the bars   
               )) +
  theme_classic() +                  # There are various plot themes you can choose from (https://ggplot2.tidyverse.org/reference/ggtheme.html)
  # Bar plot
  geom_bar(position ='fill',       # Position of bars.  Dodge means the bars are next to each other.
           stat = 'identity'        # Height of bars represent values in the data
  )  +  
  labs(x = "Cell Type",                 # x-axis label
       y = "Fraction of Epithelial Cells") +
  # Size of bars
  theme(text = element_text(size = 30),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "right")+
  scale_fill_manual(values=tdtCols)+
  scale_x_discrete(limits = c('Progenitor','LE 1','LE 2','GE','Epithelial - Foxj1+',
                              'Cycling','Unidentified'
  ))
plot # View plot

ggsave('./plots/FinalSamples/Bar_ContEPI_TdTomato_final.pdf', plot, device='pdf', width=6, height=8, units='in', dpi=300)

####Supp. Fig 5B - tdTomato Expression and positive cells in Pre-dysplastic Samples####
#Feature Plot
Feat_Tdt <- FeaturePlot(EarEpi, features='TdTomato-UTR', pt.size=0.25)+
  scale_color_viridis(option = 'A', begin=0.9, end=0,na.value=gray(0.8),limits=c(10^-30,2.5))
Feat_Tdt
ggsave('./plots/FinalSamples/Feature_EarEpi_TdTomato.pdf', last_plot(), device='pdf', width=6, height=4, units='in',dpi=300)

#Bar Chart
EarTab <- table(EarEpi$seurat_clusters2)
EarTab
EarTab<-table(EarEpi$seurat_clusters2, EarEpi$TdTomato)
EarTab

Freq<-data.frame(EarTab)
Freq.list<-split(Freq, Freq$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab

#Stacked Bar Plot
plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(x = Var1,    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = Var2  # Variable to fill the bars   
               )) +
  theme_classic() +                  # There are various plot themes you can choose from (https://ggplot2.tidyverse.org/reference/ggtheme.html)
  # Bar plot
  geom_bar(position ='fill',       # Position of bars.  Dodge means the bars are next to each other.
           stat = 'identity'        # Height of bars represent values in the data
  )  +  
  labs(x = "Cell Type",                 # x-axis label
       y = "Fraction of Epithelial Cells") +
  # Size of bars
  theme(text = element_text(size = 30),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "right")+
  scale_fill_manual(values=tdtCols)+
  scale_x_discrete(limits = c('Progenitor 1','Progenitor 2','LE 1','LE 2',
                              'GE','Epithelial - Foxj1+',
                              'Cycling','Unidentified'
  ))
plot # View plot

ggsave('./plots/FinalSamples/Bar_EarEPI_TdTomato_final.pdf', plot, device='pdf', width=6, height=8, units='in', dpi=300)

####Supp. Fig 5C - tdTomato Expression and positive cells in Dysplastic Samples####
#Feature Plot
Feat_Tdt <- FeaturePlot(LEpi, features='TdTomato-UTR', pt.size=0.25)+
  scale_color_viridis(option = 'A', begin=0.9, end=0,na.value=gray(0.8),limits=c(10^-30,2.5))
Feat_Tdt
ggsave('./plots/FinalSamples/Feature_LateEpi_TdTomato.pdf', last_plot(), device='pdf', width=6, height=4, units='in',dpi=300)

#Bar Chart
LTab <- table(LEpi$seurat_clusters2)
LTab
LTab<-table(LEpi$seurat_clusters2, LEpi$TdTomato)
LTab

Freq<-data.frame(LTab)
Freq.list<-split(Freq, Freq$Var1)

#Divide by total cells per condition
for(i in 1:length(Freq.list)) {
  Freq.list[[i]]$FractionCells <- Freq.list[[i]]$Freq/sum(Freq.list[[i]]$Freq)
}

#Combine, must have same column titles
NewTab<-data.table::rbindlist(Freq.list, use.names=TRUE)
NewTab

#Stacked Bar Plot
plot <- ggplot(data = NewTab,            # Dataset to use for plot.  Needs to be a data.frame      
               aes(x = Var1,    # Variable to plot on the x-axis
                   y = FractionCells,   # Variable to plot on the y-axis
                   fill = Var2  # Variable to fill the bars   
               )) +
  theme_classic() +                  # There are various plot themes you can choose from (https://ggplot2.tidyverse.org/reference/ggtheme.html)
  # Bar plot
  geom_bar(position ='fill',       # Position of bars.  Dodge means the bars are next to each other.
           stat = 'identity'        # Height of bars represent values in the data
  )  +  
  labs(x = "Cell Type",                 # x-axis label
       y = "Fraction of Epithelial Cells") +
  # Size of bars
  theme(text = element_text(size = 30),                                        # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),    # Text color, angle, and horizontal adjustment on x-axis
        axis.text.y = element_text(color = 'black', hjust = 1),                # Text color and horizontal adjustment on y-axis
        legend.position = "right")+
  scale_fill_manual(values=tdtCols)+
  scale_x_discrete(limits = c('Progenitor 1','Progenitor 2','Progenitor 3','LE 1','LE 2',
                              'GE','Epithelial - Foxj1+',
                              'Cycling'
  ))
plot # View plot

ggsave('./plots/Bar_LateEPI_TdTomato_final.pdf', plot, device='pdf', width=6, height=8, units='in', dpi=300)

