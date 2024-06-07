################################################################################
####################################Figure 5####################################
###################CellChat Comparisons between all SEC stages##################
################################################################################

# Load the required libraries
library(CellChat)
library(patchwork)
library(Seurat)

#Set working directory
setwd("/workdir/mgp73/Studies/MouseSampleAnalysis/Diestrus_mU7_mU30_fixedDF/scripts")

#Load CellChat Objects####
Int.chat_cont <- readRDS('./data/Control_CellChat_02202024_noMeso.rds')
Int.chat_ear <- readRDS('./data/Early_CellChat_02202024_noMeso.rds')
Int.chat_late <- readRDS('./data/Late_CellChat_02202024_noMeso.rds')

#Load Color Palette ####
CCmycols=c('Epi'='#568cb5',
           'Fib'='#8c448c',
           'SM'='#f2d9f2',
           'LV'='#fff7b3',
           'BV'='#ffa07a',
           'Im'='#ffb6c1' ) #order matters

####Data preparation####
#Combine Data
object.list <- list(Normal = Int.chat_cont, PreDysplastic = Int.chat_ear, Dysplastic=Int.chat_late) #create list of CellChat objects
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix=TRUE) #merge CellChat objects based on names outlined above
cellchat #view object

####Figure 5A - barchart of number of inferred interactions####
gg1 <- compareInteractions(cellchat, show.legend = T, group = c(1,2,3),color.use=c('#000000','#FF7FB2','#87BFBF'))
gg1

ggsave('./plots/AllCompare_cellchat_NumBarplot.pdf', last_plot(), device='pdf', width=10, height=6, units='in', dpi=300)

####Figure 5B - barchart of strength of inferred interactions####
gg2 <- compareInteractions(cellchat, show.legend = T, group = c(1,2,3), measure = "weight", color.use=c('#000000','#FF7FB2','#87BFBF'))
gg2

ggsave('./plots/AllCompare_cellchat_StrengthBarplot.pdf', last_plot(), device='pdf', width=10, height=6, units='in', dpi=300)

####Figure 5C - Differential number of interactions - circle plots####
#Pre-dysplastic to Normal comparison
par(mfrow = c(1,2), xpd=TRUE)
pdf(file='./plots/EarCont_cellchat_circleplot.pdf', width=7, height=5, compress=FALSE)

netVisual_diffInteraction(cellchat, weight.scale = T, comparison=c(1,2),
                          vertex.label.cex = 2, title.name='Differential Number of Interactions - Pre-Dysplastic to Normal', color.use =CCmycols)

netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison=c(1,2),
                          vertex.label.cex = 2, title.name='Differential Strength of Interactions - Pre-Dysplastic to Normal', color.use = CCmycols)

#Dysplastic to Normal comparison
par(mfrow = c(1,2), xpd=TRUE)
pdf(file='./plots/FinalSamples/LateCont_cellchat_circleplot.pdf', width=7, height=5, compress=FALSE)

netVisual_diffInteraction(cellchat, weight.scale = T, comparison=c(1,3),
                          vertex.label.cex = 2, title.name='Differential Number of Interactions - Dysplastic to Normal', color.use =CCmycols)

netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison=c(1,3),
                          vertex.label.cex = 2, title.name='Differential Strength of Interactions - Dysplastic to Normal', color.use = CCmycols)

#Pre-dysplastic to Dysplastic comparison
par(mfrow = c(1,2), xpd=TRUE)
pdf(file='./plots/FinalSamples/LateEar_cellchat_circleplot.pdf', width=7, height=5, compress=FALSE)

netVisual_diffInteraction(cellchat, weight.scale = T, comparison=c(2,3),
                          vertex.label.cex = 2, title.name='Differential Number of Interactions - Dysplastic to Pre-Dysplastic', color.use =CCmycols)

netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison=c(2,3),
                          vertex.label.cex = 2, title.name='Differential Strength of Interactions - Dysplastic to Pre-Dysplastic', color.use = CCmycols)
