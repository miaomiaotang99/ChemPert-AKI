
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())


#
# install.packages("ggalluvial")
# install.packages('NMF')
# devtools::install_github("jokergoo/circlize")
# devtools::install_github("sqjin/CellChat")
# pak::pak("jokergoo/circlize")
# pak::pak("jokergoo/ComplexHeatmap")
# pak::pak("sqjin/CellChat")


library(data.table)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(CellChat)
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(data.table)
library(ggsci)


####################################################################################################################################################

setwd("C:/user")

perturb <- scdata[, scdata$group == "Perturbation"]
table(perturb@meta.data$group)

identity <- subset(perturb@meta.data, select = "cellType")
data.input <- GetAssayData(perturb,  layer = "data")
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "cellType")

CellChatDB <- CellChatDB.human
pdf("Figure 4_CellChatDB_categories.pdf", width = 8, height = 6)
showDatabaseCategory(CellChatDB)
dev.off()

##
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)

unique(CellChatDB$interaction$annotation)
# use Secreted Signaling for cell-cell communication analysis
# [1] "Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact"
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use # set the used database in the object


cellchat <- subsetData(cellchat)
options(future.globals.maxSize = 64 * 1024^3)  # 64 GiB
future::plan("multicore", workers = 10)


cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)


cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

cellchat <- filterCommunication(cellchat, min.cells = 3)


df.net.1 <- subsetCommunication(cellchat,slot.name = "netP")
df.net.2 <- subsetCommunication(cellchat )
##
df.net <- subsetCommunication(cellchat, sources.use = c(1,2,3), targets.use = c(4,5))
df.net <- subsetCommunication(cellchat, signaling = c("BMP", "NRG"))
colnames(df.net)
unique(df.net$pathway_name)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

head(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))


my36colors <-c('#625D9E', '#53A85F', '#E39A35', '#F3B1A0', '#B53E2B', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#F1BB72', '#C1E6F3', '#6778AE', '#91D0BE', '#E5D2DD',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#D6E7A3', '#68A180', '#3A6963',
               '#968175'
)
my_colors <- colorRampPalette(my36colors)(41)
my8colors <- c('
               '
               '
               '
               '
               '
               '
               '
)
my8colors1 <- c('#E39A35','#F3B1A0','#625D9E', '#53A85F', '#E95C59', '#57C3F3', '#F1BB72', '#AB3282', '#E4C755', '#476D87')

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions",color.use = my8colors1)

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use = my8colors1)

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}




levels(cellchat@idents)
vertex.receiver = c(1: 11)
cellchat@netP$pathways
pathways.show <- "NRG"
# [1] "SPP1"     "EGF"      "VISFATIN" "VEGF"     "NRG"      "FGF"      "BMP"      "GDF"      "KIT"      "PARs"     "GRN"
# [12] "SEMA3"    "GAS"      "CALCR"    "PROS"

?netVisual_aggregate
vertex.receiver = seq(1:11)
netVisual_aggregate(cellchat, signaling = "SPP1",
                    vertex.receiver = c(2,3),layout="hierarchy")
在层次图中，实体圆和空心圆分别表示源和目标。圆的大小与每个细胞组的细胞数成比例。线越粗，互作信号越强。
左图中间的target是我们选定的靶细胞。右图是选中的靶细胞之外的另外一组放在中间看互作。

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="BMP", layout = "circle")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="BMP", vertex.receiver = c(1 ),
                    layout = "chord", vertex.size = groupSize)

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "SPP1", color.heatmap = "Reds")

pathways.show <- "SPP1"
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.MK <- extractEnrichedLR(cellchat, signaling = "SPP1", geneLR.return = FALSE)
LR.show <- pairLR.MK[2,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,2) # a numeric vector
netVisual_individual(cellchat, signaling ="SPP1"  ,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout="hierarchy")
netVisual_individual(cellchat, signaling ="SPP1"  , pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling ="SPP1"  , pairLR.use = LR.show, layout = "chord")


pathway.show.all=cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver=c(1:11)
getwd()
for (i in 1:length(pathway.show.all)) {

  netVisual(cellchat,signaling = pathway.show.all[i],out.format = c("pdf"),
            vertex.receiver=vertex.receiver,layout="circle")
  plot=netAnalysis_contribution(cellchat,signaling = pathway.show.all[i])
  ggsave(filename = paste0(pathway.show.all[i],".contribution.pdf"),
         plot=plot,width=6,height=4,dpi=300,units="in")

}

#####################################################################################
levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:11), remove.isolate = FALSE)

#############################  Figure 6A ###########################
netVisual_bubble(cellchat, sources.use =c(1:11), targets.use = c(2,11), remove.isolate = FALSE)
bubble_plot <- netVisual_bubble(
  cellchat,
  sources.use = c(1:11),
  targets.use = c(2,11),
  remove.isolate = FALSE
)

class(bubble_plot)

bubble_plot_rotated <- bubble_plot +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

print(bubble_plot_rotated)

ggsave("Figure6 A.pdf", plot = bubble_plot_rotated, width = 15, height = 6, dpi = 300)



cellchat@netP$pathways
netVisual_bubble(cellchat, sources.use =c(1,11), targets.use =c(1:11),
                 signaling =  c("SPP1","EGF","VISFATIN","NRG","VEGF",
                                "FGF", "BMP","GDF","KIT","PARs","GRN",
                                "SEMA3","GAS","CALCR","PROS"), remove.isolate = FALSE)


pairLR  <- extractEnrichedLR(cellchat, signaling =c("NRG","PTN","PSAP","VISFATIN","IGF","PDGF", "BMP","SPP1","SEMA3","CX3C"), geneLR.return = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,3), targets.use =c(1:5),pairLR.use =pairLR , remove.isolate = FALSE)

netVisual_chord_gene(cellchat, sources.use = 3, targets.use = c(1:8), lab.cex = 0.5,legend.pos.y = 30)

plotGeneExpression(cellchat, signaling = "NRG")
plotGeneExpression(cellchat, signaling = "NRG", enriched.only = FALSE)
plotGeneExpression(cellchat, signaling = "NRG",type = "dot")

###################################################################################################
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling =  c("SPP1","EGF","VISFATIN","NRG","VEGF",
                                                           "FGF", "BMP","GDF","KIT","PARs","GRN",
                                                           "SEMA3","GAS","CALCR","PROS"),
                                  width = 8, height = 2.5, font.size = 10)


gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat,
                                         signaling =  c("SPP1","EGF","VISFATIN","NRG","VEGF",
                                                        "FGF", "BMP","GDF","KIT","PARs","GRN",
                                                        "SEMA3","GAS","CALCR","PROS"))
gg1 + gg2


ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#471d82", "#8c564b", "#e377c2", "#7f7f7f")

cell_annotation_colors <- my36colors[15:25]
heatmap_gradient <- colorRampPalette(my36colors)(100)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",
                                         color.use = cell_annotation_colors, color.heatmap = "RdPu")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                         color.use = cell_annotation_colors, color.heatmap = "RdPu")
ht1 + ht2
########################################################################################################

library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")


nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
##river plot
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

selectK(cellchat, pattern = "incoming")
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
##river plot
netAnalysis_river(cellchat, pattern = "incoming")

#################################################################################################

##reticulate::py_install(packages = 'umap-learn')

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

#############################################################################################
table(scRNA_harmony@meta.data$orig.ident )
sc.sp=SplitObject(scRNA_harmony,split.by = "orig.ident")
sc.11=scRNA_harmony[,sample(colnames(sc.sp[["sample2"]]),1000)]
sc.3=scRNA_harmony[,sample(colnames(sc.sp[["sample21"]]),1000)]



cellchat.sc11 <- createCellChat(object =sc.11@assays$RNA@data,
                                meta =sc.11@meta.data,  group.by ="celltype")
cellchat.sc3 <- createCellChat(object =sc.3@assays$RNA@data, meta =sc.3@meta.data,  group.by ="celltype")

dir.create("compare")
setwd("C:/user")

cellchat=cellchat.sc11
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc11 = cellchat
#################################
cellchat=cellchat.sc3
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc3 = cellchat
##############################################
cc.list=list(SC11=cc.sc11,SC3=cc.sc3)
cellchat=mergeCellChat(cc.list,cell.prefix = T,add.names = names(cc.list))
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "count")
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "weight")

netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat,weight.scale = T,measure = "weight")

netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat,measure = "weight")

rankNet(cellchat,mode = "comparison",stacked = T,do.stat = T)
rankNet(cellchat,mode = "comparison",stacked =F,do.stat = T)

weight.max=getMaxWeight(cc.list,attribute = c("idents","count"))
netVisual_circle(cc.list[[1]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )

netVisual_circle(cc.list[[2]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )


table(scRNA_harmony@active.ident)
s.cell=c( "Macrophage", "Tissue_stem_cells","Monocyte")
count1=cc.list[[1]]@net$count[s.cell,s.cell]
count2=cc.list[[2]]@net$count[s.cell,s.cell]

netVisual_circle(count1,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )

netVisual_circle(count2,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )



