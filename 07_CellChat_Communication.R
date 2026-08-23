library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)

load("scRNA_processed_object.rda")

head(scdata@meta.data)


cellchat = createCellChat(object = scdata,
                          group.by = "cellType")



levels(cellchat@idents)
# group <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human

showDatabaseCategory(CellChatDB)


unique(CellChatDB$interaction$annotation)

# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
# # set the used database in the object
# cellchat@DB <- CellChatDB.use
# use Secreted Signaling for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use



##This step is necessary even if using the whole database
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)

{
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
}
{
  cellchat <- projectData(cellchat, PPI.human)

  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

  #cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)

  cellchat <- filterCommunication(cellchat, min.cells = 5)
}

save(cellchat,file = "cellchat_object.rda")
load("cellchat_object.rda")


df.net <- subsetCommunication(cellchat)
head(df.net)
write.csv(df.net,"df.net.csv")

#df.net1 <- subsetCommunication(cellchat,slot.name = "netP")
levels(cellchat@idents)

df.net1 <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
head(df.net1)
df.net2 <- subsetCommunication(cellchat, sources.use = c("T cells"), targets.use = c("Mast cells" ,"Mast cells"))
head(df.net2)
df.net3 <- subsetCommunication(cellchat, signaling = c("EGF"))
head(df.net3)


cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)




groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F, title.name = "Interaction weights/strength")


p3 <- netVisual_heatmap(cellchat)
p3
p4 <- netVisual_heatmap(cellchat, measure = "weight")
p4
p3 + p4

mat <- cellchat@net$weight

{
  mat <- cellchat@net$count
  par(mfrow = c(2,5), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
}

cellchat@netP$pathways
pathways.show <- c("TGFb")

levels(cellchat@idents)
vertex.receiver = c(1,2)

netVisual_aggregate(cellchat, signaling = "TGFb",
                    vertex.receiver = vertex.receiver,layout="hierarchy")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="TGFb", layout = "circle")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="TGFb", layout = "chord", vertex.size = groupSize)

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "TGFb", color.heatmap = "Reds")

levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 5,
                 targets.use = c(1,2,3), remove.isolate = FALSE)

cellchat@netP$pathways
netVisual_bubble(cellchat, sources.use = c(3,5), targets.use = c(1,2,4,6),
                 signaling = c("TGFb","SPP1"), remove.isolate = FALSE)


plotGeneExpression(cellchat, signaling = "SPP1")

