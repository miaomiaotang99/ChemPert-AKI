library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(Seurat)
library(tidyverse)
library(dplyr)
library(cowplot)
#devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)

load("scRNA_processed_object.rda")

unique(scdata@meta.data$group)

perturb <- scdata[, scdata$group == "Perturbation"]
table(perturb@meta.data$group)

markers <- FindAllMarkers(perturb, only.pos = FALSE,
                          min.pct = 0.25,
                          logfc.threshold = 0)


save(markers,file = "marker_results.rda")
load("marker_results.rda")

my_colors <- rainbow(11)
# my_colors <- RColorBrewer::brewer.pal(11, "Spectral")

p1 <- jjVolcano(diffData = markers, palette = my_colors)
p1=jjVolcano(diffData = markers)
p1

my36colors <-c('#625D9E', '#53A85F', '#E39A35', '#F3B1A0', '#B53E2B', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#F1BB72', '#C1E6F3', '#6778AE', '#91D0BE', '#E5D2DD',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#D6E7A3', '#68A180', '#3A6963',
               '#968175'
)
my_colors <- rainbow(11)
p1 <- jjVolcano(diffData = markers, tile.col = my36colors)
p1


mygene <- c("TLR4", "STAT3", "PTGS2", "NLRP3", "CASP1",
            "NFE2L2", "KEAP1", "ULK1", "VDAC1", "TNF",
            "IL6", "IL1B", "FAS", "CASP3", "BAX")

{
  p2=jjVolcano(diffData = markers,
               myMarkers = mygene, tile.col = my36colors)
}
p2

p3=jjVolcano(diffData = markers,
             aesCol = c('aquamarine3','lightgoldenrod3'))

p3

library(RColorBrewer)
p3 <- jjVolcano(diffData = markers,
                aesCol = c('aquamarine3','lightgoldenrod3'),
                tile.col = brewer.pal(11, "Set3"))
p3


{
  p4=markerVolcano(markers = markers,
                   topn = 5,
                   labelCol = ggsci::pal_npg()(9))
}

p4
