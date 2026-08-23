
library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(harmony)
library(plyr)
library(randomcoloR)
# devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)

scdata=readRDS("seurat_processed_object.rds")

scdata <- FindClusters(scdata, resolution = 0.1)

genes <- list(
  "Podo"   = c("NPHS1", "NPHS2", "PODXL",  "WT1"),
  "PT"     = c("LRP2", "CUBN", "ALDOB"),
  "tL"     = c("AQP1", "CLCNKA"),
  "TAL"    = c("SLC12A1", "UMOD"),
  "DCT"    = c("SLC12A3"),
  "CNT"    = c("CALB1"),
  "CD-PC"  = c("AQP2", "FXYD4", "SCNN1G"),
  "CD-IC"  = c("FOXI1", "ATP6V0D2"),
  "EC"     = c("PECAM1", "EMCN"),
  "Leuk"   = c("PTPRC"),
  "IntC"   = c("ACTA2")
)

DotPlot(object = scdata,features = genes,colors.use = c("yellow","red")
           , font.size =10)
DotPlot(object = scdata, features =genes, scale = T,group.by = "seurat_clusters", dot.scale = 5) + ##celltype_l3. ###seurat_clusters
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(10, 5, 20, 5))
ggsave("Figure7 B.pdf", plot = P4, width = 15, height = 4, dpi = 300)

DotPlot(scdata, features = genes,cols = "RdYlBu") +
  RotatedAxis()

#Manual cell type annotation
Idents(scdata) <- scdata$seurat_clusters
table(scdata@active.ident)
ann.ids <- c(
  "tL",      # 0
  "TAL",       # 1
  "PT",       # 2
  "PT",       # 3
  "CNT",
  "PT",       # 5
  "EC",       # 6
  "IntC",     # 7
  "TAL",    # 8
  "CD-IC",    # 9
  "EC",       # 10
  "Leuk",     # 11
  "DCT",      # 12
  "CD-PC",
  "PT",       # 14
  "PT",       # 15
  "Podo",     # 16
  "Podo",     # 17
  "EC",       # 18
  "PT",    # 19
  "CD-IC"
)

afidens=mapvalues(Idents(scdata), from = levels(Idents(scdata)), to = ann.ids)
Idents(scdata)=afidens
scdata$cellType=Idents(scdata)

scdata$group <- ifelse(scdata$Type %in% c("GSM6433700", "GSM6433701", "GSM6433702",
                                      "GSM6433703", "GSM6433704", "GSM6433705"),
                     "Control", "Perturbation")

scdata$group <- factor(scdata$group, levels = c("Control", "Perturbation"))



DimPlot(scdata,
        reduction = "umap", label = F, label.size = 3.5,split.by = "group")+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black",
                                    size=0.5, linetype="solid"),
        legend.position = "right")


DimPlot(scdata, reduction = "tsne", label = T, label.size = 3.5,split.by = "group")+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")

library(ggsci)
cors <- pal_igv()(12)

DimPlot(scdata,reduction = "tsne",label = F,pt.size = 0.8, group.by="cellType",split.by = "group",cols = cors)#+NoLegend()
DimPlot(scdata,reduction = "umap",label = F,pt.size = 0.8, group.by="cellType",split.by = "group",cols = cors)#+NoLegend()


clusterCornerAxes(object = scdata, reduction = 'umap',groupFacet = 'group',
                  clusterCol = "cellType",
                  addCircle = F, cicAlpha = 0.03, nbin = 100,
                  cellLabel = T,cellLabelSize = 3.5) +
  scale_color_igv() + scale_fill_igv()

my.color <- c("#AF4034", "#55967e", '#006a8e', "#6a60a9",
              "#D55E00", "#E39E3E", "#0072B2",
              "#CC79A7", '#f9a11b', "#4F86C6", "#fdc23e",
              "#e3632d")

p1 <- DimPlot(object = scdata, cols = my.color, pt.size = 0.5, split.by = "group",
              label = F, reduction = "umap", label.size = 4, repel = TRUE) +
  theme(text = element_text(family = "sans"))
p1
ggsave("Umap3.pdf", plot = p1, width = 12, height = 6, dpi = 300)

p2 <- DimPlot(object = scdata, cols = my.color, pt.size = 1, split.by = "group",
              label = F, reduction = "tsne", label.size = 4, repel = TRUE) +
  theme(text = element_text(family = "sans"))#+NoLegend()
p2
ggsave("tsne3.pdf", plot = p2, width = 12, height = 6, dpi = 300)

perturb <- scdata[, scdata$group == "Perturbation"]
table(perturb@meta.data$group)

Control <- scdata[, scdata$group == "Control"]
table(Control@meta.data$group)

genes <- list(
  "Podo"   = c("NPHS1", "NPHS2", "PODXL",  "WT1"),
  "PT"     = c("LRP2", "CUBN", "ALDOB"),
  "tL"     = c("AQP1", "CLCNKA"),
  "TAL"    = c("SLC12A1", "UMOD"),
  "DCT"    = c("SLC12A3"),
  "CNT"    = c("CALB1"),
  "CD-PC"  = c("AQP2", "FXYD4", "SCNN1G"),
  "CD-IC"  = c("FOXI1", "ATP6V0D2"),
  "EC"     = c("PECAM1", "EMCN"),
  "Leuk"   = c("PTPRC"),
  "IntC"   = c("ACTA2")
)

P4 <- DotPlot(object = perturb, features =genes, scale = T,group.by = "cellType", dot.scale = 5) + ##celltype_l3. ###seurat_clusters
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(10, 5, 20, 5))
ggsave("marker gene11.pdf", plot = P4, width = 10, height = 4, dpi = 300)

P5 <- DotPlot(object = scdata, features =genes, scale = T,group.by = "cellType", dot.scale = 5) + ##celltype_l3. ###seurat_clusters
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(10, 5, 20, 5))
ggsave("marker gene21.pdf", plot = P5, width = 10, height = 4, dpi = 300)




my36colors <-c('#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#D6E7A3', '#68A180', '#3A6963',
               '#625D9E', '#53A85F', '#E39A35', '#F3B1A0', '#B53E2B', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#F1BB72', '#C1E6F3', '#6778AE', '#91D0BE', '#E5D2DD',
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
my8colors1 <- c('#625D9E', '#53A85F', '#E95C59', '#57C3F3', '#F1BB72', '#AB3282', '#E4C755', '#476D87')

DimPlot(Control,reduction = "tsne",group.by = "cellType",cols = my36colors,
        label.size = 3.5,label = F)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,
                                    color="black", size=0.5, linetype="solid"),
        legend.position = "right")#+NoLegend()

DimPlot(perturb, reduction = "tsne",
        group.by = "cellType",
        cols = my_colors,
        label = F,
        label.size = 3.5)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,
                                    color="black",
                                    size=0.5,
                                    linetype="solid"),
        legend.position = "right")+NoLegend()

DimPlot(object = Control, cols = my36colors, pt.size = 0, group.by = "cellType",
        label = F, reduction = "umap", repel = TRUE) +
  theme(text = element_text(family = "sans"))

DimPlot(object = perturb, cols = my36colors, pt.size = 0, group.by = "cellType",
        label = F, reduction = "umap", repel = TRUE) +
  theme(text = element_text(family = "sans"))

DimPlot(object = scdata, pt.size = 0, group.by = "cellType",
        label = F, reduction = "umap", repel = TRUE) +
  theme(text = element_text(family = "sans"))

DimPlot(object = scdata, cols = my36colors, pt.size = 0, group.by = "cellType",
        label = F, reduction = "umap", repel = TRUE) +
  theme(text = element_text(family = "sans"))

