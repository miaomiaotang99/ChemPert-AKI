library(Seurat)
library(tidyverse)
library(harmony)
library(clustree)
library(patchwork)
library(ggplot2)
library(cowplot)


AKI <- readRDS("seurat_raw_object.rds")

colnames(AKI@meta.data)

table(AKI$Type)


AKI <- NormalizeData(
  AKI,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


AKI <- FindVariableFeatures(
  AKI,
  selection.method = "vst",
  nfeatures = 2000
)

top15 <- head(VariableFeatures(AKI), 15)


plot1 <- VariableFeaturePlot(AKI, pt.size = 2, raster = TRUE)

plot2 <- LabelPoints(
  plot = plot1,
  points = top15,
  xnudge = 0,
  ynudge = 0
)

plot1 + plot2

## （Scaling）

AKI <- ScaleData(AKI, features = rownames(AKI))

##  PCA Dimensional reduction

AKI <- RunPCA(AKI, features = VariableFeatures(AKI))


VizDimLoadings(AKI, dims = 1:2, reduction = "pca")

DimPlot(AKI, reduction = "pca")

DimHeatmap(AKI, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(AKI, ndims = 20, reduction = "pca")


AKI <- RunHarmony(AKI, group.by.vars = "Type")

PC <- 1:10
AKI <- FindNeighbors(AKI, reduction = "harmony", dims = PC)



AKI <- FindClusters(AKI, resolution = seq(0.2, 1.2, 0.1))

clustree(AKI@meta.data, prefix = "RNA_snn_res.")



AKI <- FindClusters(AKI, resolution = 0.5)


AKI <- RunUMAP(AKI, reduction = "harmony", dims = PC)


DimPlot(AKI, reduction = "umap", label = TRUE, repel = TRUE)


AKI <- RunTSNE(AKI, reduction = "harmony", dims = PC)

DimPlot(AKI, reduction = "tsne")

table(AKI$seurat_clusters)


markers <- FindAllMarkers(
  AKI,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


write.csv(markers, file = "cluster_marker_genes.csv")

saveRDS(AKI, "seurat_processed_object.rds")
