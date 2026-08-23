
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
library(GSEABase)
library(harmony)
library(plyr)

setwd("C:/user")

################Create Seurat object from 10X expression matrix################
data_dir <- paste0(getwd(),"/data")
samples = list.files(data_dir)
dir = file.path(data_dir, samples)
afdata <- Read10X(data.dir = dir)
AKI <- CreateSeuratObject(counts = afdata,
                           project = "SeuratObject",
                           min.cells = 3,
                           min.features = 200)

afidens=mapvalues(Idents(AKI), from = levels(Idents(AKI)), to = samples)
Idents(AKI)=afidens
AKI$Type=Idents(AKI)

################Quality control and visualization################
AKI[["percent.mt"]] <- PercentageFeatureSet(AKI, pattern = "^MT-")
AKI[["percent.rb"]] <- PercentageFeatureSet(AKI, pattern = "^RP")
VlnPlot(AKI, features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
        ncol = 3,pt.size = 0)#,pt.size = 0

plot2 <- FeatureScatter(AKI, feature1 = "nCount_RNA", feature2 = "percent.rb")+ RotatedAxis()
plot2
plot3 <- FeatureScatter(AKI, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ RotatedAxis()
plot3

plot2 + plot3

sub1 <- AKI$nCount_RNA >= 1000
sub2 <- AKI$nFeature_RNA >= 200 & yjsl$nFeature_RNA <= 10000
sub3 <- AKI$percent.mt <= 20
sub4<-AKI$percent.rb<= 20
sub <-sub2 & sub3 & sub4
AKI <- AKI[, sub]

VlnPlot(AKI, features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
        ncol = 3,pt.size = 0)

saveRDS(AKI,"seurat_raw_object.rds")
