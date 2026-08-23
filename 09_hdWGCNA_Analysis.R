library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

{
  theme_set(theme_cowplot())

  set.seed(12345)

  enableWGCNAThreads(nThreads = 7)
}


load('scRNA_processed_object.rda')
DimPlot(scdata, group.by='cellType', label=TRUE)

unique(scdata@meta.data$group)

perturb <- scdata[, scdata$group == "Perturbation"]

seurat_obj<-perturb

rm(scdata)
gc()

{
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = "tutorial"
  )
}

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cellType", "seurat_clusters"),
  reduction = 'harmony',
  k = 25,
  max_shared = 10,
  ident.group = 'cellType'
)

{
  seurat_obj <- NormalizeMetacells(seurat_obj)
}

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "PT",
  group.by='cellType',
  assay = 'RNA',
  slot = 'data'
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'
)


setwd("C:/user")
plot_list <- PlotSoftPowers(seurat_obj)

P1=wrap_plots(plot_list, ncol=2)
P1

{
  power_table <- GetPowerTable(seurat_obj)
  head(power_table)
}

{
  seurat_obj <- ConstructNetwork(
    seurat_obj, soft_power=6,
    setDatExpr=FALSE,
    tom_name = 'PT'
  )
}
P2=PlotDendrogram(seurat_obj, main='PT hdWGCNA Dendrogram')
P2

seurat_obj <- ModuleEigengenes(
  seurat_obj,
  assay = "RNA"
)


seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cellType',
  group_name = 'PT'
)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "PT"
)



P3=PlotKMEs(seurat_obj, ncol=5, n_hubs = 20)
P3
ggsave("03_KMEs.pdf", P3, width = 12, height = 8)

{
  modules <- GetModules(seurat_obj)%>%subset(module!="grey")
}
hub_df <- GetHubGenes(seurat_obj, n_hubs = 100)

write.csv(hub_df,"hub_genes.csv")

saveRDS(seurat_obj, file = 'hdWGCNA_object.rds')



seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'#AddModuleScore
)


plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs',
  order=TRUE
)

p4=wrap_plots(plot_list, ncol=3)
p4
ggsave("04_hMEs_FeaturePlot.pdf", p4, width = 12, height = 8)

{
  ModuleCorrelogram(seurat_obj)
}
correlogram_data <- ModuleCorrelogram(seurat_obj)
write.csv(correlogram_data$corr, "module_correlation_matrix.csv")

pdf("05_ModuleCorrelogram.pdf", width = 8, height = 6)
ModuleCorrelogram(seurat_obj)
dev.off()

{
  # get hMEs from seurat object
  MEs <- GetMEs(seurat_obj, harmonized=TRUE)
  mods <- colnames(MEs); mods <- mods[mods != 'grey']

  # add hMEs to Seurat meta-data:
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
}
# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'cellType')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p6 <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p6
ggsave("06_Module_CellType_DotPlot.pdf", p6, width = 12, height = 8)


p <- VlnPlot(
  seurat_obj,
  features = 'PT1',
  group.by = 'cellType',
  pt.size = 0 # don't show actual data points
)

# add box-and-whisker plots on top:
p <- p + geom_boxplot(width=.25, fill='white')

# change axis labels and remove legend:
p7 <- p + xlab('') + ylab('hME') + NoLegend()


p7
ggsave("07_Module_VlnPlot.pdf", p7, width = 12, height = 8)

library(igraph)
ModuleNetworkPlot(seurat_obj)

pdf("p8.pdf", width = 15, height = 12)

HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10,
  n_other = 15,
  edge_prop = 0.55,
  mods = "all"
)

dev.off()

write.csv(
  power_table,
  "Supplementary_SoftPower_selection.csv",
  row.names = FALSE
)

write.csv(
  modules,
  "Supplementary_Module_assignment.csv",
  row.names = FALSE
)

# Hub genes
write.csv(
  hub_df,
  "Supplementary_Hub_genes_top100.csv",
  row.names = FALSE
)

MEs <- GetMEs(
  seurat_obj,
  harmonized=TRUE
)


write.csv(
  MEs,
  "Supplementary_ModuleEigengene_matrix.csv"
)

write.csv(
  correlogram_data$corr,
  "Supplementary_Module_correlation.csv"
)



