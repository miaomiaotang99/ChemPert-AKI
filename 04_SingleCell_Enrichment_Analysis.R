
library(Seurat)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(stringr)
library(patchwork)
library(forcats)

unique(yjsl@meta.data$group)

perturb <- yjsl[, yjsl$group == "Perturbation"]
table(perturb@meta.data$group)

table(perturb@meta.data[["cellType"]])

Idents(perturb) <- perturb@meta.data$cellType
DefaultAssay(yjsl) <- "RNA"

species_use <- "Homo sapiens"


{
  epi_markers <- FindMarkers(
    object = perturb,
    ident.1 = "PT",
    ident.2 = NULL,
    only.pos = FALSE,
    min.pct = 0.10,
    logfc.threshold = 0,
    test.use = "wilcox"
  )

  epi_markers$gene <- rownames(epi_markers)

  epi_markers <- epi_markers %>%
    dplyr::filter(!is.na(avg_log2FC), !is.na(gene))

  epi_core_genes_df <- epi_markers %>%
    dplyr::filter(
      p_val_adj < 0.05,
      avg_log2FC > 0.5,
      pct.1 > 0.25,
      (pct.1 - pct.2) > 0.20
    ) %>%
    dplyr::arrange(desc(avg_log2FC))

  epi_core_genes <- epi_core_genes_df$gene
}
write.csv(epi_core_genes_df, "PT_core_genes.csv", row.names = FALSE)


{
  geneList_df <- epi_markers %>%
    dplyr::select(gene, avg_log2FC) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::arrange(desc(avg_log2FC))

  geneList <- geneList_df$avg_log2FC
  names(geneList) <- geneList_df$gene
  geneList <- sort(geneList, decreasing = TRUE)

}

## Hallmark
msig_h <- msigdbr(species = species_use, category = "H")

## KEGG
msig_kegg <- msigdbr(
  species = species_use,
  category = "C2",
  subcategory = "CP:KEGG_MEDICUS"
)

## Reactome
msig_reactome <- msigdbr(
  species = species_use,
  category = "C2",
  subcategory = "CP:REACTOME"
)


{
  gsea_hallmark <- GSEA(
    geneList = geneList,
    TERM2GENE = msig_h %>% dplyr::select(gs_name, gene_symbol),
    pvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE,
    seed = TRUE
  )

  gsea_kegg <- GSEA(
    geneList = geneList,
    TERM2GENE = msig_kegg %>% dplyr::select(gs_name, gene_symbol),
    pvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE,
    seed = TRUE
  )

  gsea_reactome <- GSEA(
    geneList = geneList,
    TERM2GENE = msig_reactome %>% dplyr::select(gs_name, gene_symbol),
    pvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE,
    seed = TRUE
  )
}

hallmark_res <- as.data.frame(gsea_hallmark@result)
kegg_res <- as.data.frame(gsea_kegg@result)
reactome_res <- as.data.frame(gsea_reactome@result)

setwd("C:/user")
write.csv(hallmark_res, "GSEA_Hallmark_all.csv", row.names = FALSE)
write.csv(kegg_res, "GSEA_KEGG_all.csv", row.names = FALSE)
write.csv(reactome_res, "GSEA_Reactome_all.csv", row.names = FALSE)


# KEGG_MEDICUS_ENV_FACTOR_E2_TO_NUCLEAR_INITIATED_ESTROGEN_SIGNALING_PATHWAY

hallmark_tgfb <- hallmark_res %>%
  dplyr::filter(grepl("XENOBIOTIC", Description, ignore.case = TRUE))

kegg_tgfb <- kegg_res %>%
  dplyr::filter(grepl("NNK_NNN_TO_RAS_ERK", Description, ignore.case = TRUE))

reactome_tgfb <- reactome_res %>%
  dplyr::filter(grepl("CELLULAR_RESPONSE_TO_CHEMICAL_STRESS", Description, ignore.case = TRUE)) %>%
  dplyr::arrange(pvalue)

write.csv(hallmark_tgfb, "GSEA_Hallmark_TGFb_related.csv", row.names = FALSE)
write.csv(kegg_tgfb, "GSEA_KEGG_TGFb_related.csv", row.names = FALSE)
write.csv(reactome_tgfb, "GSEA_Reactome_TGFb_related.csv", row.names = FALSE)



if (nrow(kegg_res) > 0 && any(grepl("NNK_NNN_TO_RAS_ERK", kegg_res$Description, ignore.case = TRUE))) {

  kegg_target <- kegg_res$Description[grep("NNK_NNN_TO_RAS_ERK", kegg_res$Description, ignore.case = TRUE)][1]

  pdf("GSEA_KEGG_TGFb_curve.pdf", width = 7, height = 5)
  print(
    gseaplot2(
      gsea_kegg,
      geneSetID = which(kegg_res$Description == kegg_target),
      title = paste0("KEGG: ", kegg_target)
    )
  )
  dev.off()
}



hallmark_targets <- c("XENOBIOTIC")

for (term in hallmark_targets) {
  if (term %in% hallmark_res$Description) {
    pdf(paste0(term, "_curve.pdf"), width = 7, height = 5)
    print(
      gseaplot2(
        gsea_hallmark,
        geneSetID = which(hallmark_res$Description == term),
        title = term
      )
    )
    dev.off()
  }
}


target_pathway <- "REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS"
if (nrow(reactome_tgfb) > 373) {
  best_term <- reactome_tgfb$Description[2]
  pdf("GSEA_Reactome_top_TGFb_related_curve.pdf", width = 7, height = 5)
  print(
    gseaplot2(
      gsea_reactome,
      geneSetID = which(reactome_res$Description == best_term),
      title = best_term
    )
  )
  dev.off()
}



kegg_tgfb_plot <- kegg_tgfb %>%
  dplyr::mutate(Database = "KEGG") %>%
  dplyr::select(Database, Description, NES, pvalue, p.adjust)

hallmark_tgfb_plot <- hallmark_tgfb %>%
  dplyr::mutate(Database = "Hallmark") %>%
  dplyr::select(Database, Description, NES, pvalue, p.adjust)

reactome_tgfb_plot <- reactome_tgfb %>%
  dplyr::mutate(Database = "Reactome") %>%
  dplyr::select(Database, Description, NES, pvalue, p.adjust)

tgfb_all_plot <- dplyr::bind_rows(
  kegg_tgfb_plot,
  hallmark_tgfb_plot,
  reactome_tgfb_plot
) %>%
  dplyr::filter(!is.na(NES), !is.na(pvalue)) %>%
  dplyr::arrange(NES) %>%
  dplyr::distinct(Database, Description, .keep_all = TRUE) %>%
  dplyr::mutate(
    label = paste0(Database, ": ", Description),
    label = factor(label, levels = label),
    sig = ifelse(pvalue < 0.05, "P < 0.05", "NS")
  )

write.csv(tgfb_all_plot, "TGFb_related_pathway_summary_pvalue.csv", row.names = FALSE)


{
  p_tgfb_summary <- ggplot(
    tgfb_all_plot,
    aes(x = NES, y = label, color = Database)
  ) +
    geom_segment(
      aes(x = 0, xend = NES, y = label, yend = label),
      linewidth = 0.8,
      alpha = 0.8
    ) +
    geom_point(
      aes(size = -log10(pvalue + 1e-300), shape = sig),
      stroke = 0.3
    ) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.5) +
    scale_color_manual(values = c(
      "KEGG" = "#4C72B0",
      "Hallmark" = "#C44E52",
      "Reactome" = "#55A868"
    )) +
    scale_shape_manual(values = c(
      "P < 0.05" = 16,
      "NS" = 1
    )) +
    theme_classic(base_size = 12) +
    labs(
      title = "TGF-beta-related pathways enriched in Epithelial cells",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      size = expression(-log[10](Pvalue))
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
}
print(p_tgfb_summary)

ggsave(
  "TGFb_related_pathway_summary_pvalue.pdf",
  p_tgfb_summary,
  width = 12,
  height = 8
)


{
  tgfb_all_plot_top20 <- tgfb_all_plot %>%
    dplyr::arrange(pvalue) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::mutate(
      label_wrap = stringr::str_wrap(as.character(label), width = 45)
    ) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(
      label_wrap = factor(label_wrap, levels = label_wrap)
    )

  p_tgfb_summary_top20 <- ggplot(
    tgfb_all_plot_top20,
    aes(x = NES, y = label_wrap, color = Database)
  ) +
    geom_segment(
      aes(x = 0, xend = NES, y = label_wrap, yend = label_wrap),
      linewidth = 0.8,
      alpha = 0.8
    ) +
    geom_point(
      aes(size = -log10(pvalue + 1e-300), shape = sig),
      stroke = 0.3
    ) +
    geom_vline(
      xintercept = 0,
      color = "grey60",
      linewidth = 0.5
    ) +
    scale_color_manual(values = c(
      "KEGG" = "#4C72B0",
      "Hallmark" = "#C44E52",
      "Reactome" = "#55A868"
    )) +
    scale_shape_manual(values = c(
      "P < 0.05" = 16,
      "NS" = 1
    )) +
    scale_size_continuous(range = c(2.5, 7)) +
    theme_classic(base_size = 12) +
    labs(
      title = "Top 20 TGF-beta-related pathways enriched in Epithelial cells",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      size = expression(-log[10](Pvalue))
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10)
    )
}
print(p_tgfb_summary_top20)

ggsave(
  "TGFb_related_pathway_summary_top20_pvalue.pdf",
  p_tgfb_summary_top20,
  width = 10,
  height = 7
)


hallmark_top20 <- hallmark_res %>%
  dplyr::arrange(pvalue) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::mutate(Database = "Hallmark") %>%
  dplyr::select(Database, Description, NES, pvalue, p.adjust)

kegg_top20 <- kegg_res %>%
  dplyr::arrange(pvalue) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::mutate(Database = "KEGG") %>%
  dplyr::select(Database, Description, NES, pvalue, p.adjust)

reactome_top20 <- reactome_res %>%
  dplyr::arrange(pvalue) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::mutate(Database = "Reactome") %>%
  dplyr::select(Database, Description, NES, pvalue, p.adjust)

all_top20 <- dplyr::bind_rows(hallmark_top20, kegg_top20, reactome_top20) %>%
  dplyr::filter(!is.na(NES), !is.na(pvalue)) %>%
  dplyr::mutate(
    label = paste0(Database, ": ", Description),
    sig = ifelse(pvalue < 0.05, "P < 0.05", "NS")
  ) %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(label = factor(label, levels = label))

write.csv(all_top20, "Top10_per_database_all.csv", row.names = FALSE)



library(ggplot2)

p_top20_all <- ggplot(all_top20, aes(x = NES, y = label, color = Database)) +
  geom_segment(aes(x = 0, xend = NES, y = label, yend = label), linewidth = 0.8, alpha = 0.8) +
  geom_point(aes(size = -log10(pvalue + 1e-300), shape = sig), stroke = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.5) +
  scale_color_manual(values = c("KEGG" = "#4C72B0", "Hallmark" = "#C44E52", "Reactome" = "#55A868")) +
  scale_shape_manual(values = c("P < 0.05" = 16, "NS" = 1)) +
  scale_size_continuous(range = c(2.5, 7)) +
  theme_classic(base_size = 12) +
  labs(
    title = "Top 20 pathways per database (by p-value)",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    size = expression(-log[10](Pvalue))
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(face = "bold")
  )

print(p_top20_all)
ggsave("Top_per_database_NES.pdf", p_top20_all, width = 15, height = 12)




{
  hallmark_top8 <- hallmark_res %>%
    dplyr::filter(!is.na(pvalue)) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::slice_head(n = 8)

  kegg_top8 <- kegg_res %>%
    dplyr::filter(!is.na(pvalue)) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::slice_head(n = 8)

  reactome_top8 <- reactome_res %>%
    dplyr::filter(!is.na(pvalue)) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::slice_head(n = 8)
}
write.csv(hallmark_top8, "Hallmark_top8_pathways_for_ridgeplot.csv", row.names = FALSE)
write.csv(kegg_top8, "KEGG_top8_pathways_for_ridgeplot.csv", row.names = FALSE)
write.csv(reactome_top8, "Reactome_top8_pathways_for_ridgeplot.csv", row.names = FALSE)

if (nrow(hallmark_top8) > 0) {
  hallmark_ids <- hallmark_top8$ID

  p_hallmark_ridge <- ridgeplot(
    gsea_hallmark,
    showCategory = hallmark_ids,
    fill = "pvalue",
    core_enrichment = TRUE,
    label_format = 45
  ) +
    labs(
      title = "Hallmark top 8 pathways in Epithelial cells",
      x = "Ranked gene logFC distribution",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10)
    )

  print(p_hallmark_ridge)
  ggsave("Hallmark_top8_ridgeplot.pdf", p_hallmark_ridge, width = 10, height = 7)
}

if (nrow(kegg_top8) > 0) {
  kegg_ids <- kegg_top8$ID

  p_kegg_ridge <- ridgeplot(
    gsea_kegg,
    showCategory = kegg_ids,
    fill = "pvalue",
    core_enrichment = TRUE,
    label_format = 45
  ) +
    labs(
      title = "KEGG top 8 pathways in Epithelial cells",
      x = "Ranked gene logFC distribution",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10)
    )

  print(p_kegg_ridge)
  ggsave("KEGG_top8_ridgeplot.pdf", p_kegg_ridge, width = 11, height = 7)
}

if (nrow(reactome_top8) > 0) {
  reactome_ids <- reactome_top8$ID

  p_reactome_ridge <- ridgeplot(
    gsea_reactome,
    showCategory = reactome_ids,
    fill = "pvalue",
    core_enrichment = TRUE,
    label_format = 50
  ) +
    labs(
      title = "Reactome top 8 pathways in Epithelial cells",
      x = "Ranked gene logFC distribution",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9)
    )

  print(p_reactome_ridge)
  ggsave("Reactome_top8_ridgeplot.pdf", p_reactome_ridge, width = 12, height = 7.5)
}


{
  hallmark_tgfb_top8 <- hallmark_tgfb %>%
    dplyr::filter(!is.na(pvalue)) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::slice_head(n = 10)

  kegg_tgfb_top8 <- kegg_tgfb %>%
    dplyr::filter(!is.na(pvalue)) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::slice_head(n = 10)

  reactome_tgfb_top8 <- reactome_tgfb %>%
    dplyr::filter(!is.na(pvalue)) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::slice_head(n = 10)
}
if (nrow(hallmark_tgfb_top8) > 0) {
  pdf("Hallmark_TGFb_top8_ridgeplot.pdf", width = 10, height = 6.5)
  print(
    ridgeplot(
      gsea_hallmark,
      showCategory = hallmark_tgfb_top8$ID,
      fill = "pvalue",
      core_enrichment = TRUE,
      label_format = 45
    ) +
      labs(title = "Hallmark TGF-beta-related top pathways")
  )
  dev.off()
}

if (nrow(kegg_tgfb_top8) > 0) {
  pdf("KEGG_TGFb_top8_ridgeplot.pdf", width = 10, height = 6.5)
  print(
    ridgeplot(
      gsea_kegg,
      showCategory = kegg_tgfb_top8$ID,
      fill = "pvalue",
      core_enrichment = TRUE,
      label_format = 45
    ) +
      labs(title = "KEGG TGF-beta-related top pathways")
  )
  dev.off()
}

if (nrow(reactome_tgfb_top8) > 0) {
  pdf("Reactome_TGFb_top8_ridgeplot.pdf", width = 12, height = 7)
  print(
    ridgeplot(
      gsea_reactome,
      showCategory = reactome_tgfb_top8$ID,
      fill = "pvalue",
      core_enrichment = TRUE,
      label_format = 50
    ) +
      labs(title = "Reactome TGF-beta-related top pathways")
  )
  dev.off()
}
