library(scTenifoldKnk)
library(scTenifoldNet)
library(Matrix)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

set.seed(123)


perturb <- scdata[, scdata$group == "Perturbation"]
table(perturb@meta.data$cellType)

# =========================
# =========================

genes <- c("STAT3", "TP53", "AKT1")
# knockout_gene <- list(genes)
knockout_gene <- c(
  combn(genes, 1, simplify = FALSE),
  combn(genes, 2, simplify = FALSE),
  combn(genes, 3, simplify = FALSE)
)

# =========================
# =========================

ko_n_hvg          <- 2000
ko_nc_nNet        <- 3
ko_nc_nCells      <- 300
ko_pval_threshold <- 0.05

knockout_params <- list(
  knockout_genes = knockout_gene,
  n_hvg = ko_n_hvg,
  nc_nNet = ko_nc_nNet,
  nc_nCells = ko_nc_nCells,
  pval_threshold = ko_pval_threshold
)

# =========================
# Create output directory
# =========================

output_dir <- "analysis_results"
knockout_output_dir <- file.path(output_dir, "11_Gene_Knockout_Analysis")

if (!dir.exists(knockout_output_dir)) {
  dir.create(knockout_output_dir, recursive = TRUE)
  cat(sprintf("创建敲除分析目录: %s\n", knockout_output_dir))
}

# =========================
# =========================

available_knockout_genes <- knockout_params$knockout_genes[
  sapply(knockout_params$knockout_genes, function(x) all(x %in% rownames(perturb)))
]

cat("可运行的联合敲除组合:\n")
print(available_knockout_genes)

if (length(available_knockout_genes) == 0) {
  stop("没有可运行的联合敲除组合，请检查基因名是否存在于 Control 对象中。")
}

# =========================
# =========================

countMat <- tryCatch(
  {
    GetAssayData(perturb, layer = "counts")
  },
  error = function(e) {
    GetAssayData(perturb, slot = "counts")
  }
)

if (!inherits(countMat, "dgCMatrix")) {
  countMat <- as(countMat, "dgCMatrix")
}

# =========================
# =========================

scObject <- FindVariableFeatures(
  object = perturb,
  selection.method = "vst",
  nfeatures = knockout_params$n_hvg
)

hvgs <- VariableFeatures(scObject)

# =========================
# =========================

run_knockout_one_gene <- function(target_gene) {

  ko_name <- paste(target_gene, collapse = "_")

  cat(sprintf("\n

  gene_output_dir <- file.path(knockout_output_dir, paste0("Gene_", ko_name))

  if (!dir.exists(gene_output_dir)) {
    dir.create(gene_output_dir, recursive = TRUE)
  }

  tryCatch({

    selected_genes <- unique(c(target_gene, hvgs))
    selected_genes <- selected_genes[selected_genes %in% rownames(countMat)]

    data_wt <- countMat[selected_genes, , drop = FALSE]

    if (!inherits(data_wt, "dgCMatrix")) {
      data_wt <- as(data_wt, "dgCMatrix")
    }

    cat(sprintf("  表达矩阵: %d 基因 x %d 细胞\n", nrow(data_wt), ncol(data_wt)))
    cat("  构建联合敲除表达矩阵...\n")

    data_ko <- data_wt
    data_ko[target_gene, ] <- 0
    data_ko <- drop0(data_ko)

    cat("  运行 scTenifoldNet 联合敲除分析...\n")

    net_args <- list(
      X = data_wt,
      Y = data_ko
    )

    scTenifoldNet_args <- names(formals(scTenifoldNet))

    if ("nc_nNet" %in% scTenifoldNet_args) {
      net_args$nc_nNet <- knockout_params$nc_nNet
    }

    if ("nc_nCells" %in% scTenifoldNet_args) {
      net_args$nc_nCells <- knockout_params$nc_nCells
    }

    if ("nNet" %in% scTenifoldNet_args) {
      net_args$nNet <- knockout_params$nc_nNet
    }

    if ("nCells" %in% scTenifoldNet_args) {
      net_args$nCells <- knockout_params$nc_nCells
    }

    result <- do.call(scTenifoldNet, net_args)

    cat("  scTenifoldNet 分析完成!\n")

    df <- result$diffRegulation

    if (is.null(df)) {
      cat("  scTenifoldNet 输出对象名称:\n")
      print(names(result))
      stop("结果中没有 diffRegulation，请根据 names(result) 检查输出结构。")
    }

    df <- df[!df$gene %in% target_gene, ]

    outTab <- df[df$p.adj < knockout_params$pval_threshold, ]

    write.csv(
      outTab,
      file = file.path(gene_output_dir, paste0(ko_name, "_sigDiff.csv")),
      row.names = FALSE
    )

    write.csv(
      df,
      file = file.path(gene_output_dir, paste0(ko_name, "_KO_allResults.csv")),
      row.names = FALSE
    )

    cat(sprintf(
      "  发现 %d 个显著差异调控基因，阈值 p.adj < %.3f\n",
      nrow(outTab),
      knockout_params$pval_threshold
    ))

    df$log_p.adj <- -log10(df$p.adj)

    df$significant <- ifelse(
      df$p.adj < knockout_params$pval_threshold,
      "Significant",
      "Not significant"
    )

    label_genes <- subset(df, p.adj < knockout_params$pval_threshold)

    y_upper <- quantile(df$log_p.adj, 0.999, na.rm = TRUE)

    scatter_colors <- colorRampPalette(brewer.pal(12, "Set3"))(2)

    p_scatter <- ggplot(df, aes(x = Z, y = log_p.adj, color = significant)) +
      geom_point(alpha = 0.7, size = 1.5) +
      scale_color_manual(
        values = c(
          "Significant" = scatter_colors[1],
          "Not significant" = "gray70"
        )
      ) +
      geom_hline(
        yintercept = -log10(knockout_params$pval_threshold),
        linetype = "dashed",
        color = scatter_colors[1]
      ) +
      geom_text_repel(
        data = label_genes,
        aes(label = gene),
        size = 3,
        max.overlaps = 50,
        color = "black",
        fontface = "italic"
      ) +
      labs(
        title = paste0(ko_name, " Knockout"),
        x = "Z-score",
        y = "-log10(p.adj)"
      ) +
      theme_classic(base_size = 14) +
      coord_cartesian(ylim = c(0, y_upper)) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.position = "none"
      )

    pdf(
      file.path(gene_output_dir, paste0(ko_name, "_KO_scatter.pdf")),
      width = 6,
      height = 5
    )
    print(p_scatter)
    dev.off()

    cat(sprintf("  已保存散点图: %s_KO_scatter.pdf\n", ko_name))

    top_genes <- head(df[order(-df$FC), ], 20)

    bar_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nrow(top_genes))

    p_bar <- ggplot(
      top_genes,
      aes(x = reorder(gene, FC), y = FC, fill = reorder(gene, FC))
    ) +
      geom_bar(stat = "identity", alpha = 0.9) +
      scale_fill_manual(values = bar_colors) +
      coord_flip() +
      labs(
        title = paste0(
          "Top 20 Differentially Regulated Genes\n(",
          ko_name,
          " Knockout)"
        ),
        x = "Gene",
        y = "FC"
      ) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 11, face = "italic"),
        axis.text.x = element_text(size = 11),
        legend.position = "none"
      )

    pdf(
      file.path(gene_output_dir, paste0(ko_name, "_KO_barplot.pdf")),
      width = 6,
      height = 5
    )
    print(p_bar)
    dev.off()

    cat(sprintf("  已保存柱状图: %s_KO_barplot.pdf\n", ko_name))

    cat(sprintf("  联合敲除 %s 分析完成!\n", ko_name))

    return(df)

  }, error = function(e) {

    cat(sprintf("  联合敲除 %s 分析失败: %s\n", ko_name, e$message))

    return(NULL)
  })
}

# =========================
# =========================

results <- lapply(available_knockout_genes, run_knockout_one_gene)

names(results) <- sapply(
  available_knockout_genes,
  function(x) paste(x, collapse = "_")
)

cat("\n全部联合敲除分析完成!\n")
print(names(results))

cat("\n成功完成的组合:\n")
print(names(results)[!sapply(results, is.null)])

cat("\n失败的组合:\n")
print(names(results)[sapply(results, is.null)])

# Split visualization to view expression by groups (replaces FeatureHeatmap)
FeaturePlot(scdata, features = c("TP53", "AKT1","STAT3"), split.by = "group")

perturb <- scdata[, scdata$group == "Perturbation"]
table(perturb@meta.data$cellType)

######################  Figure 2F   ###############################
unique(scdata@meta.data$group)
gene_list <- c("TRPM3", "LRP2", "SLC28A1", "MSRA", "FRMD4B", "ADAMTS9-AS1", "CUBN", "AGXT2", "FUT6", "ACSM2B", "ACSF2", "DPYS", "SLC13A3", "CDH9", "NLGN1", "LINC02027", "ACSM2A", "TINAG", "PTPRD", "SLC6A13", "CDHR3", "SLC17A1", "AC087762.1")

P4 <- DotPlot(object = perturb, features =gene_list, scale = T,group.by = "cellType", dot.scale = 5) + ##celltype_l3. ###seurat_clusters
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(10, 5, 20, 5))
ggsave("Figure7 B.pdf", plot = P4, width = 15, height = 4, dpi = 300)

# DoHeatmap now shows a grouping bar, splitting the heatmap into groups or clusters. This can
# be changed with the `group.by` parameter
DoHeatmap(perturb, features = gene_list, cells = 1:5000, size = 4,
          angle = 90) + NoLegend()
DoHeatmap(perturb, features = gene_list) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033',
                       name = 'Z-score')
DoHeatmap(subset(perturb, downsample = 1000), features = gene_list, group.by ="cellType",
          angle = 45,size = 3)+
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#69afd8', mid = "#ffffbf", high = '#e483a4', name = 'Z-score')+NoLegend()

DoHeatmap(subset(perturb, downsample =500), features = gene_list, group.by ="cellType",
          angle = 45,size = 3)+
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(
    low = '
    mid = '
    high = '
    name = 'Z-score'
  )#+NoLegend()
