############################################################
# AKI Disease Genes and DEHP Target Genes Intersection Analysis
#
# Description:
#   Identify overlapping genes between high-confidence
#   acute kidney injury (AKI)-associated genes and
#   di(2-ethylhexyl) phthalate (DEHP) target genes.
#
# Input:
#   1. High-confidence AKI disease genes
#   2. DEHP target genes
#
# Output:
#   - Venn diagram
#   - AKI-DEHP overlapping genes
#   - Summary statistics
#
# Author:
#   Weidong Huang
#
############################################################


# ============================================================
# 1. Environment setup and package loading
# ============================================================

rm(list = ls())


workDir <- "C:/user"


if (!dir.exists(workDir)) {
  dir.create(
    workDir,
    recursive = TRUE
  )
}


setwd(workDir)



# Create output directory

outputDir <- file.path(
  workDir,
  "Results"
)


if (!dir.exists(outputDir)) {
  dir.create(
    outputDir,
    recursive = TRUE
  )
}



# Load packages

library(dplyr)

library(ggplot2)

library(VennDiagram)

library(ggVennDiagram)

library(RColorBrewer)



message(
  "Working directory: ",
  getwd()
)



# ============================================================
# 2. Read input files
# ============================================================


message(
  "\n========== Loading input gene sets =========="
)



# ------------------------------------------------------------
# 2.1 High-confidence AKI disease genes
# ------------------------------------------------------------


disease_file <- file.path(
  workDir,
  "04_high_confidence_genes.csv"
)


if (file.exists(disease_file)) {


  disease_data <- read.csv(
    disease_file,
    stringsAsFactors = FALSE
  )


  if ("Gene" %in% colnames(disease_data)) {


    genes_AKI <- unique(
      disease_data$Gene
    )


  } else {


    genes_AKI <- unique(
      disease_data[,1]
    )

  }


  genes_AKI <- genes_AKI[
    genes_AKI != "" &
      !is.na(genes_AKI)
  ]


  message(
    sprintf(
      "AKI disease genes: %d",
      length(genes_AKI)
    )
  )


} else {


  stop(
    "AKI disease gene file not found!"
  )

}




# ------------------------------------------------------------
# 2.2 DEHP target genes
# ------------------------------------------------------------


dehp_file <- file.path(
  workDir,
  "DEHP_target_genes_filtered.csv"
)



if (file.exists(dehp_file)) {


  dehp_data <- read.csv(
    dehp_file,
    stringsAsFactors = FALSE
  )


  if ("Gene" %in% colnames(dehp_data)) {


    genes_DEHP <- unique(
      dehp_data$Gene
    )


  } else {


    genes_DEHP <- unique(
      dehp_data[,1]
    )

  }



  genes_DEHP <- genes_DEHP[
    genes_DEHP != "" &
      !is.na(genes_DEHP)
  ]



  message(
    sprintf(
      "DEHP target genes: %d",
      length(genes_DEHP)
    )
  )


} else {


  stop(
    "DEHP target gene file not found!"
  )

}




# ============================================================
# 3. Calculate gene intersection
# ============================================================


message(
  "\n========== Calculating AKI-DEHP overlap =========="
)



gene_list <- list(

  "AKI Disease Genes" = genes_AKI,

  "DEHP Target Genes" = genes_DEHP

)



# Shared genes

common_genes <- intersect(
  genes_AKI,
  genes_DEHP
)



AKI_specific <- setdiff(
  genes_AKI,
  genes_DEHP
)



DEHP_specific <- setdiff(
  genes_DEHP,
  genes_AKI
)



message(
  sprintf(
    "Common genes: %d",
    length(common_genes)
  )
)


message(
  sprintf(
    "AKI-specific genes: %d",
    length(AKI_specific)
  )
)


message(
  sprintf(
    "DEHP-specific genes: %d",
    length(DEHP_specific)
  )
)



# ============================================================
# 4. Venn diagram using ggVennDiagram
# ============================================================


message(
  "\n========== Drawing Venn diagram =========="
)



p_venn <- ggVennDiagram(
  gene_list,
  label = "count",
  label_alpha = 0,
  edge_size = 1
) +
  
  scale_fill_gradient(
    low = "white",
    high = "steelblue"
  ) +
  
  theme_void() +
  
  theme(
    legend.position = "right",
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 12
    )
  ) +
  
  labs(
    title = "AKI Disease Genes vs DEHP Target Genes",
    subtitle = paste0(
      "Intersection: ",
      length(common_genes),
      " genes"
    ),
    fill = "Gene Count"
  )



ggsave(
  filename = file.path(
    outputDir,
    "01_AKI_DEHP_Venn_Diagram.pdf"
  ),
  plot = p_venn,
  width = 8,
  height = 6
)



message(
  "ggVennDiagram saved."
)




# ============================================================
# 5. Classic Venn diagram
# ============================================================


pdf(
  file.path(
    outputDir,
    "01_AKI_DEHP_Venn_Diagram_Classic.pdf"
  ),
  width = 8,
  height = 8
)



venn_plot <- venn.diagram(

  x = gene_list,

  filename = NULL,

  fill = c(
    "#E41A1C",
    "#377EB8"
  ),

  alpha = 0.5,

  cex = 2,

  cat.cex = 1.5,

  cat.fontface = "bold",

  margin = 0.1,

  scaled = FALSE,

  euler.d = TRUE

)



grid::grid.draw(
  venn_plot
)



dev.off()



message(
  "Classic Venn diagram saved."
)



# ============================================================
# 6. Save overlapping genes
# ============================================================


message(
  "\n========== Saving results =========="
)



if (
  length(common_genes) > 0
) {


  overlap_table <- data.frame(

    Gene = sort(
      common_genes
    ),

    stringsAsFactors = FALSE

  )


  write.csv(

    overlap_table,

    file.path(
      outputDir,
      "02_AKI_DEHP_Common_Genes.csv"
    ),

    row.names = FALSE

  )

}



# ============================================================
# 7. Save summary statistics
# ============================================================


summary_table <- data.frame(

  Category = c(
    "AKI Disease Genes",
    "DEHP Target Genes",
    "Intersection",
    "AKI Specific Genes",
    "DEHP Specific Genes"
  ),


  Gene_Count = c(

    length(genes_AKI),

    length(genes_DEHP),

    length(common_genes),

    length(AKI_specific),

    length(DEHP_specific)

  )

)



write.csv(

  summary_table,

  file.path(
    outputDir,
    "03_AKI_DEHP_Intersection_Summary.csv"
  ),

  row.names = FALSE

)



# ============================================================
# 8. Final summary
# ============================================================


message(
  "\n=================================================="
)


message(
  "Analysis completed successfully."
)


message(
  "AKI genes: ",
  length(genes_AKI)
)


message(
  "DEHP target genes: ",
  length(genes_DEHP)
)


message(
  "AKI-DEHP overlapping genes: ",
  length(common_genes)
)



message(
  "\nOutput files:"
)


message(
  "01_AKI_DEHP_Venn_Diagram.pdf"
)


message(
  "01_AKI_DEHP_Venn_Diagram_Classic.pdf"
)


message(
  "02_AKI_DEHP_Common_Genes.csv"
)


message(
  "03_AKI_DEHP_Intersection_Summary.csv"
)



message(
  "=================================================="
)