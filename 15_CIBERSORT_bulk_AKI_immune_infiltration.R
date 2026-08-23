############################################################
# CIBERSORT immune infiltration analysis in bulk RNA-seq
#
# Description:
#   Estimate immune cell fractions from bulk RNA-seq
#   expression profiles using the CIBERSORT algorithm.
#
# Dataset:
#   Bulk RNA-seq AKI samples
#
# Input:
#   normalize.txt       : normalized expression matrix
#   ref.txt             : CIBERSORT signature matrix
#   Gene25.CIBERSORT.R  : CIBERSORT R script
#
# Output:
#   - CIBERSORT immune fractions
#   - Differential immune infiltration between Control and AKI
#   - Visualization plots
#
# Author:
#   Weidong Huang
#
############################################################



# ============================================================
# 1. Load packages and set working directory
# ============================================================


rm(list = ls())


library(limma)
library(tidyverse)



working_directory <- "C:/user"


setwd(
  working_directory
)



expression_file <- "normalize.txt"



# ============================================================
# 2. Load expression matrix
# ============================================================


expression_matrix <- read.table(
  expression_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)



expression_matrix <- as.matrix(
  expression_matrix
)



rownames(expression_matrix) <- expression_matrix[,1]


expression_matrix <- expression_matrix[,2:ncol(expression_matrix), drop = FALSE]



expression_matrix <- matrix(
  as.numeric(
    as.matrix(expression_matrix)
  ),
  nrow = nrow(expression_matrix),
  dimnames = list(
    rownames(expression_matrix),
    colnames(expression_matrix)
  )
)



# ============================================================
# 3. Define sample groups
# ============================================================


sample_names <- colnames(
  expression_matrix
)



sample_group <- ifelse(
  grepl("_Control$", sample_names),
  "Control",
  ifelse(
    grepl("_Treat$", sample_names),
    "AKI",
    "Unknown"
  )
)



print(
  table(sample_group)
)



if(any(sample_group == "Unknown")){
  
  warning(
    "Unknown samples detected:"
  )
  
  print(
    sample_names[
      sample_group == "Unknown"
    ]
  )
  
}



sample_annotation <- data.frame(
  
  Sample = sample_names,
  
  Group = sample_group,
  
  stringsAsFactors = FALSE
  
)



write.table(
  
  sample_annotation,
  
  file = "Sample_Group.txt",
  
  sep = "\t",
  
  quote = FALSE,
  
  row.names = FALSE
  
)



# ============================================================
# 4. Expression preprocessing
# ============================================================


# Merge duplicated gene symbols

expression_matrix <- avereps(
  expression_matrix
)



# Remove genes with zero mean expression

expression_matrix <- expression_matrix[
  rowMeans(
    expression_matrix,
    na.rm = TRUE
  ) > 0,
]



# Replace missing values

expression_matrix[
  is.na(expression_matrix)
] <- 0



dim(
  expression_matrix
)



# Normalization is skipped because input data are already normalized

cibersort_input <- expression_matrix



cibersort_input <- rbind(
  
  ID = colnames(cibersort_input),
  
  cibersort_input
  
)



write.table(
  
  cibersort_input,
  
  file = "CIBERSORT_input_expression.txt",
  
  sep = "\t",
  
  quote = FALSE,
  
  col.names = FALSE
  
)



# ============================================================
# 5. Run CIBERSORT
# ============================================================


source(
  "Gene25.CIBERSORT.R"
)



cibersort_result <- CIBERSORT(
  
  "ref.txt",
  
  "CIBERSORT_input_expression.txt",
  
  perm = 1000,
  
  QN = FALSE
  
)



cibersort_result <- as.data.frame(
  cibersort_result,
  check.names = FALSE
)



# ============================================================
# 6. Add sample annotation
# ============================================================


cibersort_result$Sample <- rownames(
  cibersort_result
)



cibersort_result$Group <- ifelse(
  
  grepl(
    "_Control$",
    cibersort_result$Sample
  ),
  
  "Control",
  
  ifelse(
    grepl(
      "_Treat$",
      cibersort_result$Sample
    ),
    
    "AKI",
    
    NA
  )
)



cibersort_result <- cibersort_result[
  ,
  c(
    "Sample",
    "Group",
    setdiff(
      colnames(cibersort_result),
      c(
        "Sample",
        "Group"
      )
    )
  )
]



write.table(
  
  cibersort_result,
  
  file = "CIBERSORT_results.txt",
  
  sep = "\t",
  
  quote = FALSE,
  
  row.names = FALSE
  
)



# ============================================================
# 7. Extract immune fractions
# ============================================================


immune_statistics <- c(
  
  "Sample",
  
  "Group",
  
  "P-value",
  
  "Correlation",
  
  "RMSE"
  
)



immune_fraction <- cibersort_result[
  ,
  c(
    "Sample",
    "Group",
    setdiff(
      colnames(cibersort_result),
      immune_statistics
    )
  )
]



write.table(
  
  immune_fraction,
  
  file = "CIBERSORT_immune_fraction.txt",
  
  sep = "\t",
  
  quote = FALSE,
  
  row.names = FALSE
  
)



# Long format

immune_fraction_long <- immune_fraction %>%
  
  pivot_longer(
    
    cols = -c(
      Sample,
      Group
    ),
    
    names_to = "CellType",
    
    values_to = "Fraction"
    
  )



write.table(
  
  immune_fraction_long,
  
  file = "CIBERSORT_immune_fraction_long.txt",
  
  sep = "\t",
  
  quote = FALSE,
  
  row.names = FALSE
  
)



# ============================================================
# 8. Differential immune infiltration
# ============================================================


immune_difference <- immune_fraction_long %>%
  
  group_by(CellType) %>%
  
  summarise(
    
    Control_mean = mean(
      Fraction[
        Group == "Control"
      ],
      na.rm = TRUE
    ),
    
    AKI_mean = mean(
      Fraction[
        Group == "AKI"
      ],
      na.rm = TRUE
    ),
    
    p_value = wilcox.test(
      Fraction[
        Group == "Control"
      ],
      Fraction[
        Group == "AKI"
      ]
    )$p.value,
    
    .groups = "drop"
    
  )



immune_difference$FDR <- p.adjust(
  
  immune_difference$p_value,
  
  method = "BH"
  
)



immune_difference <- immune_difference[
  order(
    immune_difference$p_value
  ),
]



write.table(
  
  immune_difference,
  
  file = "CIBERSORT_Control_vs_AKI_difference.txt",
  
  sep = "\t",
  
  quote = FALSE,
  
  row.names = FALSE
  
)



# ============================================================
# 9. Visualization
# ============================================================


p <- ggplot(
  
  immune_fraction_long,
  
  aes(
    x = CellType,
    y = Fraction,
    fill = Group
  )
  
) +
  
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  labs(
    x = NULL,
    y = "Immune cell fraction",
    fill = NULL
  ) +
  
  theme(
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
    
  )



ggsave(
  
  "CIBERSORT_AKI_immune_infiltration_boxplot.pdf",
  
  plot = p,
  
  width = 10,
  
  height = 6
  
)



ggsave(
  
  "CIBERSORT_AKI_immune_infiltration_boxplot.png",
  
  plot = p,
  
  width = 10,
  
  height = 6,
  
  dpi = 600
  
)



message(
  "\nCIBERSORT immune infiltration analysis completed."
)