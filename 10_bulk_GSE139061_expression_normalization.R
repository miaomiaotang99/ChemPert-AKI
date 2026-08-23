############################################################
# Bulk RNA-seq expression preprocessing and normalization
#
# Dataset:
#   GEO accession: GSE139061
#
# Description:
#   Raw expression matrix preprocessing, log2 transformation,
#   and quantile normalization using limma.
#
# Input:
#   geneMatrix.txt       : raw expression matrix
#   s1.txt               : control sample list
#   s2.txt               : treatment sample list
#
# Output:
#   GSE139061_normalized_expression.txt
#
# Author:
#   Weidong Huang
#
############################################################


# ==========================
# Load packages
# ==========================

library(limma)


# ==========================
# File configuration
# ==========================

expression_file <- "geneMatrix.txt"

control_samples <- "s1.txt"

treatment_samples <- "s2.txt"


# Set working directory
working_directory <- "C:/user"

setwd(working_directory)


# ==========================
# Load expression matrix
# ==========================

expression_matrix <- read.table(
  expression_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

expression_matrix <- as.matrix(expression_matrix)

rownames(expression_matrix) <- expression_matrix[, 1]

expression_matrix <- expression_matrix[, 2:ncol(expression_matrix)]

dimnames(expression_matrix) <- list(
  rownames(expression_matrix),
  colnames(expression_matrix)
)

expression_matrix <- matrix(
  as.numeric(as.matrix(expression_matrix)),
  nrow = nrow(expression_matrix),
  dimnames = dimnames(expression_matrix)
)

# Remove duplicated gene IDs by averaging expression values
expression_matrix <- avereps(expression_matrix)



# ==========================
# Extract control samples
# ==========================

control_sample_list <- read.table(
  control_samples,
  header = FALSE,
  sep = "\t",
  check.names = FALSE
)

control_sample_names <- as.vector(
  control_sample_list[, 1]
)

control_expression <- expression_matrix[, control_sample_names]



# ==========================
# Extract treatment samples
# ==========================

treatment_sample_list <- read.table(
  treatment_samples,
  header = FALSE,
  sep = "\t",
  check.names = FALSE
)

treatment_sample_names <- as.vector(
  treatment_sample_list[, 1]
)

treatment_expression <- expression_matrix[, treatment_sample_names]



# ==========================
# Merge expression data
# ==========================

merged_expression <- cbind(
  control_expression,
  treatment_expression
)



# ==========================
# Log2 transformation check
# ==========================

quantile_values <- as.numeric(
  quantile(
    merged_expression,
    c(0, 0.25, 0.5, 0.75, 0.99, 1.0),
    na.rm = TRUE
  )
)

need_log2 <- (
  quantile_values[5] > 100 ||
    (
      quantile_values[6] - quantile_values[1] > 50 &&
        quantile_values[2] > 0
    )
)

if (need_log2) {
  
  merged_expression[merged_expression < 0] <- 0
  
  merged_expression <- log2(
    merged_expression + 1
  )
}



# ==========================
# Expression normalization
# ==========================

normalized_expression <- normalizeBetweenArrays(
  merged_expression
)



# ==========================
# Add sample group information
# ==========================

n_control <- ncol(control_expression)

n_treatment <- ncol(treatment_expression)


group_label <- c(
  rep("Control", n_control),
  rep("Treat", n_treatment)
)


normalized_output <- rbind(
  id = paste0(
    colnames(normalized_expression),
    "_",
    group_label
  ),
  normalized_expression
)



# ==========================
# Save normalized matrix
# ==========================

write.table(
  normalized_output,
  file = "GSE139061_normalized_expression.txt",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)


############################################################
# End of script
############################################################