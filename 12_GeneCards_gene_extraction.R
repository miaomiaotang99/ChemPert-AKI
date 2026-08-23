############################################################
# GeneCards disease-associated gene extraction
#
# Description:
#   Extract disease-related genes from GeneCards search results
#   based on Relevance Score and gene category filtering.
#
# Input:
#   GeneCards-SearchResults.csv
#
# Output:
#   GeneCards_extracted_genes.csv
#   GeneCards_gene_symbols.csv
#
# Author:
#   Weidong Huang
#
############################################################


# ==========================
# Clear environment
# ==========================

rm(list = ls())


# ==========================
# Set working directory
# ==========================

setwd("C:/user")


# ==========================
# User configuration
# ==========================

# Minimum GeneCards Relevance Score
RELEVANCE_THRESHOLD <- 10


# Gene categories to retain
GENE_CATEGORIES <- c(
  "Protein Coding",
  "RNA Gene"
)


# Input and output files

input_file <- "GeneCards-SearchResults.csv"

output_file <- "GeneCards_extracted_genes.csv"

symbol_output_file <- "GeneCards_gene_symbols.csv"



# ==========================
# Load GeneCards results
# ==========================

cat("\nLoading GeneCards results...\n")


if (!file.exists(input_file)) {
  stop(
    "Input file not found. Please check the working directory."
  )
}


gene_cards <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "Loaded records:",
  nrow(gene_cards),
  "\n\n"
)



# ==========================
# Standardize column names
# ==========================

colnames(gene_cards) <- gsub(
  " ",
  ".",
  colnames(gene_cards)
)



# ==========================
# Filter by Relevance Score
# ==========================

cat("Filtering by Relevance Score...\n")


if (!"Relevance.score" %in% colnames(gene_cards)) {
  stop(
    "Column 'Relevance.score' was not found."
  )
}


gene_cards$Relevance.score <- as.numeric(
  gene_cards$Relevance.score
)


filtered_genes <- gene_cards[
  gene_cards$Relevance.score >= RELEVANCE_THRESHOLD,
]


cat(
  "Remaining genes after score filtering:",
  nrow(filtered_genes),
  "\n\n"
)



# ==========================
# Filter by gene category
# ==========================

cat("Filtering by gene category...\n")


if (!"Category" %in% colnames(filtered_genes)) {
  stop(
    "Column 'Category' was not found."
  )
}


filtered_genes <- filtered_genes[
  filtered_genes$Category %in% GENE_CATEGORIES,
]


cat(
  "Remaining genes after category filtering:",
  nrow(filtered_genes),
  "\n\n"
)



# ==========================
# Extract gene information
# ==========================

cat("Extracting gene information...\n")


required_columns <- c(
  "Gene.Symbol",
  "Description",
  "Category",
  "Uniprot.ID",
  "Relevance.score"
)


gene_table <- filtered_genes[
  ,
  required_columns
]


# Sort by relevance score

gene_table <- gene_table[
  order(
    -gene_table$Relevance.score
  ),
]


rownames(gene_table) <- NULL



# ==========================
# Save extracted gene table
# ==========================

write.csv(
  gene_table,
  file = output_file,
  row.names = FALSE
)


cat(
  "Gene table saved:",
  output_file,
  "\n"
)



# ==========================
# Save gene symbols only
# ==========================

gene_symbols <- data.frame(
  Gene.Symbol = gene_table$Gene.Symbol
)


write.csv(
  gene_symbols,
  file = symbol_output_file,
  row.names = FALSE
)


cat(
  "Gene symbol file saved:",
  symbol_output_file,
  "\n\n"
)



# ==========================
# Summary statistics
# ==========================

cat("================ Summary ================\n")

cat(
  "Input records:",
  nrow(gene_cards),
  "\n"
)

cat(
  "Final genes:",
  nrow(gene_table),
  "\n"
)

cat(
  "Relevance Score threshold: >=",
  RELEVANCE_THRESHOLD,
  "\n"
)


cat(
  "Categories:",
  paste(
    GENE_CATEGORIES,
    collapse = ", "
  ),
  "\n"
)


cat(
  "Score range:",
  round(
    min(gene_table$Relevance.score),
    2
  ),
  "-",
  round(
    max(gene_table$Relevance.score),
    2
  ),
  "\n"
)


cat(
  "Mean score:",
  round(
    mean(gene_table$Relevance.score),
    2
  ),
  "\n"
)


cat(
  "Median score:",
  round(
    median(gene_table$Relevance.score),
    2
  ),
  "\n"
)


cat("=========================================\n\n")



# ==========================
# Display top genes
# ==========================

cat("Top 10 genes by Relevance Score:\n")

print(
  head(
    gene_table[
      ,
      c(
        "Gene.Symbol",
        "Relevance.score",
        "Description"
      )
    ],
    10
  )
)


cat(
  "\nGene extraction completed successfully.\n"
)