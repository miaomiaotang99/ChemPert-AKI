############################################################
# EBI Search API disease-associated gene retrieval
#
# Description:
#   Download disease-associated genes from EBI Search API
#   and convert Ensembl IDs into gene symbols using Ensembl REST API.
#
# Example:
#   Disease: Acute kidney injury
#
# Input:
#   Disease name
#
# Output:
#   disease_associated_genes.csv
#
# Author:
#   Weidong Huang
#
############################################################


# ==========================
# Load packages
# ==========================

required_packages <- c(
  "httr",
  "jsonlite"
)


for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}


library(httr)
library(jsonlite)



# ==========================
# User configuration
# ==========================


disease_name <- "acute kidney injury"


working_directory <- "C:/user"


output_file <- file.path(
  working_directory,
  "EBI_GWAS_acute_kidney_injury_genes.csv"
)



# ==========================
# Download disease genes
# ==========================


download_ebi_genes <- function(
    disease_name,
    batch_size = 100
) {
  
  
  message(
    "\n============================================"
  )
  
  message(
    "EBI Disease Gene Retrieval"
  )
  
  message(
    "============================================\n"
  )
  
  
  message(
    "Searching disease: ",
    disease_name
  )
  
  
  api_url <- paste0(
    "https://www.ebi.ac.uk/ebisearch/ws/rest/",
    "geneDiseaseAssociations"
  )
  
  
  # Obtain total number of records
  
  count_url <- paste0(
    api_url,
    "?query=",
    URLencode(disease_name),
    "&size=1&start=0&format=json"
  )
  
  
  response <- GET(count_url)
  
  
  if (status_code(response) != 200) {
    stop(
      "Unable to connect to EBI Search API."
    )
  }
  
  
  result <- content(
    response,
    as = "text",
    encoding = "UTF-8"
  )
  
  
  result <- fromJSON(result)
  
  
  total_records <- result$hitCount
  
  
  message(
    "Total records found: ",
    total_records
  )
  
  
  if (total_records == 0) {
    return(NULL)
  }
  
  
  # Download records in batches
  
  all_records <- list()
  
  start_position <- 0
  
  
  while (
    start_position < total_records
  ) {
    
    
    message(
      sprintf(
        "Downloading %d - %d / %d",
        start_position + 1,
        min(
          start_position + batch_size,
          total_records
        ),
        total_records
      )
    )
    
    
    query_url <- paste0(
      api_url,
      "?query=",
      URLencode(disease_name),
      "&size=",
      batch_size,
      "&start=",
      start_position,
      "&format=json"
    )
    
    
    response <- GET(query_url)
    
    
    if (
      status_code(response) != 200
    ) {
      break
    }
    
    
    data <- content(
      response,
      as = "text",
      encoding = "UTF-8"
    )
    
    
    data <- fromJSON(
      data,
      flatten = TRUE
    )
    
    
    if (
      is.null(data$entries)
    ) {
      break
    }
    
    
    all_records[
      [length(all_records) + 1]
    ] <- data$entries
    
    
    start_position <- start_position + batch_size
    
    
    Sys.sleep(0.3)
  }
  
  
  output <- do.call(
    rbind,
    all_records
  )
  
  
  message(
    "Downloaded records: ",
    nrow(output)
  )
  
  
  return(output)
}




# ==========================
# Ensembl ID conversion
# ==========================


convert_ensembl_to_symbol <- function(
    ensembl_ids,
    batch_size = 200
) {
  
  
  message(
    "\nConverting Ensembl IDs..."
  )
  
  
  unique_ids <- unique(
    ensembl_ids
  )
  
  
  mapping <- data.frame(
    Ensembl_ID = unique_ids,
    Gene_Symbol = NA_character_,
    stringsAsFactors = FALSE
  )
  
  
  total <- length(unique_ids)
  
  
  for (
    start in seq(
      1,
      total,
      by = batch_size
    )
  ) {
    
    
    end <- min(
      start + batch_size - 1,
      total
    )
    
    
    batch_ids <- unique_ids[start:end]
    
    
    response <- tryCatch(
      
      POST(
        "https://rest.ensembl.org/lookup/id",
        body = toJSON(
          list(ids = batch_ids)
        ),
        content_type(
          "application/json"
        ),
        accept(
          "application/json"
        )
      ),
      
      error = function(e)
        NULL
    )
    
    
    if (
      !is.null(response) &&
      status_code(response) == 200
    ) {
      
      result <- fromJSON(
        content(
          response,
          as = "text"
        )
      )
      
      
      for (
        id in names(result)
      ) {
        
        if (
          !is.null(
            result[[id]]$display_name
          )
        ) {
          
          mapping$Gene_Symbol[
            mapping$Ensembl_ID == id
          ] <- result[[id]]$display_name
          
        }
      }
    }
    
    
    Sys.sleep(0.5)
  }
  
  
  return(mapping)
}




# ==========================
# Data processing
# ==========================


process_gene_data <- function(
    df,
    disease_name
) {
  
  
  result <- data.frame(
    ID = df$id,
    Disease = disease_name,
    stringsAsFactors = FALSE
  )
  
  
  result$Ensembl_ID <- sapply(
    strsplit(
      result$ID,
      "-"
    ),
    function(x)
      x[1]
  )
  
  
  if (
    "source" %in% colnames(df)
  ) {
    
    result$Source <- df$source
    
  } else {
    
    result$Source <- NA
    
  }
  
  
  mapping <- convert_ensembl_to_symbol(
    result$Ensembl_ID
  )
  
  
  result$Gene_Symbol <- mapping$Gene_Symbol[
    match(
      result$Ensembl_ID,
      mapping$Ensembl_ID
    )
  ]
  
  
  result <- result[
    ,
    c(
      "ID",
      "Ensembl_ID",
      "Gene_Symbol",
      "Disease",
      "Source"
    )
  ]
  
  
  result <- result[
    order(
      result$Gene_Symbol
    ),
  ]
  
  
  return(result)
}




# ==========================
# Main workflow
# ==========================


main <- function() {
  
  
  raw_data <- download_ebi_genes(
    disease_name
  )
  
  
  if (
    is.null(raw_data)
  ) {
    
    stop(
      "No records retrieved."
    )
    
  }
  
  
  processed_data <- process_gene_data(
    raw_data,
    disease_name
  )
  
  
  write.csv(
    processed_data,
    output_file,
    row.names = FALSE
  )
  
  
  message(
    "\nSaved file: ",
    output_file
  )
  
  
  message(
    "Total genes: ",
    nrow(processed_data)
  )
  
  
  print(
    head(
      processed_data[
        ,
        c(
          "Ensembl_ID",
          "Gene_Symbol"
        )
      ],
      20
    )
  )
  
  
  return(processed_data)
}



# Execute pipeline

gene_results <- main()