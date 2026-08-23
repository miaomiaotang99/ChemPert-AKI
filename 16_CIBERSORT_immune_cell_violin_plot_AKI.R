############################################################
# CIBERSORT immune cell fraction comparison between Control
# and AKI groups
#
# Description:
#   Visualize immune cell infiltration differences estimated
#   by CIBERSORT and identify significantly altered immune
#   cell populations between Control and AKI samples.
#
# Input:
#   CIBERSORT_results.txt
#
# Output:
#   - Violin plot of immune cell fractions
#   - Differential immune cell analysis table
#   - Significant immune cell list
#
# Author:
#   Weidong Huang
#
############################################################



# ============================================================
# 1. Load packages and set working directory
# ============================================================


rm(list = ls())


library(vioplot)



working_directory <- "C:/user"


setwd(
  working_directory
)



input_file <- "CIBERSORT_results.txt"




# ============================================================
# 2. Read CIBERSORT results
# ============================================================


cibersort_data <- read.table(

  input_file,

  header = TRUE,

  sep = "\t",

  check.names = FALSE,

  row.names = 1

)



dim(cibersort_data)


head(cibersort_data)




# ============================================================
# 3. Define sample groups
# ============================================================


sample_group <- ifelse(

  grepl(
    "_Control$",
    rownames(cibersort_data)
  ),

  "Control",

  ifelse(

    grepl(
      "_Treat$",
      rownames(cibersort_data)
    ),

    "AKI",

    "Unknown"

  )

)



print(
  table(sample_group)
)



if(any(sample_group=="Unknown")){

  warning(
    "Unknown samples detected:"
  )

  print(
    rownames(cibersort_data)[
      sample_group=="Unknown"
    ]
  )

}




# ============================================================
# 4. Evaluate CIBERSORT fitting quality
# ============================================================


summary(
  cibersort_data$`P-value`
)


table(
  cibersort_data$`P-value` < 0.05
)



# ============================================================
# 5. Extract immune cell fractions
# ============================================================


remove_columns <- c(

  "P-value",

  "Correlation",

  "RMSE"

)



immune_fraction <- cibersort_data[

  ,

  setdiff(
    colnames(cibersort_data),
    remove_columns
  ),

  drop = FALSE

]



immune_fraction <- as.data.frame(

  lapply(
    immune_fraction,
    as.numeric
  )

)



rownames(immune_fraction) <- rownames(
  cibersort_data
)



dim(
  immune_fraction
)



colnames(
  immune_fraction
)




# ============================================================
# 6. Separate Control and AKI samples
# ============================================================


control_fraction <- immune_fraction[
  sample_group=="Control",
  ,
  drop=FALSE
]


aki_fraction <- immune_fraction[
  sample_group=="AKI",
  ,
  drop=FALSE
]



cat(
  "Control samples:",
  nrow(control_fraction),
  "\n"
)


cat(
  "AKI samples:",
  nrow(aki_fraction),
  "\n"
)


cat(
  "Immune cell types:",
  ncol(immune_fraction),
  "\n"
)




# ============================================================
# 7. Violin plot and statistical analysis
# ============================================================


output_statistics <- data.frame()



pdf(

  "CIBERSORT_AKI_immune_violin_plot.pdf",

  width = 12,

  height = 6

)



par(

  las = 1,

  mar = c(
    11,
    6,
    3,
    3
  )

)



max_value <- max(
  immune_fraction,
  na.rm = TRUE
)



plot(

  1,

  type = "n",

  xlim = c(
    -1,
    3*ncol(immune_fraction)
  ),

  ylim = c(
    0,
    max_value+0.08
  ),

  xlab = "",

  ylab = "Immune cell fraction",

  xaxt = "n"

)



for(i in 1:ncol(immune_fraction)){


  control_values <- as.numeric(
    control_fraction[,i]
  )


  aki_values <- as.numeric(
    aki_fraction[,i]
  )


  if(sd(control_values)==0){

    control_values[1] <-
      control_values[1]+0.00001

  }


  if(sd(aki_values)==0){

    aki_values[1] <-
      aki_values[1]+0.00001

  }



  vioplot(

    control_values,

    at = 3*(i-1),

    col = "#F1BB72",

    add = TRUE

  )



  vioplot(

    aki_values,

    at = 3*(i-1)+1,

    col = "#6778AE",

    add = TRUE

  )



  test_result <- wilcox.test(

    control_values,

    aki_values,

    exact = FALSE

  )



  p_value <- test_result$p.value



  output_statistics <- rbind(

    output_statistics,

    data.frame(

      CellType =
        colnames(immune_fraction)[i],

      Control_mean =
        mean(control_values),

      AKI_mean =
        mean(aki_values),

      p_value =
        p_value

    )

  )



  ymax <- max(
    c(
      control_values,
      aki_values
    ),
    na.rm = TRUE
  )


  lines(

    c(
      3*(i-1)+0.2,
      3*(i-1)+0.8
    ),

    c(
      ymax,
      ymax
    )

  )


  label <- ifelse(

    p_value < 0.001,

    "p<0.001",

    paste0(
      "p=",
      sprintf(
        "%.03f",
        p_value
      )
    )

  )


  text(

    x =
      3*(i-1)+0.5,

    y =
      ymax+0.025,

    labels = label,

    cex = 0.7

  )

}



legend(

  "topright",

  legend = c(
    "Control",
    "AKI"
  ),

  col = c(
    "#F1BB72",
    "#6778AE"
  ),

  lwd = 3,

  bty = "n"

)



text(

  seq(
    0.5,
    3*(ncol(immune_fraction)-1)+0.5,
    3
  ),

  -0.025,

  labels =
    colnames(immune_fraction),

  xpd = NA,

  cex = 0.8,

  srt = 45,

  pos = 2

)



dev.off()




# ============================================================
# 8. Multiple testing correction
# ============================================================


output_statistics$FDR <- p.adjust(

  output_statistics$p_value,

  method = "BH"

)



output_statistics$Significance <- ifelse(

  output_statistics$p_value < 0.001,

  "***",

  ifelse(

    output_statistics$p_value < 0.01,

    "**",

    ifelse(

      output_statistics$p_value < 0.05,

      "*",

      "ns"

    )

  )

)



output_statistics <- output_statistics[

  order(
    output_statistics$p_value
  ),

]



# Save complete statistics

write.table(

  output_statistics,

  file =
    "CIBERSORT_AKI_immune_cell_difference.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)




# ============================================================
# 9. Save significant immune cells
# ============================================================


significant_cells <- output_statistics[

  output_statistics$p_value < 0.05,

  ,

  drop = FALSE

]



write.table(

  significant_cells,

  file =
    "CIBERSORT_significant_immune_cells.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)



message(
  "\nCIBERSORT immune cell comparison completed successfully."
)