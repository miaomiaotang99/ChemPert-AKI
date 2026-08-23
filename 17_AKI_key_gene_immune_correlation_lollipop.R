############################################################
# AKI key gene and immune cell correlation lollipop plot
#
# Description:
#   Visualize correlations between key AKI-associated genes
#   and immune cell fractions estimated by CIBERSORT.
#
# Input:
#   Gene correlation result files:
#   - TP53_cor.result.txt
#   - STAT3_cor.result.txt
#   - AKT1_cor.result.txt
#
# Required columns:
#   Cell
#   cor
#   pvalue
#
# Output:
#   Lollipop correlation plot
#
# Author:
#   Weidong Huang
#
############################################################



# ============================================================
# 1. Load packages and configuration
# ============================================================


rm(list = ls())



working_directory <- "C:/user"


setwd(
  working_directory
)



# Select gene correlation result

gene_name <- "AKT1"


input_file <- paste0(
  gene_name,
  "_cor.result.txt"
)



output_file <- paste0(
  gene_name,
  "_immune_correlation_lollipop.pdf"
)




# ============================================================
# 2. Read correlation results
# ============================================================


cor_data <- read.table(

  input_file,

  header = TRUE,

  sep = "\t",

  check.names = FALSE

)



head(cor_data)




# ============================================================
# 3. Define visualization parameters
# ============================================================


# Color scale based on p-value

p_colors <- c(
  "gold",
  "pink",
  "orange",
  "LimeGreen",
  "darkgreen"
)



get_point_color <- function(
    pvalue,
    colors
){

  color <- ifelse(
    
    pvalue > 0.8,
    colors[1],
    
    ifelse(
      
      pvalue > 0.6,
      colors[2],
      
      ifelse(
        
        pvalue > 0.4,
        colors[3],
        
        ifelse(
          
          pvalue > 0.2,
          colors[4],
          
          colors[5]
        )
      )
    )
  )


  return(color)

}





# Point size based on correlation strength


point_sizes <- seq(
  2.5,
  5.5,
  length = 5
)



get_point_size <- function(
    correlation
){

  correlation <- abs(
    correlation
  )


  size <- ifelse(

    correlation < 0.1,

    point_sizes[1],

    ifelse(

      correlation < 0.2,

      point_sizes[2],

      ifelse(

        correlation < 0.3,

        point_sizes[3],

        ifelse(

          correlation < 0.4,

          point_sizes[4],

          point_sizes[5]

        )
      )
    )
  )


  return(size)

}




# Add visualization parameters

cor_data$point_color <- get_point_color(
  cor_data$pvalue,
  p_colors
)



cor_data$point_size <- get_point_size(
  cor_data$cor
)



# Sort by correlation

cor_data <- cor_data[
  order(cor_data$cor),
]




# ============================================================
# 4. Draw lollipop plot
# ============================================================


x_limit <- ceiling(
  max(
    abs(cor_data$cor)
  )*10
)/10




pdf(

  output_file,

  width = 9,

  height = 7

)



layout(

  mat = matrix(
    c(
      1,1,1,1,1,
      0,2,0,3,0
    ),
    nc = 2
  ),

  width = c(
    8,
    2.2
  )

)



par(

  bg = "white",

  las = 1,

  mar = c(
    5,
    18,
    2,
    4
  ),

  cex.axis = 1.5,

  cex.lab = 2

)



plot(

  1,

  type = "n",

  xlim = c(
    -x_limit,
    x_limit
  ),

  ylim = c(
    0.5,
    nrow(cor_data)+0.5
  ),

  xlab = "Correlation coefficient",

  ylab = "",

  yaxt = "n",

  axes = FALSE

)



rect(

  par("usr")[1],

  par("usr")[3],

  par("usr")[2],

  par("usr")[4],

  col = "#F5F5F5",

  border = "#F5F5F5"

)



grid(

  ny = nrow(cor_data),

  col = "white",

  lty = 1,

  lwd = 2

)



# Correlation line

segments(

  x0 = cor_data$cor,

  y0 = 1:nrow(cor_data),

  x1 = 0,

  y1 = 1:nrow(cor_data),

  lwd = 4

)



# Correlation point

points(

  x = cor_data$cor,

  y = 1:nrow(cor_data),

  col = cor_data$point_color,

  pch = 16,

  cex = cor_data$point_size

)




# Cell names

text(

  par("usr")[1],

  1:nrow(cor_data),

  cor_data$Cell,

  adj = 1,

  xpd = TRUE,

  cex = 1.5

)




# P values

p_label <- ifelse(

  cor_data$pvalue < 0.001,

  "<0.001",

  sprintf(
    "%.03f",
    cor_data$pvalue
  )

)



text(

  par("usr")[2],

  1:nrow(cor_data),

  p_label,

  adj = 0,

  xpd = TRUE,

  cex = 1.5,

  col = ifelse(
    
    abs(cor_data$cor)>0 &
      cor_data$pvalue<0.05,
    
    "red",
    
    "black"
    
  )

)



axis(

  1,

  tick = FALSE

)



# ============================================================
# 5. Point size legend
# ============================================================


par(
  mar=c(
    0,
    4,
    3,
    4
  )
)



plot(

  1,

  type="n",

  axes=FALSE

)



legend(

  "left",

  legend=c(
    0.1,
    0.2,
    0.3,
    0.4,
    0.5
  ),

  col="black",

  pt.cex=point_sizes,

  pch=16,

  bty="n",

  cex=2,

  title="|Correlation|"

)



# ============================================================
# 6. P-value color legend
# ============================================================


par(
  mar=c(
    0,
    6,
    4,
    6
  )
)



barplot(

  rep(
    1,
    5
  ),

  horiz = TRUE,

  space = 0,

  border = NA,

  col = p_colors,

  xaxt = "n",

  yaxt = "n",

  main = "p-value"

)



axis(

  4,

  at = 0:5,

  labels = c(
    1,
    0.8,
    0.6,
    0.4,
    0.2,
    0
  ),

  tick = FALSE

)



dev.off()



message(
  "Lollipop plot completed: ",
  output_file
)