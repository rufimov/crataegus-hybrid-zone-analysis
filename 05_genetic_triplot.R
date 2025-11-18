###############################################################################
# Script: 05_genetic_triplot.R
# 
# Purpose:
#   Creates a ternary plot (triplot) showing admixture proportions of three
#   Crataegus species (C. monogyna, C. laevigata, C. rhipidophylla) based on
#   SSR genetic data. Each sample is plotted as a point in the ternary space
#   with colors and shapes indicating genetic groups.
#
# Input Files Required:
#   - Crataegus_ssr_triplot.csv: CSV file with admixture proportions
#     Columns: Sample name, K6 (genetic group), M (C. monogyna proportion),
#              L (C. laevigata proportion), R (C. rhipidophylla proportion)
#
# Output Files Generated:
#   - triplot.pdf: Ternary plot showing admixture proportions
#
# Dependencies (with versions):
#   - Base R graphics (no additional packages required)
#
# Author: Soňa Píšová
# Date: July 2025
# 
# Usage:
#   Place all input files in the same directory as this script, then run:
#   source("05_genetic_triplot.R")
#   or
#   Rscript 05_genetic_triplot.R
#
#
###############################################################################


# ============================================================================
# SECTION 1: Setup and File Paths
# ============================================================================

#' Detect script directory (works in RStudio and Rscript)
#' This makes the script portable - it will work regardless of where it's run from
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(sub("^--file=", "", file_arg)))
  } else if (exists("rstudioapi") && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  } else {
    return(getwd())
  }
}

# Set working directory to script location
script_dir <- get_script_dir()
setwd(script_dir)

# Define input file name
triplot_data_file <- "Crataegus_ssr_triplot.csv"


# ============================================================================
# SECTION 2: Load Data
# ============================================================================

# Load triplot data (admixture proportions)
triplot_data <- read.table(triplot_data_file, sep = ",", header = TRUE)

# Extract genetic group information
genetic_groups <- as.factor(triplot_data$K6)


# ============================================================================
# SECTION 3: Create Ternary Grid
# ============================================================================

#' Transform ternary coordinates to Cartesian coordinates
#' 
#' Converts coordinates from ternary space (three proportions that sum to 1)
#' to Cartesian (x, y) coordinates for plotting in an equilateral triangle.
#' 
#' @param coord Vector of length 3: (x, y, z) proportions
#' @return Vector of length 2: (x1, y1) Cartesian coordinates
tern2cart <- function(coord) {
  x <- coord[1]
  y <- coord[2]
  z <- coord[3]
  tot <- x + y + z
  x <- x / tot
  y <- y / tot
  z <- z / tot
  x1 <- (2 * y + z) / (2 * (x + y + z))
  y1 <- sqrt(3) * z / (2 * (x + y + z))
  return(c(x1, y1))
}

# Create grid lines for ternary plot (10% intervals)
# Grid coordinates: three sets of lines (parallel to each side)
a <- seq(0.9, 0.1, by = -0.1)  # Decreasing from 0.9 to 0.1
b <- rep(0, 9)                  # Constant at 0
c <- seq(0.1, 0.9, by = 0.1)   # Increasing from 0.1 to 0.9

# Combine into grid data frame
ternary_grid <- data.frame(
  x = c(a, b, c, a, c, b),
  y = c(b, c, a, c, b, a),
  z = c(c, a, b, b, a, c)
)

# Transform grid coordinates to Cartesian
grid_cartesian <- t(apply(ternary_grid, 1, tern2cart))

# Prepare grid segments (pairs of points for line drawing)
grid_segments <- cbind(grid_cartesian[1:27, ], grid_cartesian[28:54, ])


# ============================================================================
# SECTION 4: Prepare Sample Data for Plotting
# ============================================================================

# Create data frame with species proportions
sample_triplot <- data.frame(
  'Cr.monogyna' = c(triplot_data[, 3]),      # M column
  'Cr.laevigata' = c(triplot_data[, 4]),     # L column
  'Cr.rhipidophylla' = c(triplot_data[, 5]), # R column
  row.names = as.character(triplot_data[, 1])
)

# Transform sample coordinates to Cartesian
sample_cartesian <- t(apply(sample_triplot, 1, tern2cart))


# ============================================================================
# SECTION 5: Define Plotting Parameters
# ============================================================================

# Point characteristics based on genetic group (K6)
# Groups: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
point_shapes <- c(22, 22, 20, 21, 8, 23, 23, 23, 23, 3)[as.numeric(triplot_data$K6)]
point_colors <- c("black", "black", "gray50", "black", "gray50", "black", 
                  "black", "black", "black", "gray50")[as.numeric(triplot_data$K6)]
point_fills <- c("orange", "red", "gray50", "green", "gray50", "blue", 
                 "yellow", "white", "brown", "black")[as.numeric(triplot_data$K6)]
point_sizes <- c(2, 2, 2, 2, 1.7, 2, 2, 2, 2, 1.7)[as.numeric(triplot_data$K6)]

# Legend labels for genetic groups
legend_labels <- c(
  expression(paste(italic("C"), ". ", italic("laevigata"), " group 1")),
  expression(paste(italic("C"), ". ", italic("laevigata"), " group 2")),
  expression(paste(italic("C"), ". ", italic("laevigata"), " 1 x 2")),
  expression(paste(italic("C"), ". x ", italic("media"))),
  expression(paste(italic("C"), ". ", italic("monogyna"))),
  expression(paste(italic("C"), ". x ", italic("subsphaerica"))),
  expression(paste(italic("C"), ". ", italic("rhipidophylla"), " group 1")),
  expression(paste(italic("C"), ". ", italic("rhipidophylla"), " group 2")),
  expression(paste(italic("C"), ". ", italic("rhipidophylla"), " 1 x 2")),
  expression(paste(italic("C"), ". ", italic("rhipidophylla"), " group 3")),
  expression(paste(italic("C"), ". x ", italic("macrocarpa")))
)

legend_colors <- c("black", "black", "black", "gray50", "black", "gray50", 
                   "black", "black", "black", "black", "gray50")
legend_fills <- c("red", "orange", "white", "gray50", "green", "gray50", 
                  "blue", "brown", "white", "yellow", "black")
legend_shapes <- c(22, 22, 22, 20, 21, 8, 23, 23, 23, 23, 3)
legend_sizes <- c(2, 2, 2, 2, 1.7, 2, 2, 2, 2, 1.7)


# ============================================================================
# SECTION 6: Create Ternary Plot
# ============================================================================

# Open PDF device for saving plot
pdf("triplot.pdf", width = 8, height = 8)

# Set up graphics device
par(mar = c(3, 3, 3, 3) + 0.1, xpd = TRUE)

# Create empty plot with correct dimensions for equilateral triangle
plot(NA, NA, 
     xlim = c(0, 1), 
     ylim = c(0, sqrt(3) / 2), 
     bty = "n", 
     asp = 1, 
     axes = FALSE, 
     xlab = "", 
     ylab = "")

# Draw triangle outline
segments(0, 0, 0.5, sqrt(3) / 2, lwd = 3)  # Left side
segments(0.5, sqrt(3) / 2, 1, 0, lwd = 3)  # Right side
segments(1, 0, 0, 0, lwd = 3)                # Bottom side

# Add percentage labels (10% to 90%)
percentage_labels <- paste(seq(10, 90, by = 10), "%")
text(grid_cartesian[9:1, ], paste(percentage_labels), 
     col = "black", cex = 1, pos = 2)      # Left side labels
text(grid_cartesian[18:10, ], paste(percentage_labels), 
     col = "black", cex = 1, pos = 4)      # Right side labels
text(grid_cartesian[27:19, ], paste(percentage_labels), 
     col = "black", cex = 1, pos = 1)      # Bottom labels

# Add species labels at triangle vertices
text(-0.03, -0.001, 
     labels = expression(paste(bolditalic("C. monogyna"))), 
     pos = 1, cex = 1)
text(1.03, -0.001, 
     labels = expression(paste(bolditalic("C. laevigata"))), 
     pos = 1, cex = 1)
text(0.5, 0.87, 
     labels = expression(paste(bolditalic("C. rhipidophylla"))), 
     cex = 1, pos = 3)

# Draw grid lines
apply(grid_segments, 1, function(x) {
  segments(x0 = x[1], y0 = x[2], x1 = x[3], y1 = x[4], 
           lty = 1, col = "gray")
})

# Plot sample points
points(sample_cartesian, 
       pch = point_shapes, 
       col = point_colors, 
       bg = point_fills, 
       cex = point_sizes)

# Add legend
legend(0.8, 0.91, 
       bty = "o", 
       legend = legend_labels,
       col = legend_colors,
       pt.bg = legend_fills, 
       pch = legend_shapes, 
       pt.cex = legend_sizes, 
       ncol = 1, 
       cex = 1)

# Close PDF device
dev.off()
