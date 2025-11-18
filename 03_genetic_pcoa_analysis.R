###############################################################################
# Script: 03_genetic_pcoa_analysis.R
# 
# Purpose:
#   Performs Principal Coordinates Analysis (PCoA) on SSR (microsatellite)
#   genetic data using Bruvo distance. Creates PCoA plots with density
#   contours and hybrid coloring based on admixture proportions.
#
# Input Files Required:
#   - Crataegus_ssr_data_sp_clones_removed_dosage_restored.gdv: GenoDive format
#     file with SSR genetic data (ploidy and repeat unit information required)
#   - reassign.tsv: Admixture and species assignment table
#
# Output Files Generated:
#   - genetic_PCoA_densities_ggplot.pdf: PCoA plot with density contours
#
# Dependencies (with versions):
#   - polysat (1.7.7): For reading GenoDive files and calculating Bruvo distance
#   - tidyverse (2.0.0): Data manipulation and plotting
#   - scales (1.4.0): Color scaling
#   - ggnewscale (0.5.2): Multiple fill scales in ggplot
#   - ape (5.8.1): For PCoA (pcoa function)
#
# Author: Roman Ufimov
# Date: July 2025
# 
# Usage:
#   Place all input files in the same directory as this script, then run:
#   source("03_genetic_pcoa_analysis.R")
#   or
#   Rscript 03_genetic_pcoa_analysis.R
#
# Note:
#   The ploidy and repeat unit (Usatnts) vectors are hardcoded and must match
#   the order of samples in the GenoDive file. Adjust these if your data
#   structure differs.
#
###############################################################################


# ============================================================================
# SECTION 0: Load Required Libraries
# ============================================================================

library(polysat)      # SSR data handling and Bruvo distance
library(tidyverse)    # Data manipulation and plotting
library(scales)       # Color scaling
library(ggnewscale)   # Multiple fill scales in ggplot


# ============================================================================
# SECTION 1: Setup and File Paths
# ============================================================================

#' Detect script directory (works in RStudio and Rscript)
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

# Define input file names
genodive_file <- "Crataegus_ssr_data_sp_clones_removed_dosage_restored.gdv"
admixture_file <- "reassign.tsv"


# ============================================================================
# SECTION 2: Load and Prepare Genetic Data
# ============================================================================

# --- 2.1 Read GenoDive file ---
genetic_data <- read.GenoDive(genodive_file)
sample_names_genetic <- Samples(genetic_data)

# --- 2.2 Set ploidy for each sample ---
# IMPORTANT: This vector must match the order of samples in the GenoDive file
# Values: 2 = diploid, 3 = triploid, 4 = tetraploid
Ploidies(genetic_data) <- c(3,3,2,2,3,2,2,2,2,2,3,2,4,3,3,2,3,3,3,3,3,4,3,3,3,3,3,3,4,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,2,3,2,2,3,3,3,3,3,3,4,3,4,3,4,4,4,3,3,3,3,3,3,4,3,3,4,3,4,3,3,3,4,4,3)

# --- 2.3 Set repeat unit lengths for each locus ---
# IMPORTANT: This vector must match the order of loci in the GenoDive file
# Values: 2 = dinucleotide repeat, 3 = trinucleotide repeat
Usatnts(genetic_data) <- c(2,2,2,3,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2)

# Display summary of data
summary(genetic_data)

# --- 2.4 Set population names and assignments ---
# Population names (11 populations)
PopNames(genetic_data) <- c("1","2","3","4","5","6","7","8","9","10","11")

# Population assignment for each sample
# IMPORTANT: This vector must match the order of samples in the GenoDive file
PopInfo(genetic_data) <- c(6,10,3,3,6,2,2,2,2,2,9,3,10,10,10,1,10,10,10,10,10,10,10,10,10,1,10,10,10,10,3,3,10,3,3,3,2,3,3,3,2,4,3,3,9,2,2,2,2,2,2,2,2,11,2,2,10,2,2,2,3,2,2,2,2,2,2,2,1,3,2,2,11,11,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,11,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,4,4,4,4,4,4,3,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,4,4,5,5,6,4,4,5,4,5,4,6,6,9,9,10,10,9,7,10,7,7,7,10,10,10,6,1,6,6,5,5,5,5,10,5,10,8,7,5,6)


# ============================================================================
# SECTION 3: Load Admixture and Species Data
# ============================================================================

# Parental species colors
parental_colors <- c(L = "#FF5980", M = "#28e21b", R = "#302790")

#' Mix two colors based on admixture weights
#' @param col1 First color (hex code)
#' @param col2 Second color (hex code)
#' @param w1 Weight for first color (admixture proportion)
#' @param w2 Weight for second color (admixture proportion)
#' @return Mixed color as hex code
mix_colors <- function(col1, col2, w1, w2) {
  w1 <- as.numeric(w1)
  w2 <- as.numeric(w2)
  if (anyNA(c(w1, w2))) return(NA_character_)
  rgb( (w1*col2rgb(col1)[1,] + w2*col2rgb(col2)[1,]) / 255,
       (w1*col2rgb(col1)[2,] + w2*col2rgb(col2)[2,]) / 255,
       (w1*col2rgb(col1)[3,] + w2*col2rgb(col2)[3,]) / 255 )
}

# Load admixture table
reassign <- read.delim(admixture_file, 
                      sep = "\t", 
                      stringsAsFactors = FALSE, 
                      na.strings = c("#N/A")) %>%
  tidyr::separate(SSRs.group_new, into = "group", sep = "_", remove = FALSE) %>%
  mutate(species = gsub("\\d", "", group))  # Extract species code


# ============================================================================
# SECTION 4: Prepare Plotting Information
# ============================================================================

# Match colors and shapes to sample order in genetic data
plot_info <- tibble(sample = sample_names_genetic) %>%
  left_join(reassign, by = c("sample" = "Name")) %>%
  mutate(across(c(L, M, R), as.numeric),
         # Create mixed colors for hybrids based on admixture proportions
         mixed_col = case_when(
           species == "L" ~ parental_colors["L"],
           species == "M" ~ parental_colors["M"],
           species == "R" ~ parental_colors["R"],
           # MAC hybrids (L x R)
           species == "MAC" ~ mix_colors(parental_colors["L"], 
                                        parental_colors["R"], 
                                        L, R),
           # MED hybrids (L x M)
           species == "MED" ~ mix_colors(parental_colors["L"], 
                                        parental_colors["M"], 
                                        L, M),
           # SUB hybrids (R x M)
           species == "SUB" ~ mix_colors(parental_colors["R"], 
                                        parental_colors["M"], 
                                        R, M),
           TRUE ~ NA_character_
         ))


# ============================================================================
# SECTION 5: Calculate Genetic Distance and PCoA
# ============================================================================

# --- 5.1 Calculate Bruvo distance matrix ---
# Bruvo distance is appropriate for microsatellite data with variable ploidy
distance_matrix <- meandistance.matrix(genetic_data)

# --- 5.2 Perform Principal Coordinates Analysis ---
# PCoA is equivalent to classical MDS on the distance matrix
pcoa_results <- ape::pcoa(distance_matrix)

# Extract PCoA coordinates (eigenvectors)
pcoa_coordinates <- pcoa_results$vectors

# --- 5.3 Add PCoA coordinates to plot information ---
plot_info <- plot_info %>%
  mutate(
    PC1 = pcoa_coordinates[sample, 1],
    PC2 = pcoa_coordinates[sample, 2]
  )


# ============================================================================
# SECTION 6: Create PCoA Plot with Density Contours
# ============================================================================

# --- 6.1 Create density contour scales ---
n_bins <- 100
density_scale_L <- scales::alpha("#FF5980", seq(0, 0.4, length.out = n_bins))
density_scale_M <- scales::alpha("#28e21b", seq(0, 0.4, length.out = n_bins))
density_scale_R <- scales::alpha("#302790", seq(0, 0.4, length.out = n_bins))

# --- 6.2 Create plot ---
pcoa_plot <- ggplot(plot_info, aes(PC1, PC2)) +
  # Density contours for L species
  geom_density_2d_filled(
    data = plot_info %>% filter(species == "L"),
    aes(fill = after_stat(level)),
    contour_var = "ndensity",
    bins = n_bins,
    show.legend = FALSE,
    color = NA, 
    linetype = 0
  ) +
  scale_fill_manual(values = density_scale_L) +
  new_scale_fill() +
  # Density contours for M species
  geom_density_2d_filled(
    data = plot_info %>% filter(species == "M"),
    aes(fill = after_stat(level)),
    contour_var = "ndensity",
    bins = n_bins,
    show.legend = FALSE,
    color = NA, 
    linetype = 0
  ) +
  scale_fill_manual(values = density_scale_M) +
  new_scale_fill() +
  # Density contours for R species
  geom_density_2d_filled(
    data = plot_info %>% filter(species == "R"),
    aes(fill = after_stat(level)),
    contour_var = "ndensity",
    bins = n_bins,
    show.legend = FALSE,
    color = NA, 
    linetype = 0
  ) +
  scale_fill_manual(values = density_scale_R) +
  # Sample points
  geom_point(aes(colour = mixed_col, shape = species),
             size = 2) +
  scale_colour_identity() +
  scale_shape_manual(
    values = c(L = 15, M = 16, R = 18, MAC = 3, MED = 4, SUB = 8),
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  # Axis limits and labels
  xlim(-0.47, 0.32) +
  ylim(-0.27, 0.35) +
  labs(
    title = "Genetic PCoA (Bruvo distance) with densities",
    x = sprintf("PCoA axis 1 (%.1f%%)", 
                pcoa_results$values$Relative_eig[1] * 100),
    y = sprintf("PCoA axis 2 (%.1f%%)", 
                pcoa_results$values$Relative_eig[2] * 100),
    shape = "Species"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

# --- 6.3 Save plot ---
ggsave("genetic_PCoA_densities_ggplot.pdf", 
       plot = pcoa_plot, 
       width = 8.4, 
       height = 6)

# Display plot
print(pcoa_plot)

