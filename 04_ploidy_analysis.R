###############################################################################
# Script: 04_ploidy_analysis.R
# 
# Purpose:
#   Creates boxplots and frequency plots for ploidy data from seeds (embryo,
#   endosperm) and leaves. Analyzes ploidy variation across Crataegus species
#   and generates a multi-panel figure showing ploidy frequencies and indices.
#
# Input Files Required:
#   - Crataegus_seed_FCCS.xlsx: Seed ploidy data (embryo and endosperm)
#   - Crataegus_leaf_FCM.xlsx: Leaf ploidy data
#
# Output Files Generated:
#   - fcm_fccs_plots.pdf: Multi-panel figure with 8 plots showing ploidy frequencies
#     and indices for embryos, endosperms, and leaves
#
# Dependencies (with versions):
#   - readxl (1.4.5): Reading Excel files
#   - cowplot (1.2.0): Arranging multiple plots
#   - ggplot2 (3.5.2): Creating plots
#   - scales (1.4.0): Color scaling
#
# Author: Soňa Píšová
# Date: July 2025
# 
# Usage:
#   Place all input files in the same directory as this script, then run:
#   source("04_ploidy_analysis.R")
#   or
#   Rscript 04_ploidy_analysis.R
#
###############################################################################


# ============================================================================
# SECTION 0: Load Required Libraries
# ============================================================================

library(readxl)    # Reading Excel files
library(cowplot)   # Arranging multiple plots
library(ggplot2)   # Creating plots
library(scales)    # Color scaling


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

# Define input file names
seed_ploidy_file <- "Crataegus_seed_FCCS.xlsx"
leaf_ploidy_file <- "Crataegus_leaf_FCM.xlsx"


# ============================================================================
# SECTION 2: Load and Prepare Data
# ============================================================================

# --- 2.1 Load seed ploidy data ---
seed_data_raw <- read_excel(seed_ploidy_file)

# Convert to data frame with proper types
seed_data <- data.frame(
  Species = as.factor(seed_data_raw$Species),
  Standard = as.factor(seed_data_raw$Standard),
  Sample = as.character(seed_data_raw$Sample),
  Em_index = as.numeric(seed_data_raw$Em_index),
  En_index = as.numeric(seed_data_raw$En_index),
  Pl_embryo = as.numeric(seed_data_raw$Pl_embryo),
  Pl_endosperm = as.numeric(seed_data_raw$Pl_endosperm)
)

# Split by standard (Carex vs Pisum)
embryo_carex <- data.frame(seed_data[seed_data$Standard == "Carex", ])
embryo_pisum <- data.frame(seed_data[seed_data$Standard == "Pisum", ])

# --- 2.2 Load leaf ploidy data ---
leaf_data_raw <- read_excel(leaf_ploidy_file)

# Convert to data frame with proper types
leaf_data <- data.frame(
  Species = as.factor(leaf_data_raw$Species),
  Sample = as.character(leaf_data_raw$Sample),
  Leaf_index = as.numeric(leaf_data_raw$Leaf_index),
  Leaf_ploidy = as.numeric(leaf_data_raw$Leaf_ploidy)
)


# ============================================================================
# SECTION 3: Calculate Point Sizes Based on Frequency
# ============================================================================

# Species order for consistent plotting
species_order <- c("C. laevigata", "C. × media", "C. monogyna", 
                   "C. × subsphaerica", "C. rhipidophylla", 
                   "C. × macrocarpa", "C. sp.")

# --- 3.1 Calculate embryo ploidy frequencies and point sizes ---
embryo_freq_table <- table(seed_data$Species, as.character(seed_data$Pl_embryo))
embryo_freq_ordered <- embryo_freq_table[species_order, ]

# Calculate point sizes based on frequency (larger points for more common ploidy values)
embryo_freq_vector <- c()
for (i in 1:ncol(embryo_freq_ordered)) {
  embryo_freq_vector <- c(embryo_freq_vector, embryo_freq_ordered[, i])
}

# Assign size to each data point based on frequency
for (i in 1:nrow(seed_data)) {
  species <- as.character(seed_data$Species[i])
  ploidy <- as.character(seed_data$Pl_embryo[i])
  frequency <- embryo_freq_ordered[species, ploidy]
  seed_data[i, "size"] <- frequency
}

# Scale sizes: divide by 8 and add 2, with manual adjustments for specific values
seed_data$size <- (seed_data$size / 8) + 2
seed_data$size <- ifelse(seed_data$size == 57, 24, seed_data$size)
seed_data$size <- ifelse(seed_data$size == 24.250, 17, seed_data$size)
seed_data$size <- ifelse(seed_data$size == 16.25, 13, seed_data$size)

# --- 3.2 Calculate endosperm ploidy frequencies and point sizes ---
endosperm_freq_table <- table(seed_data$Species, as.character(seed_data$Pl_endosperm))
endosperm_freq_ordered <- endosperm_freq_table[species_order, ]

endosperm_freq_vector <- c()
for (i in 1:ncol(endosperm_freq_ordered)) {
  endosperm_freq_vector <- c(endosperm_freq_vector, endosperm_freq_ordered[, i])
}

# Assign size to each data point
for (i in 1:nrow(seed_data)) {
  species <- as.character(seed_data$Species[i])
  ploidy_en <- as.character(seed_data$Pl_endosperm[i])
  frequency_en <- endosperm_freq_ordered[species, ploidy_en]
  seed_data[i, "size_en"] <- frequency_en
}

# Scale sizes with manual adjustments
seed_data$size_en <- (seed_data$size_en / 8) + 2
seed_data$size_en <- ifelse(seed_data$size_en == 55.875, 16.7, seed_data$size_en)
seed_data$size_en <- ifelse(seed_data$size_en == 16, 11, seed_data$size_en)
seed_data$size_en <- ifelse(seed_data$size_en == 11.875, 3.5, seed_data$size_en)
seed_data$size_en <- ifelse(seed_data$size_en == 6.875, 3, seed_data$size_en)
seed_data$size_en <- ifelse(seed_data$size_en == 4.375, 2, seed_data$size_en)
seed_data$size_en <- ifelse(seed_data$size_en == 3.750, 2, seed_data$size_en)
seed_data$size_en <- ifelse(seed_data$size_en == 5.875, 4.5, seed_data$size_en)
seed_data$size_en <- ifelse(seed_data$size_en == 6.375, 4.1, seed_data$size_en)
seed_data$size_en <- ifelse(seed_data$size_en == 5.750, 3.2, seed_data$size_en)
seed_data$size_en <- ifelse(seed_data$size_en == 2.125, 2, seed_data$size_en)
seed_data$size_en <- ifelse(seed_data$size_en == 5.250, 2, seed_data$size_en)

# --- 3.3 Calculate leaf ploidy frequencies and point sizes ---
leaf_freq_table <- table(leaf_data$Species, as.character(leaf_data$Leaf_ploidy))
leaf_freq_ordered <- leaf_freq_table[species_order, ]

leaf_freq_vector <- c()
for (i in 1:ncol(leaf_freq_ordered)) {
  leaf_freq_vector <- c(leaf_freq_vector, leaf_freq_ordered[, i])
}

# Assign size to each data point
for (i in 1:nrow(leaf_data)) {
  species <- as.character(leaf_data$Species[i])
  ploidy_l <- as.character(leaf_data$Leaf_ploidy[i])
  frequency_l <- leaf_freq_ordered[species, ploidy_l]
  leaf_data[i, "size_l"] <- frequency_l
}

# Scale sizes with manual adjustments
leaf_data$size_l <- (leaf_data$size_l / 8) + 2
leaf_data$size_l <- ifelse(leaf_data$size_l == 57, 24, leaf_data$size_l)
leaf_data$size_l <- ifelse(leaf_data$size_l == 24.250, 17, leaf_data$size_l)
leaf_data$size_l <- ifelse(leaf_data$size_l == 16.25, 13, leaf_data$size_l)


# ============================================================================
# SECTION 4: Define Plotting Aesthetics
# ============================================================================

# Species colors and shapes (matching other scripts)
species_colors <- c("#7f7f7f", "black", "#7f7f7f", "black", "black", "black", "black")
species_fills <- c("#7f7f7f", "#7f7f7f", "gray70", "#eb5a82", "#28e21b", "#302790", "white")
species_shapes <- c(3, 21, 8, 22, 21, 23, 21)

# Common theme for all plots
common_theme <- theme(
  plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"),
  legend.position = "none",
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
  panel.background = element_blank(),
  text = element_text(color = "black"),
  axis.text.x = element_text(hjust = 0.95, vjust = 0.3, face = "bold.italic", 
                             color = "black", size = 10, angle = 90),
  axis.text.y = element_text(face = "bold", color = "black", size = 10, angle = 90),
  axis.title.y = element_text(vjust = +5),
  axis.title.x = element_blank()
)


# ============================================================================
# SECTION 5: Create Individual Plots
# ============================================================================

# --- 5.1 Embryo ploidy frequency plot ---
plot_embryo_ploidy <- ggplot(seed_data, aes(x = Species, y = Pl_embryo, group = Species)) +
  geom_point(aes(shape = Species, color = Species, fill = Species), 
             size = seed_data$size) +
  scale_shape_manual(values = species_shapes) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_fills) +
  scale_x_discrete(limits = species_order) +
  scale_y_continuous(name = "Ploidy of embryo", limits = c(1.8, 6)) +
  common_theme

# --- 5.2 Endosperm ploidy frequency plot ---
plot_endosperm_ploidy <- ggplot(seed_data, aes(x = Species, y = Pl_endosperm, group = Species)) +
  geom_point(aes(shape = Species, color = Species, fill = Species), 
             size = seed_data$size_en) +
  scale_shape_manual(values = species_shapes) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_fills) +
  scale_x_discrete(limits = species_order) +
  scale_y_continuous(name = "Ploidy of endosperm", limits = c(2.8, 16), 
                     breaks = seq(3, 16, by = 1)) +
  common_theme

# --- 5.3 Embryo index (Carex standard) plot ---
plot_embryo_carex_index <- ggplot(embryo_carex, 
                                   aes(x = Species, y = Em_index, group = Species)) +
  geom_jitter(width = 0.3, height = 0, 
              aes(shape = Species, color = Species, fill = Species), 
              size = 3) +
  scale_shape_manual(values = species_shapes) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_fills) +
  scale_x_discrete(limits = species_order) +
  scale_y_continuous(name = "Embryo index, Carex", limits = c(1.2, 4.2), 
                     breaks = seq(1.2, 4.2, by = 0.5)) +
  common_theme

# --- 5.4 Embryo index (Pisum standard) plot ---
plot_embryo_pisum_index <- ggplot(embryo_pisum, 
                                   aes(x = Species, y = Em_index, group = Species)) +
  geom_jitter(width = 0.3, height = 0, 
              aes(shape = Species, color = Species, fill = Species), 
              size = 3) +
  scale_shape_manual(values = species_shapes) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_fills) +
  scale_x_discrete(limits = species_order) +
  scale_y_continuous(name = "Embryo index, Pisum", limits = c(0.15, 0.6), 
                     breaks = seq(0.1, 0.6, by = 0.1)) +
  common_theme

# --- 5.5 Endosperm index (Carex standard) plot ---
plot_endosperm_carex_index <- ggplot(embryo_carex, 
                                      aes(x = Species, y = En_index, group = Species)) +
  geom_jitter(width = 0.3, height = 0, 
              aes(shape = Species, color = Species, fill = Species), 
              size = 3) +
  scale_shape_manual(values = species_shapes) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_fills) +
  scale_x_discrete(limits = species_order) +
  scale_y_continuous(name = "Endosperm index, Carex", limits = c(1.8, 11), 
                     breaks = seq(1.8, 11, by = 1)) +
  common_theme

# --- 5.6 Endosperm index (Pisum standard) plot ---
plot_endosperm_pisum_index <- ggplot(embryo_pisum, 
                                      aes(x = Species, y = En_index, group = Species)) +
  geom_jitter(width = 0.3, height = 0, 
              aes(shape = Species, color = Species, fill = Species), 
              size = 3) +
  scale_shape_manual(values = species_shapes) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_fills) +
  scale_x_discrete(limits = species_order) +
  scale_y_continuous(name = "Endosperm index, Pisum", limits = c(0.2, 1.6), 
                     breaks = seq(0.2, 1.6, by = 0.2)) +
  common_theme

# --- 5.7 Leaf ploidy frequency plot ---
plot_leaf_ploidy <- ggplot(leaf_data, 
                            aes(x = Species, y = Leaf_ploidy, group = Species)) +
  geom_point(aes(shape = Species, color = Species, fill = Species), 
             size = leaf_data$size_l) +
  scale_shape_manual(values = species_shapes) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_fills) +
  scale_x_discrete(limits = species_order) +
  scale_y_continuous(name = "Ploidy of trees", limits = c(1.95, 4), 
                     breaks = seq(2, 4, by = 1)) +
  common_theme

# --- 5.8 Leaf ploidy index plot ---
plot_leaf_index <- ggplot(leaf_data, 
                          aes(x = Species, y = Leaf_index, group = Species)) +
  geom_jitter(width = 0.3, height = 0, 
              aes(shape = Species, color = Species, fill = Species), 
              size = 3) +
  scale_shape_manual(values = species_shapes) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_fills) +
  scale_x_discrete(limits = species_order) +
  scale_y_continuous(name = "Leaf index", limits = c(0.18, 0.4), 
                     breaks = seq(0.18, 0.4, by = 0.03)) +
  common_theme


# ============================================================================
# SECTION 6: Combine and Save Plots
# ============================================================================

# Combine all plots into a single figure
combined_plot <- plot_grid(
  plot_embryo_ploidy,
  plot_endosperm_ploidy,
  plot_embryo_carex_index,
  plot_embryo_pisum_index,
  plot_endosperm_carex_index,
  plot_endosperm_pisum_index,
  plot_leaf_ploidy,
  plot_leaf_index,
  labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
  ncol = 2,
  nrow = 4
)

# Save to PDF
pdf(file = "fcm_fccs_plots.pdf", width = 12, height = 25, paper = 'special')
print(combined_plot)
dev.off()

