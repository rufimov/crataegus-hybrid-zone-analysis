###############################################################################
# Script: 01_relative_warp_analysis.R
# 
# Purpose:
#   Performs Relative Warp (RW) analysis on geometric morphometric data from
#   TPS landmark files. Analyzes flowering shoots and short shoots separately,
#   then combines them for a joint analysis. Also includes fruit PCA analysis.
#
# Input Files Required:
#   - Hybrid_Zones_LM_fl.tps: TPS file with landmark coordinates for flowering shoots
#   - Hybrid_Zones_LM_sh.tps: TPS file with landmark coordinates for short shoots
#   - reassign.tsv: Admixture and species assignment table (columns: Name, species, L, M, R)
#   - fruits_avg.tsv: Averaged fruit measurements (for section 7)
#
# Output Files Generated:
#   - relative_warp_morphospace_MorphoRW_combined_paired.pdf: Combined morphospace plot
#   - TPS_grid_combined_flowering_paired_tiled.pdf: Tiled TPS grids for combined flowering portion
#   - TPS_grid_combined_short_paired_tiled.pdf: Tiled TPS grids for combined short portion
#   - fruits_PCA_cent.pdf: Fruit PCA plot with hybrid coloring
#
# Dependencies (with versions):
#   - tidyverse (2.0.0), geomorph (4.0.10), ggrepel (0.9.6), prWarp (1.0.1)
#   - abind (1.4.8), here (1.0.1), ggnewscale (0.5.2)
#   - Morpho (2.13) - for relWarps function
#   - FactoMineR (2.12), factoextra (1.0.7) - for fruit PCA
#
# Author: Roman Ufimov
# Date: July 2025
# 
# Usage:
#   Place all input files in the same directory as this script, then run:
#   source("01_relative_warp_analysis.R")
#   or
#   Rscript 01_relative_warp_analysis.R
#
###############################################################################


# ============================================================================
# SECTION 0: Load Required Libraries
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)    # Data manipulation and plotting
  library(geomorph)     # Geometric morphometrics
  library(ggrepel)      # Text labels in plots
  library(prWarp)       # Procrustes warping
  library(abind)        # Array binding
  library(here)         # Path management
  library(ggnewscale)   # Multiple fill scales in ggplot
})


# ============================================================================
# SECTION 1: Helper Functions
# ============================================================================

#' Find convex hull for plotting
#' @param df Data frame with RW1 and RW2 columns
#' @return Data frame with convex hull vertices
find_hull <- function(df) {
  df[chull(df$RW1, df$RW2), ]
}

#' Mix two colors based on weights (for hybrid coloring)
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

#' Remove outlier leaves within each sample
#' Uses robust outlier detection (median + MAD) to flag leaves that deviate
#' too much from the sample mean shape.
#' @param coords Array of coordinates (landmarks x dimensions x specimens)
#' @param sample_id Vector of sample IDs for each specimen
#' @param tol Tolerance multiplier for MAD (default 2.5)
#' @return List with filtered coords and sample_id
rm_leaf_outliers <- function(coords, sample_id, tol = 2.5) {
  keep <- rep(TRUE, dim(coords)[3])  # Logical vector: keep all initially
  for (s in unique(sample_id)) {
    idx <- which(sample_id == s)
    if (length(idx) < 3) next  # Need at least 3 leaves to detect outliers
    # Calculate mean shape for this sample
    mshape <- apply(coords[ , , idx, drop = FALSE], c(1, 2), mean)
    # Calculate distances from mean
    dists <- sapply(idx, \(j) sqrt(sum((coords[ , , j] - mshape)^2)))
    # Robust threshold: median + tolerance * MAD
    cut <- median(dists) + tol * mad(dists)
    keep[idx[dists > cut]] <- FALSE  # Flag outliers
  }
  list(coords = coords[ , , keep, drop = FALSE],
       sample_id = sample_id[keep])
}


# ============================================================================
# SECTION 2: Setup and File Paths
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
tps_file_flowering <- "Hybrid_Zones_LM_fl.tps"
tps_file_short <- "Hybrid_Zones_LM_sh.tps"
admixture_file <- "reassign.tsv"


# ============================================================================
# SECTION 3: Load Admixture and Species Data
# ============================================================================

# Read admixture table and extract species codes
reassign <- read.delim(admixture_file, sep = "\t", stringsAsFactors = FALSE) %>%
  tidyr::separate(SSRs.group_new, into = "group", sep = "_", remove = FALSE) %>%
  mutate(species = gsub("\\d", "", group))  # Remove numbers to get species code


# ============================================================================
# SECTION 4: Main Relative Warp Analysis Function
# ============================================================================

#' Perform Relative Warp analysis on a TPS file
#' 
#' This function:
#' 1. Reads TPS landmark file
#' 2. Filters samples (excludes problematic samples and unknowns)
#' 3. Performs Generalized Procrustes Analysis (GPA) on all leaves
#' 4. Optionally removes outlier leaves within samples
#' 5. Averages coordinates per sample
#' 6. Performs GPA on sample means
#' 7. Runs Relative Warp analysis using Morpho::relWarps
#' 8. Creates morphospace plot with density contours
#' 9. Calculates mean shapes for each species
#'
#' @param tps_file Path to TPS file
#' @param dataset_label Label for plots (e.g., "flowering shoots")
#' @param parental_cols Named vector of colors for parental species
#' @param rm_outliers Logical: remove outlier leaves within samples?
#' @param tol_out Tolerance for outlier detection (MAD multiplier)
#' @return List containing plot, coordinates, metadata, RW results, mean shapes
run_relative_warp_analysis <- function(tps_file,
                                       dataset_label = "dataset",
                                       parental_cols = c(L="#FF5980", M="#28e21b", R="#302790"),
                                       rm_outliers = FALSE,
                                       tol_out = 2) {
  
  # --- 4.1 Read TPS file ---
  raw_coords <- geomorph::readland.tps(tps_file, specID = "ID", warnmsg = FALSE)
  n_landmarks <- dim(raw_coords)[1]
  full_labels <- trimws(dimnames(raw_coords)[[3]])
  sample_id <- sub("_.*", "", full_labels)  # Extract sample ID (before underscore)
  
  # --- 4.2 Filter samples ---
  # Exclude problematic samples and those without species assignment
  excluded_samples <- c("L101", "HMAC611", "L707", "R707", "R905", "HMED210", "H405")
  species_lookup <- reassign$species[match(sample_id, reassign$Name)]
  keep_idx <- !(sample_id %in% excluded_samples | 
                is.na(species_lookup) | 
                species_lookup == "UNKN")
  raw_coords <- raw_coords[ , , keep_idx, drop = FALSE]
  full_labels <- full_labels[keep_idx]
  sample_id <- sample_id[keep_idx]
  
  # --- 4.3 Generalized Procrustes Analysis on all leaves ---
  gpa_all <- gpagen(raw_coords, print.progress = FALSE)
  
  # --- 4.4 Remove within-sample outliers (optional) ---
  if (rm_outliers) {
    filt <- rm_leaf_outliers(gpa_all$coords, sample_id, tol = tol_out)
    gpa_all$coords <- filt$coords
    sample_id <- filt$sample_id
    full_labels <- full_labels[match(sample_id, sample_id)]  # Re-sync labels
  }
  
  # --- 4.5 Average GPA-aligned coordinates per sample ---
  # Multiple leaves per sample are averaged to get one shape per sample
  unique_samples <- unique(sample_id)
  coords_averaged <- array(NA, 
                           dim = c(n_landmarks, 2, length(unique_samples)),
                           dimnames = list(paste0("LM", seq_len(n_landmarks)), 
                                         c("x","y"), 
                                         unique_samples))
  for (i in seq_along(unique_samples)) {
    idx <- which(sample_id == unique_samples[i])
    coords_averaged[ , , i] <- apply(gpa_all$coords[ , , idx, drop = FALSE], 
                                     c(1,2), mean, na.rm = TRUE)
  }
  
  # --- 4.6 GPA on sample means ---
  gpa <- gpagen(coords_averaged, print.progress = FALSE)
  coords_aligned <- gpa$coords  # Final aligned coordinates: landmarks x dim x samples
  
  # --- 4.7 Prepare metadata ---
  metadata <- tibble(sample = unique_samples) %>%
    left_join(reassign %>% dplyr::select(Name, species), 
              by = c("sample" = "Name")) %>%
    mutate(species = if_else(is.na(species), 
                            substring(sample, 1, 1),  # Fallback: use first letter
                            species))
  
  # --- 4.8 Relative Warp analysis ---
  # Morpho::relWarps performs RW analysis (equivalent to PCA on Procrustes coordinates)
  rw_results <- Morpho::relWarps(coords_aligned)
  
  # Extract first 3 relative warp scores
  n_components <- 3
  rw_scores <- rw_results$bescores[, 1:n_components, drop = FALSE]
  colnames(rw_scores) <- paste0("RW", 1:n_components)
  
  # --- 4.9 Prepare scores dataframe with colors ---
  scores_df <- as_tibble(rw_scores, .name_repair = "minimal") %>%
    mutate(sample = unique_samples,
           species = metadata$species[match(unique_samples, metadata$sample)]) %>%
    left_join(reassign %>% dplyr::select(Name, L, M, R), 
              by = c("sample" = "Name")) %>%
    mutate(across(c(L, M, R), as.numeric),
           # Create mixed colors for hybrids based on admixture proportions
           mixed_col = case_when(
             species == "L" ~ parental_cols["L"],
             species == "M" ~ parental_cols["M"],
             species == "R" ~ parental_cols["R"],
             TRUE ~ mix_colors(parental_cols["L"], parental_cols["R"], L, R)
           ))
  
  # --- 4.10 Calculate species centroids ---
  centroids <- scores_df %>%
    filter(species %in% c("L", "M", "R")) %>%
    group_by(species) %>%
    summarise(RW1 = mean(RW1, na.rm = TRUE),
              RW2 = mean(RW2, na.rm = TRUE),
              .groups = "drop")
  
  # --- 4.11 Create density contour scales ---
  n_bins <- 100
  density_scale_L <- scales::alpha("#FF5980", seq(0, 0.4, length.out = n_bins))
  density_scale_M <- scales::alpha("#28e21b", seq(0, 0.4, length.out = n_bins))
  density_scale_R <- scales::alpha("#302790", seq(0, 0.4, length.out = n_bins))
  
  # --- 4.12 Create morphospace plot ---
  morphospace_plot <- ggplot(scores_df, aes(RW1, RW2)) +
    # Density contours for L species
    geom_density_2d_filled(
      data = scores_df %>% filter(species == "L"),
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
      data = scores_df %>% filter(species == "M"),
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
      data = scores_df %>% filter(species == "R"),
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
               size = 3,
               stroke = ifelse(scores_df$species %in% names(parental_cols), 1.2, 0.6)) +
    scale_colour_identity() +
    scale_shape_manual(values = c(L = 15, M = 16, R = 18, MAC = 3, MED = 4, SUB = 8),
                       guide = guide_legend(override.aes = list(size = 3))) +
    labs(title = paste("Morphospace –", dataset_label, "(Morpho relWarps)"),
         x = sprintf("RW1 (%.1f%%)", rw_results$Var$exVar[1] * 100),
         y = sprintf("RW2 (%.1f%%)", rw_results$Var$exVar[2] * 100)) +
    # Species centroids
    geom_point(data = centroids,
               aes(x = RW1, y = RW2, shape = species)) +
    theme_bw()
  
  # --- 4.13 Calculate mean shapes for each species ---
  species_list <- unique(metadata$species)
  mean_shapes_by_species <- lapply(species_list, function(sp) {
    idx <- which(metadata$species == sp)
    apply(coords_aligned[ , , idx, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  })
  names(mean_shapes_by_species) <- species_list
  
  # Grand mean shape (all samples)
  grand_mean_shape <- apply(coords_aligned, c(1, 2), mean, na.rm = TRUE)
  
  # Return results
  list(plot = morphospace_plot,
       coords = coords_averaged,
       metadata = metadata,
       relwarps = rw_results,
       mean_shapes = mean_shapes_by_species,
       grand_mean = grand_mean_shape)
}


# ============================================================================
# SECTION 5: Prepare Data for Combined Analysis
# ============================================================================

# Run RW analysis for both datasets to get coordinates for combining
# (These are needed for the combined analysis but separate outputs are not saved)
results_flowering <- run_relative_warp_analysis(
  tps_file_flowering, 
  dataset_label = "flowering shoots",
  rm_outliers = TRUE,
  tol_out = 2
)

results_short <- run_relative_warp_analysis(
  tps_file_short, 
  dataset_label = "short shoots",
  rm_outliers = TRUE,
  tol_out = 2
)


# ============================================================================
# SECTION 6: Combined Dataset Analysis
# ============================================================================

# Parental species colors
parental_colors <- c(L = "#FF5980", M = "#28e21b", R = "#302790")

#' Pad array with NA values for missing samples
#' Used when combining datasets with different sample sets
pad_array_with_na <- function(arr, sample_union) {
  missing_samples <- setdiff(sample_union, dimnames(arr)[[3]])
  if (length(missing_samples) > 0) {
    na_array <- array(NA,
                     dim = c(dim(arr)[1], 2, length(missing_samples)),
                     dimnames = list(dimnames(arr)[[1]], c("x","y"), missing_samples))
    arr <- abind::abind(arr, na_array, along = 3)
  }
  arr[ , , sample_union, drop = FALSE]  # Reorder to match sample_union
}

# --- 6.1 Combine datasets ---
# Option: use only samples present in both datasets (paired_only = TRUE)
# or use all samples with missing data imputation (paired_only = FALSE)
paired_only <- TRUE

if (paired_only) {
  # Use only samples present in both datasets
  matched_samples <- intersect(
    dimnames(results_flowering$coords)[[3]],
    dimnames(results_short$coords)[[3]]
  )
  coords_flowering_combined <- results_flowering$coords[ , , matched_samples, drop = FALSE]
  coords_short_combined <- results_short$coords[ , , matched_samples, drop = FALSE]
  file_label <- "_paired"
  
  # Rename landmarks to distinguish flowering vs short
  dimnames(coords_flowering_combined)[[1]] <- paste0(dimnames(coords_flowering_combined)[[1]], "_fl")
  dimnames(coords_short_combined)[[1]] <- paste0(dimnames(coords_short_combined)[[1]], "_sh")
  
  # Combine: flowering landmarks (1-13) + short landmarks (14-26)
  combined_coords <- abind::abind(coords_flowering_combined, 
                                  coords_short_combined, 
                                  along = 1)
  
  # Prepare metadata
  combined_metadata <- tibble(sample = dimnames(coords_flowering_combined)[[3]]) %>%
    left_join(reassign %>% dplyr::select(Name, species, L, M, R), 
              by = c("sample" = "Name")) %>%
    mutate(across(c(L, M, R), as.numeric)) %>%
    mutate(species = if_else(is.na(species), 
                            substring(sample, 1, 1), 
                            species))
} else {
  # Use all samples, pad with NA where missing
  all_samples <- union(
    dimnames(results_flowering$coords)[[3]],
    dimnames(results_short$coords)[[3]]
  )
  coords_flowering_combined <- pad_array_with_na(results_flowering$coords, all_samples)
  coords_short_combined <- pad_array_with_na(results_short$coords, all_samples)
  file_label <- "_all"
  
  # Rename landmarks
  dimnames(coords_flowering_combined)[[1]] <- paste0(dimnames(coords_flowering_combined)[[1]], "_fl")
  dimnames(coords_short_combined)[[1]] <- paste0(dimnames(coords_short_combined)[[1]], "_sh")
  
  # Combine
  combined_coords <- abind::abind(coords_flowering_combined, 
                                  coords_short_combined, 
                                  along = 1)
  
  # Prepare metadata
  combined_metadata <- tibble(sample = dimnames(coords_flowering_combined)[[3]]) %>%
    left_join(reassign %>% select(Name, species, L, M, R), 
              by = c("sample" = "Name")) %>%
    mutate(across(c(L, M, R), as.numeric)) %>%
    mutate(species = if_else(is.na(species), 
                            substring(sample, 1, 1), 
                            species))
  
  # Impute missing coordinates using TPS
  combined_coords <- geomorph::estimate.missing(combined_coords, method = "TPS")
}

# --- 6.2 GPA on combined coordinates ---
gpa_combined <- gpagen(combined_coords, print.progress = FALSE)

# --- 6.3 Calculate mean shapes for combined dataset ---
species_list_combined <- unique(combined_metadata$species)
mean_shapes_combined <- lapply(species_list_combined, function(sp) {
  idx <- which(combined_metadata$species == sp)
  apply(gpa_combined$coords[ , , idx, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
})
names(mean_shapes_combined) <- species_list_combined

grand_mean_combined <- apply(gpa_combined$coords, c(1, 2), mean, na.rm = TRUE)

# --- 6.4 Add species centroids to coordinate array for RW analysis ---
# This allows us to visualize species centroids in the morphospace
coords_with_centroids <- abind::abind(
  gpa_combined$coords,
  simplify2array(mean_shapes_combined),
  along = 3
)

dimnames(coords_with_centroids)[[3]] <- c(
  dimnames(gpa_combined$coords)[[3]],
  names(mean_shapes_combined)
)

# --- 6.5 Relative Warp analysis on combined dataset ---
rw_combined_with_centroids <- Morpho::relWarps(coords_with_centroids)
rw_scores_all <- rw_combined_with_centroids$bescores

# Extract centroid scores
centroid_names <- names(mean_shapes_combined)
centroid_scores <- rw_scores_all[centroid_names, , drop = FALSE]
centroid_scores <- as.data.frame(centroid_scores)
centroid_scores <- centroid_scores[, 1:3, drop = FALSE]
colnames(centroid_scores) <- paste0("RW", 1:3)
centroid_scores$species <- centroid_names

# Extract sample scores (original specimens)
original_sample_names <- dimnames(gpa_combined$coords)[[3]]
rw_combined_scores <- rw_scores_all[original_sample_names, , drop = FALSE]
rw_combined_scores <- rw_combined_scores[, 1:3, drop = FALSE]
colnames(rw_combined_scores) <- paste0("RW", 1:3)

# Prepare scores dataframe
scores_combined <- as_tibble(rw_combined_scores, .name_repair = "minimal") %>%
  mutate(sample = rownames(rw_combined_scores)) %>%
  left_join(combined_metadata, by = "sample") %>%
  mutate(mixed_col = case_when(
    species == "L" ~ parental_colors["L"],
    species == "M" ~ parental_colors["M"],
    species == "R" ~ parental_colors["R"],
    TRUE ~ mix_colors(parental_colors["L"], parental_colors["R"], L, R)
  ))

# --- 6.6 Create combined morphospace plot ---
n_bins <- 100
density_scale_L <- scales::alpha("#FF5980", seq(0, 0.4, length.out = n_bins))
density_scale_M <- scales::alpha("#28e21b", seq(0, 0.4, length.out = n_bins))
density_scale_R <- scales::alpha("#302790", seq(0, 0.4, length.out = n_bins))

plot_combined <- ggplot(scores_combined, aes(RW1, RW2)) +
  geom_density_2d_filled(
    data = scores_combined %>% filter(species == "L"),
    aes(fill = after_stat(level)),
    contour_var = "ndensity",
    bins = n_bins,
    show.legend = FALSE,
    color = NA,
    linetype = 0
  ) +
  scale_fill_manual(values = density_scale_L) +
  new_scale_fill() +
  geom_density_2d_filled(
    data = scores_combined %>% filter(species == "M"),
    aes(fill = after_stat(level)),
    contour_var = "ndensity",
    bins = n_bins,
    show.legend = FALSE,
    color = NA,
    linetype = 0
  ) +
  scale_fill_manual(values = density_scale_M) +
  new_scale_fill() +
  geom_density_2d_filled(
    data = scores_combined %>% filter(species == "R"),
    aes(fill = after_stat(level)),
    contour_var = "ndensity",
    bins = n_bins,
    show.legend = FALSE,
    color = NA,
    linetype = 0
  ) +
  scale_fill_manual(values = density_scale_R) +
  geom_point(aes(colour = mixed_col, shape = species), 
             size = 4,
             stroke = ifelse(scores_combined$species %in% names(parental_colors), 1.2, 0.6)) +
  scale_colour_identity() +
  scale_shape_manual(values = c(L = 15, M = 16, R = 18, MAC = 3, MED = 4, SUB = 8),
                     guide = guide_legend(override.aes = list(size = 4))) +
  labs(title = paste("Morphospace – combined (Morpho relWarps)"),
       x = sprintf("RW1 (%.1f%%)", rw_combined_with_centroids$Var$exVar[1] * 100),
       y = sprintf("RW2 (%.1f%%)", rw_combined_with_centroids$Var$exVar[2] * 100)) +
  geom_point(data = centroid_scores,
             aes(x = RW1, y = RW2, shape = species)) +
  xlim(-0.2, 0.2) +
  ylim(-0.075, 0.125) +
  theme_bw()

ggsave(paste0("relative_warp_morphospace_MorphoRW_combined", file_label, ".pdf"),
       plot = plot_combined, width = 12, height = 6)

# --- 6.7 Generate TPS grids for combined dataset (tiled) ---
# Define landmark links for visualization
landmark_links <- rbind(cbind(1:12, 2:13), c(13, 2))

#' Calculate optimal grid dimensions for tiled layout
#' @param n Number of plots
#' @return Vector with nrow, ncol
calculate_grid_dims <- function(n) {
  if (n <= 0) return(c(1, 1))
  ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n / ncol)
  c(nrow, ncol)
}
n_species_combined <- length(species_list_combined)
grid_dims_combined <- calculate_grid_dims(n_species_combined)

# Flowering shoot portion (landmarks 1-13)
pdf(file = paste0("TPS_grid_combined_flowering", file_label, "_tiled.pdf"), 
    width = 4 * grid_dims_combined[2], 
    height = 4 * grid_dims_combined[1])
par(mfrow = grid_dims_combined, mar = c(2, 2, 2, 2))
for (species in species_list_combined) {
  plotRefToTarget(
    grand_mean_combined[1:13, ],
    mean_shapes_combined[[species]][1:13, ],
    method = "TPS",
    gridPars = gridPar(grid.lwd = 0.1, link.lwd = 0.5),
    links = landmark_links
  )
  title(main = species, cex.main = 1.2)
}
dev.off()

# Short shoot portion (landmarks 14-26)
pdf(file = paste0("TPS_grid_combined_short", file_label, "_tiled.pdf"), 
    width = 4 * grid_dims_combined[2], 
    height = 4 * grid_dims_combined[1])
par(mfrow = grid_dims_combined, mar = c(2, 2, 2, 2))
for (species in species_list_combined) {
  plotRefToTarget(
    grand_mean_combined[14:26, ],
    mean_shapes_combined[[species]][14:26, ],
    method = "TPS",
    gridPars = gridPar(grid.lwd = 0.1, link.lwd = 0.5),
    links = landmark_links
  )
  title(main = species, cex.main = 1.2)
}
dev.off()


# ============================================================================
# SECTION 7: Fruit PCA Analysis
# ============================================================================

library(FactoMineR)
library(factoextra)

# --- 7.1 Load and prepare fruit data ---
fruits_data <- read.delim("fruits_avg.tsv", sep = "\t", stringsAsFactors = FALSE)

# Filter excluded samples
excluded_samples_fruits <- c("L101", "HMAC611", "L707", "R707", "R905", 
                             "HMED210", "H405", "R904")
fruits_data <- fruits_data %>% 
  filter(!ID %in% excluded_samples_fruits, !Taxon %in% c("UNKN", NA))

# Add admixture data
fruits_data <- fruits_data %>%
  left_join(reassign %>% dplyr::select(Name, species, L, M, R), 
            by = c("ID" = "Name"))

# Calculate species means in original space
species_means_original <- fruits_data %>%
  group_by(species) %>%
  summarise(across(DR:P, \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("species")

# --- 7.2 Prepare colors for hybrids ---
parental_colors <- c(L = "#FF5980", M = "#28e21b", R = "#302790")
fruits_data <- fruits_data %>%
  mutate(across(c(L.y, M, R), as.numeric),
         mixed_col = case_when(
           species == "L" ~ parental_colors["L"],
           species == "M" ~ parental_colors["M"],
           species == "R" ~ parental_colors["R"],
           species == "MAC" ~ mix_colors(parental_colors["L"], parental_colors["R"], L.y, R),
           species == "MED" ~ mix_colors(parental_colors["L"], parental_colors["M"], L.y, M),
           species == "SUB" ~ mix_colors(parental_colors["R"], parental_colors["M"], R, M),
           TRUE ~ NA_character_
         ))

# --- 7.3 Perform PCA ---
# Select numeric variables and scale
fruit_variables <- fruits_data %>% 
  dplyr::select(where(is.numeric)) %>% 
  dplyr::select(-Population, -L.y, -M, -R) %>% 
  mutate(across(everything(), scale))

# Calculate species means in scaled space
species_means_scaled <- fruit_variables %>%
  mutate(species = fruits_data$species) %>%
  group_by(species) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("species")

# Add row names for matching
rownames(fruit_variables) <- fruits_data$ID

# Combine samples and species means for PCA
fruit_data_augmented <- rbind(fruit_variables, species_means_scaled)

# Perform PCA
fruit_pca <- PCA(fruit_data_augmented, scale.unit = FALSE, graph = FALSE)

# --- 7.4 Extract scores ---
all_pca_scores <- as.data.frame(fruit_pca$ind$coord)

# Extract centroid scores (species means)
centroid_names_fruits <- rownames(species_means_scaled)
centroids_pca <- all_pca_scores[centroid_names_fruits, 1:2, drop = FALSE]
colnames(centroids_pca) <- c("PC1", "PC2")
centroids_pca$species <- centroid_names_fruits

# Extract sample scores
scores_pca_fruits <- all_pca_scores[fruits_data$ID, 1:2, drop = FALSE]
colnames(scores_pca_fruits) <- c("PC1", "PC2")
scores_pca_fruits$ID <- fruits_data$ID
scores_pca_fruits$species <- fruits_data$species
scores_pca_fruits$mixed_col <- fruits_data$mixed_col

# --- 7.5 Create PCA plot ---
n_bins <- 100
density_scale_L <- scales::alpha("#FF5980", seq(0, 0.4, length.out = n_bins))
density_scale_M <- scales::alpha("#28e21b", seq(0, 0.4, length.out = n_bins))
density_scale_R <- scales::alpha("#302790", seq(0, 0.4, length.out = n_bins))

plot_fruits_pca <- ggplot(scores_pca_fruits, aes(PC1, PC2)) +
  geom_density_2d_filled(
    data = scores_pca_fruits %>% filter(species == "L"),
    aes(fill = after_stat(level)),
    contour_var = "ndensity", 
    bins = n_bins,
    show.legend = FALSE, 
    color = NA, 
    linetype = 0
  ) +
  scale_fill_manual(values = density_scale_L) +
  new_scale_fill() +
  geom_density_2d_filled(
    data = scores_pca_fruits %>% filter(species == "M"),
    aes(fill = after_stat(level)),
    contour_var = "ndensity", 
    bins = n_bins,
    show.legend = FALSE, 
    color = NA, 
    linetype = 0
  ) +
  scale_fill_manual(values = density_scale_M) +
  new_scale_fill() +
  geom_density_2d_filled(
    data = scores_pca_fruits %>% filter(species == "R"),
    aes(fill = after_stat(level)),
    contour_var = "ndensity", 
    bins = n_bins,
    show.legend = FALSE, 
    color = NA, 
    linetype = 0
  ) +
  scale_fill_manual(values = density_scale_R) +
  geom_point(aes(colour = mixed_col, shape = species), 
             size = 2, 
             stroke = 1) +
  geom_point(data = centroids_pca,
             aes(x = PC1, y = PC2, shape = species)) +
  scale_colour_identity() +
  scale_shape_manual(values = c(L = 15, M = 16, R = 18, MAC = 3, MED = 4, SUB = 8),
                     guide = guide_legend(override.aes = list(size = 2))) +
  labs(title = "Fruits PCA with hybrid and parental coloring + density contours + centroids",
       x = paste0("PC1 (", round(fruit_pca$eig[1, 2], 1), "%)"),
       y = paste0("PC2 (", round(fruit_pca$eig[2, 2], 1), "%)")) +
  xlim(-6, 5.75) +
  ylim(-3.75, 3.25) +
  theme_bw() +
  theme(legend.position = "right")

ggsave("fruits_PCA_cent.pdf", plot = plot_fruits_pca, width = 8.4, height = 6)

