###############################################################################
# Script: 02_canonical_discriminant_analysis.R
# 
# Purpose:
#   Performs Canonical Discriminant Analysis (CDA) on morphometric data from
#   flowering shoots, short shoots, fruits, and combined morphology datasets.
#   Removes highly correlated variables, performs CDA, and creates plots with
#   hybrid coloring. Also includes LDA classification accuracy assessment.
#
# Input Files Required:
#   - Hybrid_Zones_LM_fl_procrustes.tsv: Procrustes-aligned coordinates for flowering shoots
#   - Hybrid_Zones_LM_sh_procrustes.tsv: Procrustes-aligned coordinates for short shoots
#   - all_measured_fruits.tsv: Raw fruit measurements
#   - completed_morphology.tsv: Combined morphology dataset
#   - reassign.tsv: Admixture and species assignment table
#
# Output Files Generated:
#   - morphometrics_avg_wo_hybrids.pdf: CDA plots for all datasets (without hybrids)
#   - loadings_avg_wo_hybrids.pdf: CDA loadings plots for all datasets
#   - hybrids.pdf: CDA plots showing hybrid groups separately
#   - flowering_shoots_avg.tsv, short_shoots_avg.tsv, fruits_avg.tsv: Averaged datasets
#   - loadings_*.csv: CDA loadings for each dataset
#
# Dependencies (with versions):
#   - tidyr (1.3.1), dplyr (1.1.4), mice (3.18.0), VIM (6.2.2)
#   - ggpubr (0.6.1), FactoMineR (2.12), MorphoTools2 (1.0.2.1)
#   - ggrepel (0.9.6), MASS (7.3.65) - for LDA
#
# Author: Roman Ufimov
# Date: July 2025
# 
# Usage:
#   Place all input files in the same directory as this script, then run:
#   source("02_canonical_discriminant_analysis.R")
#   or
#   Rscript 02_canonical_discriminant_analysis.R
#
###############################################################################


# ============================================================================
# SECTION 0: Load Required Libraries
# ============================================================================

library(tidyr)
library(dplyr)
library(mice)          # Multiple imputation (loaded but may not be used)
library(VIM)           # Visualization and imputation (loaded but may not be used)
library(ggpubr)        # Publication-ready plots
library(FactoMineR)    # Multivariate analysis
library(MorphoTools2)  # Morphometric analysis tools
library(ggrepel)       # Text labels in plots
library(MASS)          # Linear Discriminant Analysis


# ============================================================================
# SECTION 1: Helper Functions
# ============================================================================

#' Find convex hull vertices for a group
#' Used for drawing convex hulls around species groups in CDA plots
#' @param df Data frame with Can1 and Can2 columns (canonical axes)
#' @return Data frame with convex hull vertices
find_hull <- function(df) {
  df[chull(df$Can1, df$Can2), ]
}

#' Mix two colors based on admixture proportions
#' Used to create gradient colors for hybrid individuals
#' @param color1 First color (hex code)
#' @param color2 Second color (hex code)
#' @param percentage1 Admixture proportion for first color
#' @param percentage2 Admixture proportion for second color
#' @return Mixed color as hex code
mix_colors <- function(color1, color2, percentage1, percentage2) {
  if (!is.na(percentage1)) {
    # Convert colors to RGB (0-1 scale)
    col1_rgb <- col2rgb(color1) / 255
    col2_rgb <- col2rgb(color2) / 255
    
    # Calculate weighted average
    mixed_rgb <- (percentage1 * col1_rgb + percentage2 * col2_rgb)
    
    # Convert back to hex
    rgb(mixed_rgb[1], mixed_rgb[2], mixed_rgb[3])
  }
}

#' Create dummy data frame for missing species
#' Used to maintain consistent factor levels in plots
#' @param missing_species Vector of species names to create dummy entries for
#' @return Data frame with NA values for Can1, Can2, and mixed_color
create_dummy_data <- function(missing_species) {
  data.frame(
    Can1 = NA,
    Can2 = NA,
    sp = missing_species,
    mixed_color = NA
  )
}


# ============================================================================
# SECTION 2: Setup and File Paths
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
procrustes_file_flowering <- "Hybrid_Zones_LM_fl_procrustes.tsv"
procrustes_file_short <- "Hybrid_Zones_LM_sh_procrustes.tsv"
fruits_file_raw <- "all_measured_fruits.tsv"
morphology_file_combined <- "completed_morphology.tsv"
admixture_file <- "reassign.tsv"


# ============================================================================
# SECTION 3: Process Flowering Shoots Data
# ============================================================================

# --- 3.1 Load and prepare flowering shoots data ---
flowering_shoots_raw <- read.csv(procrustes_file_flowering, sep = '\t')

# Separate sample ID from leaf number
flowering_shoots_raw <- flowering_shoots_raw %>% 
  separate(Id, into = c('Sample', 'leaf'), sep = '_', remove = FALSE)

# Rename columns to standard format
colnames(flowering_shoots_raw) <- c('Id', 'Sample', 'leaf', 
                                    'LM1_x', 'LM1_y', 'LM2_x', 'LM2_y', 
                                    'LM3_x', 'LM3_y', 'LM4_x', 'LM4_y', 
                                    'LM5_x', 'LM5_y', 'LM6_x', 'LM6_y',
                                    'LM7_x', 'LM7_y', 'LM8_x', 'LM8_y', 
                                    'LM9_x', 'LM9_y', 'LM10_x', 'LM10_y',
                                    'LM11_x', 'LM11_y', 'LM12_x', 'LM12_y', 
                                    'LM13_x', 'LM13_y')

# Flip y-coordinates (standard morphometric convention)
flowering_shoots_raw <- flowering_shoots_raw %>% 
  mutate(across(ends_with("_y"), ~ . * -1))

# --- 3.2 Load admixture data ---
reassign <- read.csv(admixture_file, sep = '\t')
reassign$pop <- gsub(".*?(\\d).*", "\\1", reassign$Name)  # Extract population number
reassign <- reassign %>% 
  separate(SSRs.group_new, into = c('group'), sep = '_', remove = FALSE)
reassign$species <- gsub("\\d", "", reassign$group)  # Extract species code

# --- 3.3 Filter and average by sample ---
excluded_samples <- c('L101', 'HMAC611', 'L707', 'R707', 'R905', 'HMED210', 'H405')
flowering_shoots_processed <- flowering_shoots_raw %>% 
  filter(!Sample %in% excluded_samples) %>%
  group_by(Sample) %>%
  summarize(across(LM1_x:LM13_y, mean, na.rm = TRUE))

# --- 3.4 Add species information ---
flowering_shoots_processed$Genetic_group_detailed <- reassign$SSRs.group_new[
  match(flowering_shoots_processed$Sample, reassign$Name)
]
flowering_shoots_processed <- flowering_shoots_processed %>% 
  separate(Genetic_group_detailed, into = c('Genetic_group'), sep = '_', remove = FALSE)
flowering_shoots_processed$species <- gsub("\\d", "", flowering_shoots_processed$Genetic_group)
flowering_shoots_processed$pop <- gsub(".*?(\\d).*", "\\1", flowering_shoots_processed$Sample)
flowering_shoots_processed <- flowering_shoots_processed %>% 
  filter(species != 'UNKN')

# --- 3.5 Save averaged dataset for MorphoTools2 ---
flowering_shoots_for_morpho <- flowering_shoots_processed %>% 
  dplyr::select(1, 31, 30, 2:27)
colnames(flowering_shoots_for_morpho)[1:3] <- c('ID', 'Population', 'Taxon')
write.table(flowering_shoots_for_morpho, 
            'flowering_shoots_avg.tsv', 
            sep = '\t', 
            row.names = FALSE)

# --- 3.6 Load into MorphoTools2 format ---
flowering_shoots_morpho <- read.morphodata('flowering_shoots_avg.tsv', sep = '\t')
flowering_shoots_data <- flowering_shoots_morpho$data

# --- 3.7 Remove highly correlated variables ---
# Calculate Spearman correlations
correlations_spearman <- cormat(flowering_shoots_morpho, method = "spearman")
correlation_threshold <- 0.90

# Prepare correlation matrix
rownames(correlations_spearman) <- correlations_spearman$Spearman
correlations_matrix <- correlations_spearman[, -1]

# Find highly correlated pairs
highly_correlated_pairs <- which(abs(correlations_matrix) > correlation_threshold, 
                                 arr.ind = TRUE)
# Remove diagonal (self-correlations)
highly_correlated_pairs <- highly_correlated_pairs[
  highly_correlated_pairs[, 1] != highly_correlated_pairs[, 2], 
]

# Convert to data frame and remove duplicates (i,j and j,i)
highly_correlated_df <- as.data.frame(highly_correlated_pairs)
highly_correlated_df <- highly_correlated_df[
  highly_correlated_df[,1] < highly_correlated_df[,2], 
]

# Calculate variances to decide which variable to remove
variances <- apply(flowering_shoots_data, 2, var)
variable_names <- colnames(flowering_shoots_data)

# Function to select variable with lower variance for removal
select_low_variance <- function(i, j, variances) {
  if (variances[i] < variances[j]) {
    return(i)  # Remove variable with lower variance
  } else {
    return(j)
  }
}

# Determine which variables to remove
variables_to_remove_indices <- apply(highly_correlated_df, 1, 
                                     function(row) select_low_variance(row[1], row[2], variances))
variables_to_remove_indices <- unique(variables_to_remove_indices)
variables_to_remove_names <- variable_names[variables_to_remove_indices]

# Remove highly correlated variables
flowering_shoots_morpho <- removeCharacter(flowering_shoots_morpho, 
                                          variables_to_remove_names)

# --- 3.8 Perform Canonical Discriminant Analysis ---
# Hybrids (MAC, MED, SUB) are set as passive samples (not used in CDA calculation)
cda_results_flowering <- cda.calc(flowering_shoots_morpho, 
                                  passiveSamples = c('MAC', 'MED', 'SUB'))

# Extract canonical scores
cda_scores_flowering <- cda_results_flowering$obj$scores
cda_scores_flowering$sample <- rownames(cda_scores_flowering)
cda_scores_flowering$species <- flowering_shoots_processed$species[
  match(flowering_shoots_processed$Sample, cda_scores_flowering$sample)
]

# Remove hybrids from plotting (they were passive in CDA)
cda_scores_flowering <- cda_scores_flowering %>% 
  filter(!species %in% c('MAC', 'MED', 'SUB'))

# --- 3.9 Extract loadings (canonical structure) ---
loadings_flowering <- data.frame(cda_results_flowering$totalCanonicalStructure)
loadings_flowering$Character <- rownames(loadings_flowering)
max_range_flowering <- max(abs(loadings_flowering$Can1), 
                          abs(loadings_flowering$Can2))

# --- 3.10 Calculate convex hulls for plotting ---
hulls_flowering <- cda_scores_flowering %>%
  group_by(species) %>%
  do(data.frame(find_hull(.)))

# --- 3.11 Create loadings plot ---
plot_loadings_flowering <- ggplot(loadings_flowering, aes(x = Can1, y = Can2)) +
  geom_segment(aes(x = 0, y = 0, xend = Can1, yend = Can2), 
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black") +
  geom_text_repel(aes(label = Character), 
                  size = 3, 
                  max.overlaps = 20) +
  labs(x = paste('Can1 ', 
                 round(cda_results_flowering$eigenvaluesAsPercentages[[1]] * 100, 
                       digits = 2), '%'), 
       y = paste('Can2 ', 
                 round(cda_results_flowering$eigenvaluesAsPercentages[[2]] * 100, 
                       digits = 2), '%')) +
  theme_bw() +
  xlim(-max_range_flowering, max_range_flowering) +
  ylim(-max_range_flowering, max_range_flowering)

# --- 3.12 Create CDA plot ---
plot_cda_flowering <- ggplot(cda_scores_flowering, aes(x = Can1, y = Can2)) + 
  geom_point(aes(color = species, shape = species), 
             stroke = 1, 
             size = 3) +
  labs(x = paste('Can1 ', 
                 round(cda_results_flowering$eigenvaluesAsPercentages[[1]] * 100, 
                       digits = 2), '%'), 
       y = paste('Can2 ', 
                 round(cda_results_flowering$eigenvaluesAsPercentages[[2]] * 100, 
                       digits = 2), '%')) +
  theme_bw() +
  # Add convex hulls
  geom_polygon(data = hulls_flowering, 
               aes(x = Can1, y = Can2, group = species, fill = species), 
               alpha = 0.2, 
               color = NA, 
               show.legend = FALSE) +
  scale_color_manual(name = "Species", 
                     values = c("L" = "#FF5980", 
                               "M" = "#28e21b", 
                               "R" = "#302790")) +
  scale_shape_manual(name = "Species", 
                    values = c("L" = 15,   # Square
                              "M" = 16,   # Circle
                              "R" = 18))  # Diamond


# ============================================================================
# SECTION 4: Process Short Shoots Data
# ============================================================================
# (Same procedure as flowering shoots)

short_shoots_raw <- read.csv(procrustes_file_short, sep = '\t')
short_shoots_raw <- short_shoots_raw %>% 
  separate(Id, into = c('Sample', 'leaf'), sep = '_', remove = FALSE)
colnames(short_shoots_raw) <- c('Id', 'Sample', 'leaf', 
                                'LM1_x', 'LM1_y', 'LM2_x', 'LM2_y', 
                                'LM3_x', 'LM3_y', 'LM4_x', 'LM4_y', 
                                'LM5_x', 'LM5_y', 'LM6_x', 'LM6_y',
                                'LM7_x', 'LM7_y', 'LM8_x', 'LM8_y', 
                                'LM9_x', 'LM9_y', 'LM10_x', 'LM10_y',
                                'LM11_x', 'LM11_y', 'LM12_x', 'LM12_y', 
                                'LM13_x', 'LM13_y')
short_shoots_raw <- short_shoots_raw %>% 
  mutate(across(ends_with("_y"), ~ . * -1))

reassign <- read.csv(admixture_file, sep = '\t')
reassign$pop <- gsub(".*?(\\d).*", "\\1", reassign$Name)
reassign <- reassign %>% 
  separate(SSRs.group_new, into = c('group'), sep = '_', remove = FALSE)
reassign$species <- gsub("\\d", "", reassign$group)

short_shoots_processed <- short_shoots_raw %>% 
  filter(!Sample %in% excluded_samples) %>%
  group_by(Sample) %>%
  summarize(across(LM1_x:LM13_y, mean, na.rm = TRUE))

short_shoots_processed$Genetic_group_detailed <- reassign$SSRs.group_new[
  match(short_shoots_processed$Sample, reassign$Name)
]
short_shoots_processed <- short_shoots_processed %>% 
  separate(Genetic_group_detailed, into = c('Genetic_group'), sep = '_', remove = FALSE)
short_shoots_processed$species <- gsub("\\d", "", short_shoots_processed$Genetic_group)
short_shoots_processed$pop <- gsub(".*?(\\d).*", "\\1", short_shoots_processed$Sample)
short_shoots_processed <- short_shoots_processed %>% 
  filter(species != 'UNKN')

short_shoots_for_morpho <- short_shoots_processed %>% 
  dplyr::select(1, 31, 30, 2:27)
colnames(short_shoots_for_morpho)[1:3] <- c('ID', 'Population', 'Taxon')
write.table(short_shoots_for_morpho, 
            'short_shoots_avg.tsv', 
            sep = '\t', 
            row.names = FALSE)

short_shoots_morpho <- read.morphodata('short_shoots_avg.tsv', sep = '\t')
short_shoots_data <- short_shoots_morpho$data

# Remove highly correlated variables
correlations_spearman <- cormat(short_shoots_morpho, method = "spearman")
rownames(correlations_spearman) <- correlations_spearman$Spearman
correlations_matrix <- correlations_spearman[, -1]
highly_correlated_pairs <- which(abs(correlations_matrix) > correlation_threshold, 
                                 arr.ind = TRUE)
highly_correlated_pairs <- highly_correlated_pairs[
  highly_correlated_pairs[, 1] != highly_correlated_pairs[, 2], 
]
highly_correlated_df <- as.data.frame(highly_correlated_pairs)
highly_correlated_df <- highly_correlated_df[
  highly_correlated_df[,1] < highly_correlated_df[,2], 
]
variances <- apply(short_shoots_data, 2, var)
variable_names <- colnames(short_shoots_data)
variables_to_remove_indices <- apply(highly_correlated_df, 1, 
                                     function(row) select_low_variance(row[1], row[2], variances))
variables_to_remove_indices <- unique(variables_to_remove_indices)
variables_to_remove_names <- variable_names[variables_to_remove_indices]
short_shoots_morpho <- removeCharacter(short_shoots_morpho, 
                                       variables_to_remove_names)

# CDA
cda_results_short <- cda.calc(short_shoots_morpho, 
                              passiveSamples = c('MAC', 'MED', 'SUB'))
cda_scores_short <- cda_results_short$obj$scores
cda_scores_short$sample <- rownames(cda_scores_short)
cda_scores_short$species <- short_shoots_processed$species[
  match(short_shoots_processed$Sample, cda_scores_short$sample)
]
cda_scores_short <- cda_scores_short %>% 
  filter(!species %in% c('MAC', 'MED', 'SUB'))

loadings_short <- data.frame(cda_results_short$totalCanonicalStructure)
loadings_short$Character <- rownames(loadings_short)
max_range_short <- max(abs(loadings_short$Can1), abs(loadings_short$Can2))

hulls_short <- cda_scores_short %>%
  group_by(species) %>%
  do(data.frame(find_hull(.)))

plot_loadings_short <- ggplot(loadings_short, aes(x = Can1, y = Can2)) +
  geom_segment(aes(x = 0, y = 0, xend = Can1, yend = Can2), 
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black") +
  geom_text_repel(aes(label = Character), 
                  size = 3, 
                  max.overlaps = 20) +
  labs(x = paste('Can1 ', 
                 round(cda_results_short$eigenvaluesAsPercentages[[1]] * 100, 
                       digits = 2), '%'), 
       y = paste('Can2 ', 
                 round(cda_results_short$eigenvaluesAsPercentages[[2]] * 100, 
                       digits = 2), '%')) +
  theme_bw() +
  xlim(-max_range_short, max_range_short) +
  ylim(-max_range_short, max_range_short)

plot_cda_short <- ggplot(cda_scores_short, aes(x = Can1, y = Can2)) + 
  geom_point(aes(color = species, shape = species), 
             stroke = 1, 
             size = 3) +
  labs(x = paste('Can1 ', 
                 round(cda_results_short$eigenvaluesAsPercentages[[1]] * 100, 
                       digits = 2), '%'), 
       y = paste('Can2 ', 
                 round(cda_results_short$eigenvaluesAsPercentages[[2]] * 100, 
                       digits = 2), '%')) +
  theme_bw() +
  geom_polygon(data = hulls_short, 
               aes(x = Can1, y = Can2, group = species, fill = species), 
               alpha = 0.2, 
               color = NA, 
               show.legend = FALSE) +
  scale_color_manual(name = "Species", 
                     values = c("L" = "#FF5980", 
                               "M" = "#28e21b", 
                               "R" = "#302790")) +
  scale_shape_manual(name = "Species", 
                    values = c("L" = 15, 
                              "M" = 16, 
                              "R" = 18))


# ============================================================================
# SECTION 5: Process Fruits Data
# ============================================================================

fruits_raw <- read.csv(fruits_file_raw, sep = '\t')
fruits_raw <- fruits_raw %>% 
  separate(Id, into = c('Sample', 'fruit'), sep = '_', remove = FALSE)
fruits_raw <- na.omit(fruits_raw)

reassign <- read.csv(admixture_file, sep = '\t')
reassign$pop <- gsub(".*?(\\d).*", "\\1", reassign$Name)
reassign <- reassign %>% 
  separate(SSRs.group_new, into = c('group'), sep = '_', remove = FALSE)
reassign$species <- gsub("\\d", "", reassign$group)

fruits_processed <- fruits_raw %>% 
  filter(!Sample %in% excluded_samples) %>%
  group_by(Sample) %>%
  summarize(across(DR:P, mean, na.rm = TRUE))

fruits_processed$Genetic_group_detailed <- reassign$SSRs.group_new[
  match(fruits_processed$Sample, reassign$Name)
]
fruits_processed <- fruits_processed %>% 
  separate(Genetic_group_detailed, into = c('Genetic_group'), sep = '_', remove = FALSE)
fruits_processed$species <- gsub("\\d", "", fruits_processed$Genetic_group)
fruits_processed$pop <- gsub(".*?(\\d).*", "\\1", fruits_processed$Sample)
fruits_processed <- fruits_processed %>% 
  filter(species != 'UNKN')

fruits_for_morpho <- fruits_processed %>% 
  dplyr::select(1, 12, 11, 2:8)
colnames(fruits_for_morpho)[1:3] <- c('ID', 'Population', 'Taxon')
write.table(fruits_for_morpho, 
            'fruits_avg.tsv', 
            sep = '\t', 
            row.names = FALSE)

fruits_morpho <- read.morphodata('fruits_avg.tsv', sep = '\t')

# CDA for fruits (no correlation removal needed - fewer variables)
cda_results_fruits <- cda.calc(fruits_morpho, 
                              passiveSamples = c('MAC', 'MED', 'SUB'))
cda_scores_fruits <- cda_results_fruits$obj$scores
cda_scores_fruits$sample <- rownames(cda_scores_fruits)
cda_scores_fruits$species <- fruits_processed$species[
  match(fruits_processed$Sample, cda_scores_fruits$sample)
]
cda_scores_fruits <- cda_scores_fruits %>% 
  filter(!species %in% c('MAC', 'MED', 'SUB'))

loadings_fruits <- data.frame(cda_results_fruits$totalCanonicalStructure)
loadings_fruits$Character <- rownames(loadings_fruits)
max_range_fruits <- max(abs(loadings_fruits$Can1), abs(loadings_fruits$Can2))

hulls_fruits <- cda_scores_fruits %>%
  group_by(species) %>%
  do(data.frame(find_hull(.)))

plot_loadings_fruits <- ggplot(loadings_fruits, aes(x = Can1, y = Can2)) +
  geom_segment(aes(x = 0, y = 0, xend = Can1, yend = Can2), 
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black") +
  geom_text_repel(aes(label = Character), 
                  size = 3, 
                  max.overlaps = 20) +
  labs(x = paste('Can1 ', 
                 round(cda_results_fruits$eigenvaluesAsPercentages[[1]] * 100, 
                       digits = 2), '%'), 
       y = paste('Can2 ', 
                 round(cda_results_fruits$eigenvaluesAsPercentages[[2]] * 100, 
                       digits = 2), '%')) +
  theme_bw() +
  xlim(-max_range_fruits, max_range_fruits) +
  ylim(-max_range_fruits, max_range_fruits)

plot_cda_fruits <- ggplot(cda_scores_fruits, aes(x = Can1, y = Can2)) + 
  geom_point(aes(color = species, shape = species), 
             stroke = 1, 
             size = 3) +
  labs(x = paste('Can1 ', 
                 round(cda_results_fruits$eigenvaluesAsPercentages[[1]] * 100, 
                       digits = 2), '%'), 
       y = paste('Can2 ', 
                 round(cda_results_fruits$eigenvaluesAsPercentages[[2]] * 100, 
                       digits = 2), '%')) +
  theme_bw() +
  geom_polygon(data = hulls_fruits, 
               aes(x = Can1, y = Can2, group = species, fill = species), 
               alpha = 0.2, 
               color = NA, 
               show.legend = FALSE) +
  scale_color_manual(name = "Species", 
                     values = c("L" = "#FF5980", 
                               "M" = "#28e21b", 
                               "R" = "#302790")) +
  scale_shape_manual(name = "Species", 
                    values = c("L" = 15, 
                              "M" = 16, 
                              "R" = 18))


# ============================================================================
# SECTION 6: Process Combined Morphology Dataset
# ============================================================================

morphology_combined <- read.morphodata(morphology_file_combined, sep = '\t')
morphology_combined_data <- morphology_combined$data

# Remove highly correlated variables (same procedure as before)
correlations_spearman <- cormat(morphology_combined, method = "spearman")
rownames(correlations_spearman) <- correlations_spearman$Spearman
correlations_matrix <- correlations_spearman[, -1]
highly_correlated_pairs <- which(abs(correlations_matrix) > correlation_threshold, 
                                 arr.ind = TRUE)
highly_correlated_pairs <- highly_correlated_pairs[
  highly_correlated_pairs[, 1] != highly_correlated_pairs[, 2], 
]
highly_correlated_df <- as.data.frame(highly_correlated_pairs)
highly_correlated_df <- highly_correlated_df[
  highly_correlated_df[,1] < highly_correlated_df[,2], 
]
variances <- apply(morphology_combined_data, 2, var)
variable_names <- colnames(morphology_combined_data)
variables_to_remove_indices <- apply(highly_correlated_df, 1, 
                                     function(row) select_low_variance(row[1], row[2], variances))
variables_to_remove_indices <- unique(variables_to_remove_indices)
variables_to_remove_names <- variable_names[variables_to_remove_indices]
morphology_combined <- removeCharacter(morphology_combined, 
                                      variables_to_remove_names)

# CDA
cda_results_combined <- cda.calc(morphology_combined, 
                                passiveSamples = c('MAC', 'MED', 'SUB'))
cda_scores_combined <- cda_results_combined$obj$scores
cda_scores_combined$sample <- rownames(cda_scores_combined)
cda_scores_combined$species <- reassign$species[
  match(cda_scores_combined$sample, reassign$Name)
]
cda_scores_combined <- cda_scores_combined %>% 
  filter(!species %in% c('MAC', 'MED', 'SUB'))

loadings_combined <- data.frame(cda_results_combined$totalCanonicalStructure)
loadings_combined$Character <- rownames(loadings_combined)
max_range_combined <- max(abs(loadings_combined$Can1), 
                         abs(loadings_combined$Can2))

hulls_combined <- cda_scores_combined %>%
  group_by(species) %>%
  do(data.frame(find_hull(.)))

plot_loadings_combined <- ggplot(loadings_combined, aes(x = Can1, y = Can2)) +
  geom_segment(aes(x = 0, y = 0, xend = Can1, yend = Can2), 
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black") +
  geom_text_repel(aes(label = Character), 
                  size = 3, 
                  max.overlaps = 20) +
  labs(x = paste('Can1 ', 
                 round(cda_results_combined$eigenvaluesAsPercentages[[1]] * 100, 
                       digits = 2), '%'), 
       y = paste('Can2 ', 
                 round(cda_results_combined$eigenvaluesAsPercentages[[2]] * 100, 
                       digits = 2), '%')) +
  theme_bw() +
  xlim(-max_range_combined, max_range_combined) +
  ylim(-max_range_combined, max_range_combined)

plot_cda_combined <- ggplot(cda_scores_combined, aes(x = Can1, y = Can2)) + 
  geom_point(aes(color = species, shape = species), 
             stroke = 1, 
             size = 3) +
  labs(x = paste('Can1 ', 
                 round(cda_results_combined$eigenvaluesAsPercentages[[1]] * 100, 
                       digits = 2), '%'), 
       y = paste('Can2 ', 
                 round(cda_results_combined$eigenvaluesAsPercentages[[2]] * 100, 
                       digits = 2), '%')) +
  theme_bw() +
  geom_polygon(data = hulls_combined, 
               aes(x = Can1, y = Can2, group = species, fill = species), 
               alpha = 0.2, 
               color = NA, 
               show.legend = FALSE) +
  scale_color_manual(name = "Species", 
                     values = c("L" = "#FF5980", 
                               "M" = "#28e21b", 
                               "R" = "#302790")) +
  scale_shape_manual(name = "Species", 
                    values = c("L" = 15, 
                              "M" = 16, 
                              "R" = 18))


# ============================================================================
# SECTION 7: Create Combined Plots and Save Outputs
# ============================================================================

# Save combined CDA plots (without hybrids)
pdf('morphometrics_avg_wo_hybrids.pdf', width = 15, height = 10.5)
print(ggarrange(plot_cda_fruits, 
                plot_cda_flowering, 
                plot_cda_short, 
                plot_cda_combined, 
                ncol = 2, 
                nrow = 2, 
                labels = c('A', 'B', 'C', 'D'), 
                common.legend = TRUE, 
                legend = 'right'))
dev.off()

# Save combined loadings plots
pdf('loadings_avg_wo_hybrids.pdf', width = 15, height = 10.5)
print(ggarrange(plot_loadings_fruits, 
                plot_loadings_flowering, 
                plot_loadings_short, 
                plot_loadings_combined, 
                ncol = 2, 
                nrow = 2, 
                labels = c('A', 'B', 'C', 'D'), 
                common.legend = TRUE, 
                legend = 'right'))
dev.off()

# Save loadings tables
write.csv(loadings_flowering, 
          'loadings_flowering_shoots_avg.csv', 
          row.names = FALSE)
write.csv(loadings_short, 
          'loadings_short_shoots_avg.csv', 
          row.names = FALSE)
write.csv(loadings_fruits, 
          'loadings_fruits_avg.csv', 
          row.names = FALSE)
write.csv(loadings_combined, 
          'loadings_completed_morphology.csv', 
          row.names = FALSE)


# ============================================================================
# SECTION 8: LDA Classification Accuracy Assessment
# ============================================================================

#' Calculate LDA classification accuracy from CDA coordinates
#' Uses leave-one-out cross-validation to assess how well CDA separates groups
#' @param df Data frame with Can1, Can2, and species columns
#' @param active_groups Vector of species codes to include in analysis
#' @param print_matrix Logical: print confusion matrix?
#' @return List with accuracy and predicted classes (invisibly)
lda_cda_accuracy <- function(df, active_groups = c("R", "M", "L"), 
                             print_matrix = TRUE) {
  # Filter to active groups and remove incomplete cases
  dat <- df[df$species %in% active_groups, ]
  dat <- dat[complete.cases(dat[, c("Can1", "Can2", "species")]), ]
  
  # Check for sufficient samples
  if (nrow(dat) < length(active_groups) * 2) {
    cat("Not enough samples in all groups for reliable LDA.\n")
    return(invisible(NULL))
  }
  
  # Run LDA with leave-one-out cross-validation
  lda_cv <- lda(species ~ Can1 + Can2, data = dat, CV = TRUE)
  accuracy <- mean(lda_cv$class == dat$species)
  
  # Print results
  cat("---- LDA cross-validated classification accuracy ----\n")
  cat(sprintf("Groups used: %s\n", paste(active_groups, collapse = ", ")))
  cat(sprintf("Accuracy: %.2f\n", accuracy))
  if (print_matrix) {
    cat("Confusion matrix:\n")
    print(table(True = dat$species, Pred = lda_cv$class))
  }
  invisible(list(accuracy = accuracy, pred = lda_cv$class))
}

# Calculate accuracy for each dataset
cat("\n=== Flowering Shoots ===\n")
lda_cda_accuracy(cda_scores_flowering, active_groups = c("R", "M", "L"))

cat("\n=== Fruits ===\n")
lda_cda_accuracy(cda_scores_fruits, active_groups = c("R", "M", "L"))

cat("\n=== Combined Morphology ===\n")
lda_cda_accuracy(cda_scores_combined, active_groups = c("R", "M", "L"))


# ============================================================================
# SECTION 9: Hybrid Analysis with Admixture Coloring
# ============================================================================

# Define species levels and plotting parameters
species_levels_all <- c("L", "R", "M", "MAC", "MED", "SUB")
shape_values_all <- c("L" = 15, "R" = 18, "M" = 16, 
                     "MAC" = 3, "MED" = 4, "SUB" = 8)
size_values_all <- c("L" = 4, "R" = 4, "M" = 4, 
                    "MAC" = 2, "MED" = 2, "SUB" = 2)

# Prepare combined morphology data with hybrids
cda_scores_with_hybrids <- cda_results_combined$obj$scores
cda_scores_with_hybrids$sample <- rownames(cda_scores_with_hybrids)
cda_scores_with_hybrids$Genetic_group <- reassign$group[
  match(cda_scores_with_hybrids$sample, reassign$Name)
]
cda_scores_with_hybrids$L <- reassign$L[
  match(cda_scores_with_hybrids$sample, reassign$Name)
]
cda_scores_with_hybrids$R <- reassign$R[
  match(cda_scores_with_hybrids$sample, reassign$Name)
]
cda_scores_with_hybrids$M <- reassign$M[
  match(cda_scores_with_hybrids$sample, reassign$Name)
]

# Create species column
cda_scores_with_hybrids <- cda_scores_with_hybrids %>%
  mutate(sp = case_when(
    Genetic_group %in% c("R1", "R2", "R3", "R") ~ "R",
    Genetic_group %in% c("L1", "L2", "L") ~ "L",
    Genetic_group %in% c("MAC") ~ "MAC",
    Genetic_group %in% c("MED") ~ "MED",
    Genetic_group %in% c("SUB") ~ "SUB",
    Genetic_group %in% c("M") ~ "M"
  ))

cda_scores_with_hybrids$L <- as.numeric(as.character(cda_scores_with_hybrids$L))
cda_scores_with_hybrids$R <- as.numeric(as.character(cda_scores_with_hybrids$R))
cda_scores_with_hybrids$M <- as.numeric(as.character(cda_scores_with_hybrids$M))
cda_scores_with_hybrids$species <- cda_scores_with_hybrids$sp

# --- 9.1 MAC hybrid plot (L x R) ---
cda_scores_mac <- cda_scores_with_hybrids %>% 
  filter(!Genetic_group %in% c('M', 'MED', 'SUB'))

missing_species_mac <- setdiff(species_levels_all, unique(cda_scores_mac$sp))
dummy_data_mac <- create_dummy_data(missing_species_mac)
cda_scores_mac <- bind_rows(cda_scores_mac, dummy_data_mac)

# Mix colors for MAC hybrids
color_L <- "#FF5980"
color_R <- "#302790"
cda_scores_mac$mixed_color <- mapply(mix_colors, 
                                     color_L, 
                                     color_R, 
                                     cda_scores_mac$L, 
                                     cda_scores_mac$R)

cda_scores_mac$sp <- factor(cda_scores_mac$sp, levels = species_levels_all)

plot_cda_mac <- ggplot(cda_scores_mac, aes(x = Can1, y = Can2)) +
  geom_point(aes(shape = sp, 
                color = mixed_color, 
                stroke = ifelse(sp == "MAC", 2, 1), 
                size = sp)) +
  labs(x = paste('Can1 ', 
                 round(cda_results_combined$eigenvaluesAsPercentages[[1]] * 100, 
                       digits = 2), '%'), 
       y = paste('Can2 ', 
                 round(cda_results_combined$eigenvaluesAsPercentages[[2]] * 100, 
                       digits = 2), '%')) +
  scale_color_identity() +
  scale_shape_manual(name = "Species", 
                     values = shape_values_all, 
                     drop = FALSE) +
  scale_size_manual(values = size_values_all, 
                   drop = FALSE) +
  theme_bw() +
  guides(size = "none",
         shape = guide_legend(override.aes = list(size = 4)))

# --- 9.2 MED hybrid plot (L x M) ---
cda_scores_med <- cda_scores_with_hybrids %>% 
  filter(!Genetic_group %in% c('R', 'R2','R1', 'R3', 'MAC', 'SUB'))

color_L <- "#FF5980"
color_M <- "#28e21b"
cda_scores_med$mixed_color <- mapply(mix_colors, 
                                     color_L, 
                                     color_M, 
                                     cda_scores_med$L, 
                                     cda_scores_med$M)

cda_scores_med$sp <- factor(cda_scores_med$sp, levels = species_levels_all)

plot_cda_med <- ggplot(cda_scores_med, aes(x = Can1, y = Can2)) +
  geom_point(aes(shape = sp, 
                color = mixed_color, 
                stroke = ifelse(sp == "MED", 2, 1), 
                size = sp)) +
  labs(x = paste('Can1 ', 
                 round(cda_results_combined$eigenvaluesAsPercentages[[1]] * 100, 
                       digits = 2), '%'), 
       y = paste('Can2 ', 
                 round(cda_results_combined$eigenvaluesAsPercentages[[2]] * 100, 
                       digits = 2), '%')) +
  scale_color_identity() +
  scale_shape_manual(name = "Species", 
                     values = shape_values_all, 
                     drop = FALSE) +
  scale_size_manual(values = size_values_all, 
                   drop = FALSE) +
  theme_bw() +
  guides(size = "none",
         shape = guide_legend(override.aes = list(size = 4)))

# --- 9.3 SUB hybrid plot (R x M) ---
cda_scores_sub <- cda_scores_with_hybrids %>% 
  filter(!Genetic_group %in% c('L', 'L2','L1', 'MAC', 'MED'))

color_R <- "#302790"
color_M <- "#28e21b"
cda_scores_sub$mixed_color <- mapply(mix_colors, 
                                    color_R, 
                                    color_M, 
                                    cda_scores_sub$R, 
                                    cda_scores_sub$M)

cda_scores_sub$sp <- factor(cda_scores_sub$sp, levels = species_levels_all)

plot_cda_sub <- ggplot(cda_scores_sub, aes(x = Can1, y = Can2)) +
  geom_point(aes(shape = sp, 
                color = mixed_color, 
                stroke = ifelse(sp == "SUB", 1.5, 1), 
                size = sp)) +
  labs(x = paste('Can1 ', 
                 round(cda_results_combined$eigenvaluesAsPercentages[[1]] * 100, 
                       digits = 2), '%'), 
       y = paste('Can2 ', 
                 round(cda_results_combined$eigenvaluesAsPercentages[[2]] * 100, 
                       digits = 2), '%')) +
  scale_color_identity() +
  scale_shape_manual(name = "Species", 
                     values = shape_values_all, 
                     drop = FALSE) +
  scale_size_manual(values = size_values_all, 
                   drop = FALSE) +
  theme_bw() +
  guides(size = "none",
         shape = guide_legend(override.aes = list(size = 4)))

# --- 9.4 Save hybrid plots ---
pdf('hybrids.pdf', width = 10.5, height = 15)
print(ggarrange(plot_cda_mac,
                plot_cda_med,
                plot_cda_sub,
                ncol = 1,
                nrow = 3,
                labels = c('A', 'B', 'C'),
                common.legend = TRUE,
                legend = 'right'))
dev.off()

# --- 9.5 LDA accuracy for hybrid groups ---
cat("\n=== Hybrid Classification Accuracy ===\n")
cat("\nMAC vs R, L:\n")
lda_cda_accuracy(cda_scores_with_hybrids, active_groups = c("MAC", "R", "L"))

cat("\nMED vs L, M:\n")
lda_cda_accuracy(cda_scores_with_hybrids, active_groups = c("MED", "L", "M"))

cat("\nSUB vs R, M:\n")
lda_cda_accuracy(cda_scores_with_hybrids, active_groups = c("SUB", "R", "M"))

cat("\nAll groups:\n")
lda_cda_accuracy(cda_scores_with_hybrids, 
                active_groups = c("L", "R", "M", "MED", "MAC", "SUB"))

