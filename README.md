# Morphometric and Genetic Analysis Scripts

This repository contains R scripts for analyzing morphometric and genetic data from a paper titeled "Wide diversity in narrow geographic space: genetic, morphological and ploidy variation in three Central European Crataegus species with emphasis on their reproductive modes" and accepted in Annals of Botany Plants in November 2025. The scripts perform geometric morphometrics, canonical discriminant analysis, and genetic principal coordinates analysis.

## Scripts Overview

### 1. `01_relative_warp_analysis.R`
**Purpose:** Performs Relative Warp (RW) analysis on geometric morphometric landmark data.

**What it does:**
- Reads TPS landmark files for flowering and short shoots
- Performs Generalized Procrustes Analysis (GPA)
- Averages coordinates per sample (multiple leaves per sample)
- Runs Relative Warp analysis using `Morpho::relWarps`
- Creates morphospace plots with density contours
- Generates thin-plate spline (TPS) grids showing shape differences
- Analyzes combined datasets (flowering + short shoots)
- Includes fruit PCA analysis

**Key outputs:**
- `relative_warp_morphospace_MorphoRW.pdf`: Morphospace plots
- `relative_warp_morphospace_MorphoRW_combined_paired.pdf`: Combined analysis
- `TPS_grid_*.pdf`: Shape deformation grids
- `fruits_PCA_cent.pdf`: Fruit PCA plot

---

### 2. `02_canonical_discriminant_analysis.R`
**Purpose:** Performs Canonical Discriminant Analysis (CDA) on morphometric data.

**What it does:**
- Processes flowering shoots, short shoots, fruits, and combined morphology
- Removes highly correlated variables (correlation > 0.9)
- Performs CDA with hybrids as passive samples
- Creates CDA plots with convex hulls
- Generates loadings plots showing variable contributions
- Calculates LDA classification accuracy
- Creates hybrid-specific plots with admixture-based coloring

**Key outputs:**
- `morphometrics_avg_wo_hybrids.pdf`: CDA plots for all datasets
- `loadings_avg_wo_hybrids.pdf`: CDA loadings plots
- `hybrids.pdf`: Hybrid-specific CDA plots
- `loadings_*.csv`: CDA loadings tables
- `*_avg.tsv`: Averaged datasets

---

### 3. `03_genetic_pcoa_analysis.R`
**Purpose:** Performs Principal Coordinates Analysis (PCoA) on SSR genetic data.

**What it does:**
- Reads GenoDive format SSR data
- Calculates Bruvo distance matrix (appropriate for variable ploidy)
- Performs PCoA on genetic distances
- Creates genetic morphospace plot with density contours
- Colors samples based on admixture proportions

**Key outputs:**
- `genetic_PCoA_densities_ggplot.pdf`: Genetic PCoA plot

---

### 4. `04_ploidy_analysis.R`
**Purpose:** Creates boxplots and frequency plots for ploidy data from seeds and leaves.

**What it does:**
- Reads ploidy data from seeds (embryo and endosperm) and leaves
- Calculates ploidy frequencies per species
- Creates frequency plots with point sizes proportional to frequency
- Generates index plots for embryo and endosperm (using Carex and Pisum standards)
- Combines all plots into a multi-panel figure

**Key outputs:**
- `fcm_fccs_plots.pdf`: Multi-panel figure with 8 plots showing ploidy frequencies and indices

**Author:** Soňa Píšová

---

### 5. `05_genetic_triplot.R`
**Purpose:** Creates a ternary plot showing admixture proportions of three Crataegus species.

**What it does:**
- Reads SSR admixture proportions (C. monogyna, C. laevigata, C. rhipidophylla)
- Transforms ternary coordinates to Cartesian coordinates
- Creates ternary plot with grid lines and percentage labels
- Colors and shapes points based on genetic groups
- Adds legend with group descriptions

**Key outputs:**
- `triplot.pdf`: Ternary plot showing admixture proportions

**Author:** Soňa Píšová

---

## Required Input Files

**Essential files:**
- `reassign.tsv`: Admixture and species assignment table (required by all scripts)
  - Must contain columns: `Name`, `SSRs.group_new`, `species`, `L`, `M`, `R`

**Script-specific files:**
- **Script 01:** `Hybrid_Zones_LM_fl.tps`, `Hybrid_Zones_LM_sh.tps`, `fruits_avg.tsv`
- **Script 02:** `Hybrid_Zones_LM_fl_procrustes.tsv`, `Hybrid_Zones_LM_sh_procrustes.tsv`, `all_measured_fruits.tsv`, `completed_morphology.tsv`
- **Script 03:** `Crataegus_ssr_data_sp_clones_removed_dosage_restored.gdv`
- **Script 04:** `Crataegus_seed_FCCS.xlsx`, `Crataegus_leaf_FCM.xlsx`
- **Script 05:** `Crataegus_ssr_triplot.csv`

---

## Installation

### Required R Packages

Install the following R packages before running the scripts:

```r
install.packages(c("tidyverse", "geomorph", "ggrepel", "abind", "here", 
                   "ggnewscale", "ggpubr", "FactoMineR", "scales", "MASS",
                   "readxl", "cowplot", "ggplot2"))
install.packages("BiocManager")
BiocManager::install("Morpho")
install.packages("devtools")
devtools::install_github("MarekSlenker/MorphoTools2")
install.packages("polysat")
install.packages("ape")
```

### Package Versions

The scripts were tested with:
- R version 4.5.1
- tidyverse 2.0.0
- geomorph (latest)
- MorphoTools2 1.0.2.1

---

## Usage

### Running Individual Scripts

Each script can be run independently:

```bash
# From command line
Rscript 01_relative_warp_analysis.R
Rscript 02_canonical_discriminant_analysis.R
Rscript 03_genetic_pcoa_analysis.R
Rscript 04_ploidy_analysis.R
Rscript 05_genetic_triplot.R
```

Or from R/RStudio:

```r
source("01_relative_warp_analysis.R")
source("02_canonical_discriminant_analysis.R")
source("03_genetic_pcoa_analysis.R")
source("04_ploidy_analysis.R")
source("05_genetic_triplot.R")
```

### Running All Scripts

To run all analyses in sequence:

```bash
Rscript 01_relative_warp_analysis.R
Rscript 02_canonical_discriminant_analysis.R
Rscript 03_genetic_pcoa_analysis.R
Rscript 04_ploidy_analysis.R
Rscript 05_genetic_triplot.R
```

---

## Script Portability

All scripts are designed to be portable and will:
- Automatically detect their own location
- Work in RStudio and when run via `Rscript`
- Use relative paths (all files in the same directory)
- Work on any operating system (Windows, macOS, Linux)

**Important:** Place all input files in the same directory as the scripts.

---

## Output Files

### PDF Files
All plots are saved as PDF files (as per user preference). Main outputs include:
- Morphospace plots (Relative Warp analysis)
- CDA plots (Canonical Discriminant Analysis)
- Loadings plots (CDA variable contributions)
- TPS deformation grids (shape differences)
- Genetic PCoA plots (Principal Coordinates Analysis)
- Ploidy frequency and index plots 
- Genetic triplot 

### CSV Files
- CDA loadings tables for each dataset

### TSV Files
- Averaged datasets (intermediate files, also used as inputs)

---

## Notes

1. **Ploidy and Repeat Units (Script 03):** The ploidy and repeat unit vectors are hardcoded and must match your GenoDive file structure. 

2. **Excluded Samples:** Several samples are excluded from analyses (e.g., "L101", "HMAC611", etc.). These are hardcoded in the scripts.

3. **Color Scheme:**
   - L species: #FF5980 (pink/red)
   - M species: #28e21b (green)
   - R species: #302790 (purple)
   - Hybrids: Mixed colors based on admixture proportions

4. **Shape Symbols:**
   - L: Square (15)
   - M: Circle (16)
   - R: Diamond (18)
   - MAC: Plus (3)
   - MED: Cross (4)
   - SUB: Asterisk (8)

---

## Citation

If you use these scripts in your research, please cite:
- The relevant R packages (see script headers)
- This repository

---

## Authors

**Scripts 01-03:** Roman Ufimov  
**Scripts 04-05:** Soňa Píšová  

July 2025

---

## License

These scripts are provided as-is for reproducibility purposes. Please cite appropriately if used in research.
