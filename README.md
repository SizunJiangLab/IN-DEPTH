# Same-Slide Spatial Multi-Omics Integration Reveals Tumor Virus-Linked Spatial Reorganization of the Tumor Microenvironment

![overview](https://github.com/SizunJiangLab/IN-DEPTH/blob/main/docs/assets/images/overview.png)



## Abstract

The advent of spatial transcriptomics and spatial proteomics has enabled profound insights into tissue organization, however, these technologies remain largely independent, and emerging same slide multi-omics approaches are limited in plex, spatial resolution, and integrative analytics. We introduce <ins>***IN-situ DEtailed Phenotyping To High-resolution transcriptomics (IN-DEPTH)***</ins>, a streamlined, resource-efficient, and commercially platform-compatible workflow that uses single-cell spatial proteomics-derived imaging to guide transcriptomic capture on the same slide without RNA signal loss. To integrate modalities beyond domain- or niche-level mapping, we developed Spectral Graph Cross-Correlation (SGCC), an integrated proteomic-transcriptomic framework that resolves spatially coordinated functional state changes across interacting cell populations. Applied to diffuse large B-cell lymphoma (DLBCL), IN-DEPTH and SGCC enabled stepwise discovery from multicellular comparisons of Epstein–Barr virus (EBV)-positive and EBV-negative tumors to single-cell analyses revealing coordinated tumor–macrophage–CD4 T-cell remodeling. This approach identified enrichment of immunosuppressive C1Q macrophage polarization and CD4 T-cell dysfunction and further localized these programs to viral oncoprotein LMP1-driven tumor states and a candidate IL-27–STAT3 signaling axis. Collectively, IN-DEPTH enables scalable spatial multi-omics to uncover clinically relevant microenvironmental mechanisms.




## Tutorials

- [Experiment Protocols](https://sizunjianglab.github.io/IN-DEPTH/protocols/): Detailed Protocols for performing IN-DEPTH on various proteomics and transcriptomics platforms.
- [Data Integration Tutorials](https://github.com/SizunJiangLab/IN-DEPTH/blob/tutorial/tutorial/indepth_codex_geomx.ipynb): Detailed tutorials on integrating proteomics and transcriptomics data via image registration. 
- [SGCC Tutorial](): Detailed tutorials on SGCC.



## Copyright & Licence

The IN-DEPTH codebase and data are released to the academic community for non-commercial academic research only. **Any commercial research use, integration into commercial products or services requires prior approvals.**




## Script for figures

An overview of the purpose of different scripts in `paper_figures`.


| File Name                              | Description                                             |
|----------------------------------------|---------------------------------------------------------|
| 01_figure1.R                           | Code for plots in figure 1                              |
| 02_figure2.Rmd                         | Code for plots and analyses in figure 2                 |
| 02_figure2E_supp2EFG.Rmd               | Code for plots and analyses in figure 2E and supplementary figure 2EFG |
| 02_figure2F_cNMF.ipynb                 | Code for running cNMF in figure 2F                      |
| 02_figure2C_and_SuppFig2H.R            | Code for plots in figure 2C and supplementary figure 2H |
| 03_figure3D_and_SuppFig3D.R            | Code for plots in figure 3D and supplementary figure 3D |
| 03_figure3E_and_SuppFig3G.R            | Code for plots in figure 3E and supplementary figure 3G |
| 03_SuppFig3A&B.R                       | Code for supplementary figures 3A and 3B                |
| 03_SuppFig3C.R                         | Code for supplementary figure 3C                        |
| 04_figure4_4G&J.Rmd                    | Code for plots in figure 4G and 4J                      |
| 04_figure4_suppFig6.Rmd                | Code for plots in figure 4 and supplementary figure 6   |
| 05_figure5A&B_SuppFig7A&B_half.R       | Code for plots in figure 5A&B and supplementary figures 7A&B (half portion) |
| 05_figure5A&B&C&D_SuppFig7A&B.R        | Code for plots in figure 5A&B&C&D and supplementary figures 7A&B |
| 05_figure5E.R                          | Code for plots in figure 5E                             |



## Overview of code

An overview of code in `src`. Scripts in this folder contain the preprocessing steps used to generate input data for scripts in `paper_figures`.

| Folder name                            | Description                                             |
|----------------------------------------|---------------------------------------------------------|
| 01_figure_1_CODEXonly_vs_postCODEX         | Code for the correlation analysis in figure 1.          |
| 02_figure_2_CODEX_GeoMx_Tonsil_run                            | Pipeline for proteomics data preprocessing and analysis in figure 2.|
| 04_figure_4_CODEX_GeoMMx                | Pipeline for proteomics data preprocessing and analysis in figure 4. |
| 05_figure_5_CODEX_GeoMx_analysis                | Pipeline for transcriptomics data preprocessing and analysis in figure 5. |
| 05_figure_5_SGCC           | Pipeline for SGCC analysis in figure 5. |
| 06_figure_6_CODEX_pipeline           | Pipeline for proteomics data preprocessing and analysis in figure 6. |
| 06_figure_6_scSGCC         | Pipeline for scSGCC analysis in figure 6. |



## SGCC

### R Package Documentation

This document provides a comprehensive guide to the functions included in the R script, designed for managing package installation, performing graph-based calculations, and analyzing spatial and multi-omics data.

### Install and Load Packages
Run the following code in your R console. This code checks if each package is installed, installs any missing packages, and loads them into the R session.

```R
# Function to install and load packages
install_and_load <- function(packages) {
  # Ensure that the 'remotes' package is installed and loaded
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  library(remotes)
  
  # Loop through each package and install from the appropriate source
  for (pkg in names(packages)) {
    package_name <- pkg
    source_url <- packages[pkg]
    
    # Check if the package is already installed
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
      # Install from GitHub if URL is provided
      if (grepl("github", source_url)) {
        cat(sprintf("Installing %s from GitHub repository %s\n", package_name, source_url))
        remotes::install_github(source_url)
      } else {  # Install from CRAN
        cat(sprintf("Installing package: %s from CRAN\n", package_name))
        install.packages(package_name)
      }
      
      # Load the package
      library(package_name, character.only = TRUE)
    } else {
      cat(sprintf("Package already installed: %s\n", package_name))
    }
  }
}

# List of packages with their installation sources
packages <- c(
  spdep = "CRAN", 
  spatialEco = "CRAN", 
  sp = "CRAN", 
  RANN = "CRAN", 
  igraph = "CRAN", 
  RSpectra = "CRAN", 
  Matrix = "CRAN", 
  patchwork = "CRAN", 
  ggpubr = "CRAN", 
  egg = "CRAN",
  kneedle = "etam4260/kneedle", 
  gasper = "fabnavarro/gasper"
)

# Run the installation and loading function
install_and_load(packages)
```

### Functions Overview

#### 1. `install_and_load`
- **Purpose:** Installs and loads the necessary R packages from CRAN or GitHub repositories.
- **Parameters:**
  - `packages`: A named vector with package names as keys and source URLs as values.

#### 2. `cal_laplacian`
- **Purpose:** Calculates the Laplacian matrix of a given adjacency matrix.
- **Parameters:**
  - `W`: Adjacency matrix.
  - `type`: Type of Laplacian to compute (`"unnormalized"`, `"normalized"`, or `"randomwalk"`).

#### 3. `FastDecompositionLap`
- **Purpose:** Performs a fast eigen decomposition of the Laplacian matrix.
- **Parameters:**
  - `laplacianMat`: Laplacian matrix to decompose.
  - `k_fold`: Factor to determine the number of eigenvalues/vectors to compute.
  - `which`: Part of the spectrum to target.
  - `sigma`: Optional shift-invert parameter.
  - `opts`: List of additional options for `eigs_sym`.

#### 4. `gft`
- **Purpose:** Computes the Graph Fourier Transform of a given signal using the eigenvectors of the Laplacian matrix.
- **Parameters:**
  - `signal`: Input signal vector or matrix.
  - `U`: Matrix of eigenvectors.

#### 5. `plot_FM`
- **Purpose:** Plots frequency modulation from Graph Fourier Transform data.
- **Parameters:**
  - `input`: Data containing spatial coordinates and frequencies.
  - `FM_idx`: Indices of frequency columns to plot.
  - `ncol`: Number of columns in the plot layout.

#### 6. `Cal_Eigen`
- **Purpose:** Calculates the eigenvalues and eigenvectors of the graph Laplacian, and determines the cutoff for low-frequency components.
- **Parameters:**
  - `data.in`: Input data with spatial coordinates.
  - `k`: Number of nearest neighbors.
  - `k_fold`: Factor to determine the number of eigenvalues/vectors to compute.
  - `sensitivity`: Sensitivity parameter for knee point detection.

#### 7. `Cal_GCC`
- **Purpose:** Calculates the Spectral Graph Cross-Correlation between two signals using a subset of low-frequency components.
- **Parameters:**
  - `data.in`: Input data frame.
  - `knee`: Knee point for low-frequency cutoff.
  - `signal1`: Column name of the first signal in `data.in`.
  - `signal2`: Column name of the second signal in `data.in`.
  - `eigenvector`: Matrix of eigenvectors.



