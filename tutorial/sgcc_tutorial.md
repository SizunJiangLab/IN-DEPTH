# SGCC

## R Package Documentation

This document provides a comprehensive guide to the functions included in the R script, designed for managing package installation, performing graph-based calculations, and analyzing spatial and multi-omics data.

## Install and Load Packages

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

## Functions Overview

### 1. `install_and_load`

- **Purpose:** Installs and loads the necessary R packages from CRAN or GitHub repositories.
- **Parameters:**
  - `packages`: A named vector with package names as keys and source URLs as values.

### 2. `cal_laplacian`

- **Purpose:** Calculates the Laplacian matrix of a given adjacency matrix.
- **Parameters:**
  - `W`: Adjacency matrix.
  - `type`: Type of Laplacian to compute (`"unnormalized"`, `"normalized"`, or `"randomwalk"`).

### 3. `FastDecompositionLap`

- **Purpose:** Performs a fast eigen decomposition of the Laplacian matrix.
- **Parameters:**
  - `laplacianMat`: Laplacian matrix to decompose.
  - `k_fold`: Factor to determine the number of eigenvalues/vectors to compute.
  - `which`: Part of the spectrum to target.
  - `sigma`: Optional shift-invert parameter.
  - `opts`: List of additional options for `eigs_sym`.

### 4. `gft`

- **Purpose:** Computes the Graph Fourier Transform of a given signal using the eigenvectors of the Laplacian matrix.
- **Parameters:**
  - `signal`: Input signal vector or matrix.
  - `U`: Matrix of eigenvectors.

### 5. `plot_FM`

- **Purpose:** Plots frequency modulation from Graph Fourier Transform data.
- **Parameters:**
  - `input`: Data containing spatial coordinates and frequencies.
  - `FM_idx`: Indices of frequency columns to plot.
  - `ncol`: Number of columns in the plot layout.

### 6. `Cal_Eigen`

- **Purpose:** Calculates the eigenvalues and eigenvectors of the graph Laplacian, and determines the cutoff for low-frequency components.
- **Parameters:**
  - `data.in`: Input data with spatial coordinates.
  - `k`: Number of nearest neighbors.
  - `k_fold`: Factor to determine the number of eigenvalues/vectors to compute.
  - `sensitivity`: Sensitivity parameter for knee point detection.

### 7. `Cal_GCC`

- **Purpose:** Calculates the Spectral Graph Cross-Correlation between two signals using a subset of low-frequency components.
- **Parameters:**
  - `data.in`: Input data frame.
  - `knee`: Knee point for low-frequency cutoff.
  - `signal1`: Column name of the first signal in `data.in`.
  - `signal2`: Column name of the second signal in `data.in`.
  - `eigenvector`: Matrix of eigenvectors.
