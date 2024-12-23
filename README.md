# Same-Slide Spatial Multi-Omics Integration Reveals Tumor Virus-Linked Spatial Reorganization of the Tumor Microenvironment

![overview](https://github.com/SizunJiangLab/IN-DEPTH/blob/main/doc/overview.png)

## Abstract

The advent of spatial transcriptomics and spatial proteomics have enabled profound insights into tissue organization and provided systems-level understanding of diseases. However, both technologies remain largely independent, and emerging same slide spatial multi-omics approaches are generally limited in plex, spatial resolution, and analysis approaches. We introduce <ins>***IN-situ DEtailed Phenotyping To High-resolution transcriptomics (IN-DEPTH)***</ins>, a streamlined and resource-effective approach compatible with various spatial platforms. This iterative approach first entails single-cell spatial proteomics and rapid analysis to guide subsequent spatial transcriptomics capture on the same slide without loss in RNA signal. We also introduce k-bandlimited Spectral Graph Cross-Correlation (SGCC) to facilitate additional insights through integrative spatial multi-omics analysis. Application of IN-DEPTH and SGCC on lymphoid tissues demonstrated precise single-cell phenotyping and cell-type specific transcriptome capture, and accurately resolved the local and global transcriptome changes associated with the cellular organization of germinal centers. We then implemented IN-DEPTH and SGCC to dissect the tumor microenvironment (TME) of Epstein-Barr Virus (EBV)-positive and EBV-negative diffuse large B-cell lymphoma (DLBCL). Our results identified a key tumor-macrophage-CD4 T-cell immunomodulatory axis differently regulated between EBV-positive and EBV-negative DLBCL, and its central role in coordinating immune dysfunction and suppression. IN-DEPTH enables scalable, resource-efficient, and comprehensive spatial multi-omics dissection of tissues to advance clinically relevant discoveries.

## Overview of Code

An overview of the purpose of different scripts in this repository.

### Script for figures


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

## SGCC (The codes are in src/SGCC_code for R users)
### R Package Documentation

This document provides a comprehensive guide to the functions included in the R script, designed for managing package installation, performing graph-based calculations, and analyzing spatial and multi-omics data.

#### Functions Overview

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
- **Purpose:** Calculates the Graph Cross-Correlation between two signals using a subset of low-frequency components.
- **Parameters:**
  - `data.in`: Input data frame.
  - `knee`: Knee point for low-frequency cutoff.
  - `signal1`: Column name of the first signal in `data.in`.
  - `signal2`: Column name of the second signal in `data.in`.
  - `eigenvector`: Matrix of eigenvectors.



