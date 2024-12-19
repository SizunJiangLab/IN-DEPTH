# Same-Slide Spatial Multi-Omics Integration Reveals Tumor Virus-Linked Spatial Reorganization of the Tumor Microenvironment

![overview](https://github.com/SizunJiangLab/IN-DEPTH/blob/main/doc/overview.png)

## Abstract

The advent of spatial transcriptomics and spatial proteomics have enabled profound insights into tissue organization and provided systems-level understanding of diseases. However, both technologies remain largely independent, and emerging same slide spatial multi-omics approaches are generally limited in plex, spatial resolution, and analysis approaches. We introduce <ins>***IN-situ DEtailed Phenotyping To High-resolution transcriptomics (IN-DEPTH)***</ins>, a streamlined and resource-effective approach compatible with various spatial platforms. This iterative approach first entails single-cell spatial proteomics and rapid analysis to guide subsequent spatial transcriptomics capture on the same slide without loss in RNA signal. We also introduce k-bandlimited Spectral Graph Cross-Correlation (SGCC) to facilitate additional insights through integrative spatial multi-omics analysis. Application of IN-DEPTH and SGCC on lymphoid tissues demonstrated precise single-cell phenotyping and cell-type specific transcriptome capture, and accurately resolved the local and global transcriptome changes associated with the cellular organization of germinal centers. We then implemented IN-DEPTH and SGCC to dissect the tumor microenvironment (TME) of Epstein-Barr Virus (EBV)-positive and EBV-negative diffuse large B-cell lymphoma (DLBCL). Our results identified a key tumor-macrophage-CD4 T-cell immunomodulatory axis differently regulated between EBV-positive and EBV-negative DLBCL, and its central role in coordinating immune dysfunction and suppression. IN-DEPTH enables scalable, resource-efficient, and comprehensive spatial multi-omics dissection of tissues to advance clinically relevant discoveries.

## Overview of Code

An overview of the purpose of different scripts in this repository.

### Script for figures

| File Name     | Description |
| ------------- | ------------- |
| 01_figure1.R  | Code for plots in figure 1  |
| 02_figure2.Rmd | Code for plots and analyses in figure 2 |
| 02_figure2E_supp2EFG.Rmd | Code for plots and analyses in figure 2E and supplementary figure 2EFG |
| 02_figure2F_cNMF.ipynb | Code for running cNMF in figure 2F |
| 04_figure4_4G&J.Rmd | Code for plots in figure 4G and 4J |
| 04_figure4_suppFig6.Rmd | Code for plots in figure 4 and supplimentary figure 6 |
