# Same-Slide Spatial Multi-Omics Integration Reveals Tumor Virus-Linked Spatial Reorganization of the Tumor Microenvironment

![overview](https://github.com/SizunJiangLab/IN-DEPTH/blob/main/docs/assets/images/overview.png)

## Abstract

The advent of spatial transcriptomics and spatial proteomics has enabled profound insights into tissue organization, however, these technologies remain largely independent, and emerging same slide multi-omics approaches are limited in plex, spatial resolution, and integrative analytics. We introduce <ins>***IN-situ DEtailed Phenotyping To High-resolution transcriptomics (IN-DEPTH)***</ins>, a streamlined, resource-efficient, and commercially platform-compatible workflow that uses single-cell spatial proteomics-derived imaging to guide transcriptomic capture on the same slide without RNA signal loss. To integrate modalities beyond domain- or niche-level mapping, we developed Spectral Graph Cross-Correlation (SGCC), an integrated proteomic-transcriptomic framework that resolves spatially coordinated functional state changes across interacting cell populations. Applied to diffuse large B-cell lymphoma (DLBCL), IN-DEPTH and SGCC enabled stepwise discovery from multicellular comparisons of Epstein–Barr virus (EBV)-positive and EBV-negative tumors to single-cell analyses revealing coordinated tumor–macrophage–CD4 T-cell remodeling. This approach identified enrichment of immunosuppressive C1Q macrophage polarization and CD4 T-cell dysfunction and further localized these programs to viral oncoprotein LMP1-driven tumor states and a candidate IL-27–STAT3 signaling axis. Collectively, IN-DEPTH enables scalable spatial multi-omics to uncover clinically relevant microenvironmental mechanisms.

## Tutorials

- [Experiment Protocols](https://sizunjianglab.github.io/IN-DEPTH/protocols/): Detailed Protocols for performing IN-DEPTH on various proteomics and transcriptomics platforms.
- [Data Integration Tutorials](https://github.com/SizunJiangLab/IN-DEPTH/blob/tutorial/tutorial/indepth_codex_geomx.ipynb): Detailed tutorials on integrating proteomics and transcriptomics data via image registration.
- [SGCC Tutorial](tutorial/sgcc_tutorial.md): R package documentation, installation instructions, and a full functions overview.

## Copyright & Licence

The IN-DEPTH codebase and data are released to the academic community for non-commercial academic research only. **Any commercial research use, integration into commercial products or services requires prior approvals.**

## Script for figures

An overview of the purpose of different scripts to reproduce the figures in the manuscript.

| File Name | Description |
| --- | --- |
| `paper_figures/01_figure1.R` | Code for plots in figure 1 |
| `paper_figures/02_figure2.Rmd` | Code for plots and analyses in figure 2 |
| `paper_figures/02_figure2E_supp2EFG.Rmd` | Code for plots and analyses in figure 2E and supplementary figure 2EFG |
| `paper_figures02_figure2F_cNMF.ipynb` | Code for running cNMF in figure 2F |
| `paper_figures02_figure2C_and_SuppFig2H.R` | Code for plots in figure 2C and supplementary figure 2H |
| `paper_figures03_figure3D_and_SuppFig3D.R` | Code for plots in figure 3D and supplementary figure 3D |
| `paper_figures03_figure3E_and_SuppFig3G.R` | Code for plots in figure 3E and supplementary figure 3G |
| `paper_figures03_SuppFig3A&B.R` | Code for supplementary figures 3A and 3B |
| `paper_figures03_SuppFig3C.R` | Code for supplementary figure 3C |
| `paper_figures04_figure4_4G&J.Rmd` | Code for plots in figure 4G and 4J |
| `paper_figures04_figure4_suppFig6.Rmd` | Code for plots in figure 4 and supplementary figure 6 |
| `paper_figures05_figure5A&B_SuppFig7A&B_half.R` | Code for plots in figure 5A, B and supplementary figures 7A, B (half portion) |
| `paper_figures05_figure5A&B&C&D_SuppFig7A&B.R` | Code for plots in figure 5A, B, C, D and supplementary figures 7A, B |
| `paper_figures05_figure5E.R` | Code for plots in figure 5E |

## Overview of code

An overview of code in `src`. Scripts in this folder contain the preprocessing steps used to generate input data for scripts in `paper_figures`.

| Folder name | Description |
| --- | --- |
| 01_figure_1_CODEXonly_vs_postCODEX | Code for the correlation analysis in figure 1. |
| 02_figure_2_CODEX_GeoMx_Tonsil_run | Pipeline for proteomics data preprocessing and analysis in figure 2. |
| 04_figure_4_CODEX_GeoMMx | Pipeline for proteomics data preprocessing and analysis in figure 4. |
| 05_figure_5_CODEX_GeoMx_analysis | Pipeline for transcriptomics data preprocessing and analysis in figure 5. |
| 05_figure_5_SGCC | Pipeline for SGCC analysis in figure 5. |
| 06_figure_6_scSGCC | Pipeline for scSGCC analysis in figure 6. |
| 07_figure_7_CODEX_pipeline | Pipeline for proteomics data preprocessing and analysis in figure 7. The `notebook/` folder contains the python notebook to perform the analysis for figure 7B. |
