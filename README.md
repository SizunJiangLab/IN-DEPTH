# Same-Slide Spatial Multi-Omics Integration Reveals Tumor Virus-Linked Spatial Reorganization of the Tumor Microenvironment

![overview](https://github.com/SizunJiangLab/IN-DEPTH/blob/main/docs/assets/images/overview.png)

## Abstract

The advent of spatial transcriptomics and spatial proteomics has enabled profound insights into tissue organization, however, these technologies remain largely independent, and emerging same slide multi-omics approaches are limited in plex, spatial resolution, and integrative analytics. We introduce <ins>***IN-situ DEtailed Phenotyping To High-resolution transcriptomics (IN-DEPTH)***</ins>, a streamlined, resource-efficient, and commercially platform-compatible workflow that uses single-cell spatial proteomics-derived imaging to guide transcriptomic capture on the same slide without RNA signal loss. To integrate modalities beyond domain- or niche-level mapping, we developed Spectral Graph Cross-Correlation (SGCC), an integrated proteomic-transcriptomic framework that resolves spatially coordinated functional state changes across interacting cell populations. Applied to diffuse large B-cell lymphoma (DLBCL), IN-DEPTH and SGCC enabled stepwise discovery from multicellular comparisons of Epstein–Barr virus (EBV)-positive and EBV-negative tumors to single-cell analyses revealing coordinated tumor–macrophage–CD4 T-cell remodeling. This approach identified enrichment of immunosuppressive C1Q macrophage polarization and CD4 T-cell dysfunction and further localized these programs to viral oncoprotein LMP1-driven tumor states and a candidate IL-27–STAT3 signaling axis. Collectively, IN-DEPTH enables scalable spatial multi-omics to uncover clinically relevant microenvironmental mechanisms.

## Tutorials

- [Experiment Protocols](https://sizunjianglab.github.io/IN-DEPTH/protocols/): Detailed Protocols for performing IN-DEPTH on various proteomics and transcriptomics platforms.
- [Data Integration Tutorials](https://github.com/SizunJiangLab/IN-DEPTH/blob/tutorial/tutorial/indepth_codex_geomx.ipynb): Detailed tutorials on integrating proteomics and transcriptomics data via image registration.
- [SGCC Tutorial](tutorial/sgcc_tutorial.md): R package documentation, installation instructions, and a full functions overview.

## Copyright & License

The IN-DEPTH codebase and data are released to the academic community for non-commercial academic research only. **Any commercial research use, integration into commercial products or services requires prior approvals.**

## Script for figures

An overview of scripts in `paper_figures/`. These scripts reproduce the figures in the manuscript.

<details>
<summary>Click to expand</summary>

| File Name | Description |
| --- | --- |
| [`01_figure_1.R`](paper_figures/01_figure_1.R) | Code for plots in figure 1 |
| [`02_figure_2.Rmd`](paper_figures/02_figure_2.Rmd) | Code for plots and analyses in figure 2 |
| [`02_figure_2C_supp_2D.R`](paper_figures/02_figure_2C_supp_2D.R) | Code for plots in figure 2C and supplementary figure 2D |
| [`02_figure_2D_supp_2H.R`](paper_figures/02_figure_2D_supp_2H.R) | Code for plots in figure 2D and supplementary figure 2H |
| [`02_figure_2E_supp_2EFG.Rmd`](paper_figures/02_figure_2E_supp_2EFG.Rmd) | Code for plots and analyses in figure 2E and supplementary figure 2EFG |
| [`02_figure_2F.ipynb`](paper_figures/02_figure_2F.ipynb) | Code for running cNMF in figure 2F |
| [`03_figure_3C.R`](paper_figures/03_figure_3C.R) | Code for plots in figure 3C |
| [`03_figure_3D_supp_3D.R`](paper_figures/03_figure_3D_supp_3D.R) | Code for plots in figure 3D and supplementary figure 3D |
| [`03_figure_3E_supp_3G.R`](paper_figures/03_figure_3E_supp_3G.R) | Code for plots in figure 3E and supplementary figure 3G |
| [`03_supp_3.R`](paper_figures/03_supp_3.R) | Code for supplementary figure 3 (revision) |
| [`03_supp_3AB.R`](paper_figures/03_supp_3AB.R) | Code for supplementary figures 3A and 3B |
| [`03_supp_3C.R`](paper_figures/03_supp_3C.R) | Code for supplementary figure 3C |
| [`04_figure_4_GJ.Rmd`](paper_figures/04_figure_4_GJ.Rmd) | Code for plots in figure 4G and 4J |
| [`04_figure_4_supp_6.Rmd`](paper_figures/04_figure_4_supp_6.Rmd) | Code for plots in figure 4 and supplementary figure 6 |
| [`05_figure_5AB_supp_7AB_half.R`](paper_figures/05_figure_5AB_supp_7AB_half.R) | Code for plots in figure 5A, B and supplementary figures 7A, B (half portion) |
| [`05_figure_5ABCD_supp_7AB.R`](paper_figures/05_figure_5ABCD_supp_7AB.R) | Code for plots in figure 5A, B, C, D and supplementary figures 7A, B |
| [`05_figure_5E.R`](paper_figures/05_figure_5E.R) | Code for plots in figure 5E |

</details>

## Overview of code

An overview of code in `src`. Scripts in this folder contain the preprocessing steps used to generate input data for scripts in `paper_figures`.

<details>
<summary>Click to expand</summary>

| Folder name | Description |
| --- | --- |
| [`01_figure_1_CODEXonly_vs_postCODEX`](src/01_figure_1_CODEXonly_vs_postCODEX) | Code for the correlation analysis in figure 1 |
| [`02_figure_2_CODEX_GeoMx_Tonsil_run`](src/02_figure_2_CODEX_GeoMx_Tonsil_run) | Pipeline for proteomics data preprocessing and analysis in figure 2 |
| [`04_figure_4_CODEX_GeoMMx`](src/04_figure_4_CODEX_GeoMx) | Pipeline for proteomics data preprocessing and analysis in figure 4 |
| [`05_figure_5_CODEX_GeoMx_analysis`](src/05_figure_5_CODEX_GeoMx_analysis) | Pipeline for transcriptomics data preprocessing and analysis in figure 5 |
| [`05_figure_5_SGCC`](src/05_figure_5_SGCC) | Pipeline for SGCC analysis in figure 5 |
| [`06_figure_6_scSGCC`](src/06_figure_6_scSGCC) | Pipeline for scSGCC analysis in figure 6 |
| [`07_figure_7_CODEX_pipeline`](src/07_figure_7_CODEX_pipeline) | Pipeline for proteomics data preprocessing and analysis in figure 7. The `notebook/` folder contains the python notebook to perform the analysis for figure 7B |

</details>
