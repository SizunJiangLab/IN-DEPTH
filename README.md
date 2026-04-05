# Same-Slide Spatial Multi-Omics Integration with IN-DEPTH Reveals Tumor Virus-Linked Spatial Reorganization of the Tumor Microenvironment

![overview](https://github.com/SizunJiangLab/IN-DEPTH/blob/main/docs/assets/images/overview.png)

> **Published in *Cancer Discovery* (2026):** [https://doi.org/10.1158/2159-8290.CD-25-0775](https://doi.org/10.1158/2159-8290.CD-25-0775)

## Updates

- **[2026-03]** Published in *Cancer Discovery*! [[Link]](https://doi.org/10.1158/2159-8290.CD-25-0775)
- **[2024-12]** Preprint available on bioRxiv. [[Link]](https://doi.org/10.1101/2024.12.20.629650)

## Abstract

Spatial transcriptomics and proteomics have enabled profound insights into tissue organization, yet these technologies remain largely disparate, and emerging same-slide multi-omics approaches are limited in plex, spatial resolution, signal retention, and integrative analytics. We introduce <ins>***IN-situ DEtailed Phenotyping To High-resolution transcriptomics (IN-DEPTH)***</ins>, a streamlined, resource-efficient, commercially compatible workflow using single-cell spatial proteomics-derived imaging to guide transcriptomic capture on the same slide without RNA signal loss. To integrate modalities beyond niche-level mapping, we developed Spectral Graph Cross-Correlation (SGCC), a proteomic-transcriptomic framework resolving spatially coordinated functional state changes across interacting cell populations. Applied to diffuse large B-cell lymphoma (DLBCL), IN-DEPTH and SGCC enabled stepwise discovery from EBV-positive and EBV-negative tumor comparisons to single-cell resolution, revealing coordinated tumor–macrophage–CD4 T-cell remodeling, immunosuppressive C1Q macrophage enrichment, CD4 T-cell dysfunction, and a candidate IL-27–STAT3 signaling axis. Collectively, IN-DEPTH enables scalable spatial multi-omics to uncover clinically relevant microenvironmental mechanisms and towards robust spatial multi-modal AI models.

## Tutorials

- [Experiment Protocols](https://sizunjianglab.github.io/IN-DEPTH/protocols/): Detailed Protocols for performing IN-DEPTH on various proteomics and transcriptomics platforms.
- [Data Integration Tutorials](https://github.com/SizunJiangLab/IN-DEPTH/blob/tutorial/tutorial/indepth_codex_geomx.ipynb): Detailed tutorials on integrating proteomics and transcriptomics data via image registration.
- [SGCC Tutorial](tutorial/sgcc_tutorial.md): R package documentation, installation instructions, and a full functions overview.

## Reproducibility

### Figure Scripts

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
| [`06_figure_6B_supp_9E.R`](paper_figures/06_figure_6B_supp_9E.R) | Code for plots in figure 6B and supplementary figure 9E |
| [`07_figure_7CD.R`](paper_figures/07_figure_7CD.R) | Code for plots in figure 7C and 7D |
| [`07_figure_7E.R`](paper_figures/07_figure_7E.R) | Code for plots in figure 7E |
| [`07_figure_7F_supp_10A.R`](paper_figures/07_figure_7F_supp_10A.R) | Code for plots in figure 7F and supplementary figure 10A |
| [`07_figure_7G.py`](paper_figures/07_figure_7G.py) | Code for plots in figure 7G |

</details>

### Analysis Pipelines

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
| [`06_figure_6_correlation`](src/06_figure_6_correlation) | Pipeline for AUCell scoring and Spearman correlation analysis in figure 6 |
| [`07_figure_7_CODEX_pipeline`](src/07_figure_7_CODEX_pipeline) | Pipeline for proteomics data preprocessing and analysis in figure 7. The `notebook/` folder contains the python notebook to perform the analysis for figure 7B |
| [`07_figure_7_neighborhood_analysis`](src/07_figure_7_neighborhood_analysis) | Pipeline for spatial neighborhood analysis, DEG, GSEA, and ligand-receptor analysis in figure 7 |

</details>

### Data Availability

The processed data required to reproduce the analyses and figures are deposited on Zenodo under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/) license.

- Part 1 ([10.5281/zenodo.14530077](https://doi.org/10.5281/zenodo.14530077))
- Part 2 ([10.5281/zenodo.18379155](https://doi.org/10.5281/zenodo.18379155))

## Reference

```bibtex
@article{yiu2026indepth,
    title={Same-Slide Spatial Multi-Omics Integration with IN-DEPTH Reveals Tumor Virus-Linked Spatial Reorganization of the Tumor Microenvironment},
    author={Yiu, Stephanie Pei Tung and Chang, Yuzhou and Yeo, Yao Yu and Qiu, Huaying and Wu, Wenrui and Michel, Hendrik A and Jin, Xiaojie and Huang, Rongting and Kure, Shoko and Parmelee, Lindsay and Luo, Shuli and Cramer, Precious and Lee, Jia Le and Wang, Yang and Zhao, Zhangxin and Yeung, Jason and El Ahmar, Nourhan and Simsek, Berkay and Mohanna, Razan and Van Orden, McKayla and Lu, Wesley S and Livak, Kenneth J and Li, Shuqiang and Gao, Ce and Burgess, Melinda and Keane, Colm and Shahryari, Jahanbanoo and Kingsley, Leandra G and Al-Humadi, Reem N and Nasr, Sahar and Nkosi, Dingani and Sadigh, Sam and Rock, Philip and Frauenfeld, Leonie and Kaufmann, Louisa and Zhu, Bokai and Basak, Ankit and Dhanikonda, Nagendra and Chan, Chi Ngai and Krull, Jordan and Cho, Ye Won and Chen, Chia-Yu and Brown, Jonathan and Wang, Hongbo and Zhao, Bo and Lee, Jia-Ying Joey and Loo, Lit-Hsin and Kim, David M and Boussiotis, Vassiliki A and Zhang, Baochun and Wei, Kevin and Shalek, Alex K and Howitt, Brooke E and Signoretti, Sabina and Sch\"urch, Christian M and Hodi, F Stephen and Burack, W Richard and Rodig, Scott J and Ma, Qin and Jiang, Sizun},
    journal={Cancer Discovery},
    year={2026},
    doi={10.1158/2159-8290.CD-25-0775},
    pmid={41874448},
    publisher={American Association for Cancer Research}
}
```

## Copyright & License

The accompanying [publication](https://doi.org/10.1158/2159-8290.CD-25-0775) is distributed under the [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)](https://creativecommons.org/licenses/by-nc-nd/4.0/) license.

The software code in this repository is provided for **non-commercial academic research purposes only**. Any commercial use — including but not limited to incorporation into commercial products, services, or clinical diagnostic tools — requires prior written approval from the corresponding authors. See the [LICENSE](LICENSE) file for full terms.

The processed data deposited on Zenodo are shared under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) license.

## Contact

If you have questions, comments, or concerns, feel free to email **Sizun Jiang** (`sjiang3@bidmc.harvard.edu`) or **Stephanie Pei Tung Yiu** (`syiu@bidmc.harvard.edu`). All requests and questions will be answered in a timely manner. Immediate responses may not be available.
