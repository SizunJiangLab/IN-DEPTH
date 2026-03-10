library(tidyverse)
source("src/CODEXonly_vs_postCODEX/visualization.R")

# Figure 1C #####
## cosmx ####
output_figure_cosmx("cosmx_1k")
output_figure_cosmx("cosmx_6k")

## geomx ####
output_figure_geomx("geomx_DLBCL")
output_figure_geomx("geomx_LN")
output_figure_geomx("geomx_RCC")
output_figure_geomx("geomx_USC")
output_figure_geomx("geomx_cycling")
output_figure_geomx("geomx_fusion")

# visualization of correlation plot in the same x-y axis #####
output_figure_cor_slide_gene()

## VisiumHD ####
output_figure_visiumhd("visiumhd_only=jiang_post=postcodex")
output_figure_visiumhd("visiumhd_only=jiang_post=postcodex_cpb10000")

## correlation heatmap ####
output_heatmap_roi_pearson_cor()
