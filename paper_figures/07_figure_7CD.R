# %% Load Libraries ==========
library(tidyverse)
library(ggrepel)
library(patchwork)

# %% Define plot_volcano ==========
#' Create volcano plot from differential expression results
#'
#' @param df_deg Data frame containing differential expression results
#' @param genes_to_label Character vector of gene names to label
#' @param title Plot title
#' @param col_p Column name for p-values (default: "pvals_adj")
#' @param col_log2fc Column name for log2 fold changes (default: "logfoldchanges")
#' @param col_names Column name for gene names (default: "names")
#' @param cutoff_log2fc Log2 fold change cutoff for significance (default: 0.25)
#' @param cutoff_pval P-value cutoff for significance (default: 0.05)
#' @param point_size Size of points (default: 1.5)
#' @param point_alpha Transparency of points (default: 1)
#' @param label_size Size of labels (default: 3.5)
#' @param color_up Color for upregulated genes (default: "#E64B35")
#' @param color_down Color for downregulated genes (default: "#4DBBD5")
#' @param color_ns Color for non-significant genes (default: "grey70")
#' @param xlim_manual Manually set x-axis limits (default: NULL, auto-calculated)
#' @param seed Random seed for reproducible label placement (default: 126)
#'
#' @return A ggplot object
plot_volcano <- function(
    df_deg,
    genes_to_label,
    title,
    col_p = "pvals_adj",
    col_log2fc = "logfoldchanges",
    col_names = "names",
    cutoff_log2fc = 0.25,
    cutoff_pval = 0.05,
    point_size = 1.5,
    point_alpha = 1,
    label_size = 3.5,
    color_up = "#E64B35",
    color_down = "#4DBBD5",
    color_ns = "grey70",
    xlim_manual = NULL,
    seed = 126) {
  library(dplyr)
  library(ggplot2)
  library(ggrepel)

  # Set seed for reproducible label placement
  set.seed(seed)

  # Prepare data
  df_volcano <- df_deg %>%
    rename(
      p = !!col_p,
      logfoldchanges = !!col_log2fc,
      names = !!col_names
    ) %>%
    mutate(
      neg_log10_pval = -log10(p),
      significant = case_when(
        p < cutoff_pval & logfoldchanges > cutoff_log2fc ~ "Up",
        p < cutoff_pval & logfoldchanges < -cutoff_log2fc ~ "Down",
        TRUE ~ "Not Sig"
      ),
      label = ifelse(names %in% genes_to_label, names, "")
    ) %>%
    arrange(desc(.data$logfoldchanges))

  # Print labeled genes for verification
  labeled_genes <- df_volcano %>%
    filter(.data$label != "") %>%
    pull(.data$label)
  cat("Labeled genes (", length(labeled_genes), "):\n", sep = "")
  print(labeled_genes)

  # Calculate x-axis limits
  if (is.null(xlim_manual)) {
    xlim <- max(abs(df_volcano$logfoldchanges), na.rm = TRUE)
  } else {
    xlim <- xlim_manual
  }

  # Create volcano plot
  p_volcano <- ggplot(
    df_volcano,
    aes(
      x = .data$logfoldchanges,
      y = .data$neg_log10_pval,
      color = .data$significant
    )
  ) +
    geom_point(alpha = point_alpha, size = point_size) +
    scale_color_manual(
      values = c("Up" = color_up, "Down" = color_down, "Not Sig" = color_ns)
    ) +
    geom_hline(
      yintercept = -log10(cutoff_pval),
      linetype = "dashed", color = "grey30"
    ) +
    geom_vline(
      xintercept = c(-cutoff_log2fc, cutoff_log2fc),
      linetype = "dashed", color = "grey30"
    ) +
    geom_text_repel(
      aes(label = .data$label),
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.5,
      segment.size = 0.3,
      min.segment.length = 0,
      size = label_size,
      fontface = "italic",
      show.legend = FALSE,
      direction = "both",
      seed = seed
    ) +
    scale_x_continuous(limits = c(-xlim, xlim)) +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 P-value"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )

  return(p_volcano)
}

# %% LMP1+ vs LMP1- Tumor ==========
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
data_dir <- file.path(dirname(script_path), "src", "07_figure_7_neighborhood_analysis")

df_deg <- read_csv(fs::path(data_dir, "02_deg_lmp1pos_vs_lmp1neg_tumors.csv"))
genes_to_label <- c(
  "CCL22", "CD40", "TRAF1", "LMP1", "TNIP1", "EBI3", "STAT5A", "BCL2L1",
  "NFKBIA", "IGHM", "IGHG1/2", "MZB1", "PTEN", "NDUFA3", "UCP2", "FBXL4",
  "COX20", "CD79B"
)
print(setdiff(genes_to_label, df_deg$names))


p_tumor_pval <- plot_volcano(
  df_deg = df_deg,
  genes_to_label = genes_to_label,
  title = "Volcano Plot: LMP1+ vs LMP1- Tumor",
  col_p = "pvals",
  cutoff_log2fc = 0.25
)


# %% Macrophages Near LMP1+ vs LMP1- Volcano ==========
df_deg <- read_csv(fs::path(data_dir, "02_deg_mac_near_lmp1pos_vs_lmp1neg.csv"))
genes_to_label <- c(
  "CTSB", "LYZ", "IDO1", "C1S", "GPNMB", "LGALS9", "CYBB", "PTEN", "COX1",
  "COX2", "NDUFA3", "MCM2", "IGHG1/2", "CD79B"
)
print(setdiff(genes_to_label, df_deg$names))

p_mac_pval <- plot_volcano(
  df_deg = df_deg,
  genes_to_label = genes_to_label,
  title = "Volcano Plot: Macrophages Near LMP1+ vs LMP1- Tumors",
  col_p = "pvals",
  cutoff_log2fc = 0.25
)

# %% Save Volcano Plots ==========
p <- patchwork::wrap_plots(p_tumor_pval, p_mac_pval, ncol = 2)

fs::dir_create(fs::path(data_dir, "figures"))
ggsave(
  fs::path(data_dir, "figures", "01_volcano_plot_deg_combined.pdf"),
  p,
  width = 16,
  height = 20
)

# %% Define plot_gsea_barplot ==========
#' Create GSEA barplot from enrichment results
#'
#' @param df_gsea Data frame containing GSEA results
#' @param terms_to_plot Character vector of terms to plot
#' @param title Plot title (default: NULL, no title)
#' @param col_term Column name for terms (default: "Term")
#' @param col_nes Column name for normalized enrichment score (default: "NES")
#' @param col_p Column name for p-values/FDR (default: "FDR q-val")
#' @param cutoff_p P-value cutoff for significance marking (default: 0.25)
#' @param sig_shift Shift factor for significance markers (default: 1.05)
#' @param sig_marker Marker for significant terms (default: "*")
#' @param sig_size Size of significance markers (default: 5)
#' @param color_low Color for low p-values (significant, default: "#e74c3c" red)
#' @param color_high Color for high p-values (not significant, default: "#3498db" blue)
#' @param bar_alpha Transparency of bars (default: 1)
#' @param xlim_manual Manually set x-axis limits (default: NULL, auto-calculated)
#' @param reverse_order Reverse the order of terms (default: FALSE, low NES at bottom)
#'
#' @return A ggplot object
plot_gsea_barplot <- function(
    df_gsea,
    terms_to_plot,
    title = NULL,
    col_term = "Term",
    col_nes = "NES",
    col_p = "FDR q-val",
    cutoff_p = 0.25,
    sig_shift = 1.05,
    sig_marker = "*",
    sig_size = 5,
    color_low = "#e74c3c",
    color_high = "#3498db",
    bar_alpha = 1,
    xlim_manual = NULL,
    reverse_order = FALSE) {
  library(dplyr)
  library(ggplot2)

  # Prepare data
  df_plot <- df_gsea %>%
    filter(!!sym(col_term) %in% terms_to_plot)

  # Arrange by NES
  if (reverse_order) {
    df_plot <- df_plot %>% arrange(desc(!!sym(col_nes)))
  } else {
    df_plot <- df_plot %>% arrange(!!sym(col_nes))
  }

  # Set factor levels for ordering
  df_plot[[col_term]] <- factor(df_plot[[col_term]], levels = df_plot[[col_term]])

  # Add p-value and significance columns
  df_plot[["p"]] <- df_plot[[col_p]]
  df_plot[["significant"]] <- df_plot[["p"]] < cutoff_p

  # Print significant terms for verification
  sig_terms <- df_plot %>%
    filter(.data$significant) %>%
    pull(!!sym(col_term))
  cat("Significant terms (", length(sig_terms), "):\n", sep = "")
  print(sig_terms)

  # Calculate x-axis limits
  if (is.null(xlim_manual)) {
    xlim <- max(abs(df_plot[[col_nes]]), na.rm = TRUE) * sig_shift
  } else {
    xlim <- xlim_manual
  }

  # Create barplot
  p <- ggplot(df_plot, aes(x = !!sym(col_nes), y = !!sym(col_term), fill = p)) +
    geom_col(alpha = bar_alpha) +
    scale_fill_gradient(
      low = color_low,
      high = color_high,
      name = col_p
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      color = "black",
      linewidth = 0.5
    ) +
    geom_text(
      aes(
        x = !!sym(col_nes) * sig_shift,
        label = ifelse(.data$significant, sig_marker, "")
      ),
      size = sig_size
    ) +
    labs(
      title = title,
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    scale_x_continuous(limits = c(-xlim, xlim)) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      panel.grid.minor = element_blank()
    )

  return(p)
}


# %% Macrophages Near LMP1+ vs LMP1- ==========

df_gsea_mac <- read_csv(fs::path(data_dir, "04_gsea_mac_all_merged.csv"))
terms_to_plot <- c(
  "Ferroptosis WP4313",
  "Cellular Response To Starvation R-HSA-9711097",
  "Terminal Pathway Of Complement R-HSA-166665",
  "ROS And RNS Production In Phagocytes R-HSA-1222556",
  "Negative Regulation Of Leukocyte Activation (GO:0002695)",
  "Negative Regulation Of T Cell Proliferation (GO:0042130)",
  "Negative Regulation Of Phosphatidylinositol 3-Kinase Signaling (GO:0014067)",
  "Negative Regulation Of Focal Adhesion Assembly (GO:0051895)",
  "Inositol Phosphate Metabolism R-HSA-1483249",
  "CD22 Mediated BCR Regulation R-HSA-5690714",
  "Scavenging Of Heme From Plasma R-HSA-2168880",
  "Electron Transport Chain OXPHOS System In Mitochondria WP111"
)
print(setdiff(terms_to_plot, df_gsea_mac$Term))

p_mac_fdr <- plot_gsea_barplot(
  df_gsea = df_gsea_mac,
  terms_to_plot = terms_to_plot,
  title = "GSEA: Macrophages Near LMP1+ vs LMP1-",
  col_p = "FDR q-val",
  cutoff_p = 0.25
)

p_mac_p <- plot_gsea_barplot(
  df_gsea = df_gsea_mac,
  terms_to_plot = terms_to_plot,
  title = "GSEA: Macrophages Near LMP1+ vs LMP1-",
  col_p = "NOM p-val",
  cutoff_p = 0.05
)

# %% LMP1+ vs LMP1- Tumor ==========

df_gsea_tumor <- read_csv(fs::path(data_dir, "04_gsea_tumor_all_merged.csv"))
terms_to_plot <- c(
  "STAT5 Activation Downstream Of FLT3 ITD Mutants R-HSA-9702518",
  "Positive Regulation Of Telomerase Activity (GO:0051973)",
  "Tumor Necrosis Factor-Mediated Signaling Pathway (GO:0033209)",
  "Interleukin-27 Signaling R-HSA-9020956",
  "Ebstein Barr Virus LMP1 Signaling WP262",
  "NF-kappa B signaling pathway",
  "Basement Membrane Organization (GO:0071711)",
  "Negative Regulation Of Calcium Ion Transmembrane Transport (GO:1903170)",
  "Glycosaminoglycan Metabolic Process (GO:0030203)",
  "CD22 Mediated BCR Regulation R-HSA-5690714",
  "Positive Regulation Of Calcium Ion Import Across Plasma Membrane (GO:1905665)",
  "Sulfatase And Aromatase Pathway WP5368"
)
print(setdiff(terms_to_plot, df_gsea_tumor$Term))

p_tumor_fdr <- plot_gsea_barplot(
  df_gsea = df_gsea_tumor,
  terms_to_plot = terms_to_plot,
  title = "GSEA: LMP1+ vs LMP1- Tumor Cells",
  col_p = "FDR q-val",
  cutoff_p = 0.25
)

p_tumor_p <- plot_gsea_barplot(
  df_gsea = df_gsea_tumor,
  terms_to_plot = terms_to_plot,
  title = "GSEA: LMP1+ vs LMP1- Tumor Cells",
  col_p = "NOM p-val",
  cutoff_p = 0.05
)


# %% Save GSEA Plots ==========
p <- patchwork::wrap_plots(p_tumor_fdr, p_tumor_p, p_mac_fdr, p_mac_p, ncol = 1)

fs::dir_create(fs::path(data_dir, "figures"))
ggsave(
  fs::path(data_dir, "figures", "02_gsea_combined.pdf"),
  p,
  width = 10,
  height = 12
)


