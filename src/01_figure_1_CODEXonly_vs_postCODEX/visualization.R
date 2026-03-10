# Utils #####
## Total count barplot #####
plot_barplot <- function(df_p, y_log = TRUE, order_bar, color_bar) {
  if (y_log) {
    df_p$y <- log1p(df_p$y)
    y_label <- "Total Count (Log1p)"
  } else {
    power_10 <- floor(log10(max(df_p$y))) - 1
    df_p$y <- df_p$y / 10 ^ power_10
    y_label <- stringr::str_glue('Total Count ({sprintf("%.0e", 10 ^ power_10)})')
  }
  p <- df_p %>% 
    dplyr::mutate(dplyr::across(x, ~ factor(.x, levels = order_bar))) %>% 
    ggplot(aes(x = x, y = y, fill = x)) +
    geom_col() + 
    geom_col(width = 0.75) +
    facet_grid(~ type, scales = "free") +
    scale_fill_manual(values = color_bar) +
    labs(x = "Group", y = y_label) +
    theme_bw() +
    theme(
      legend.position = "none", 
      strip.background = element_blank(), 
      strip.text = element_text(size = 11)
    )
  return(p)
}

plot_total_count <- function(df_data) {
  color_bar <- c("only" = "#168D46", "post" = "#653575")
  df_p <- df_data %>% 
    dplyr::group_by(group, type) %>% 
    dplyr::summarise(n = sum(n), .groups = "drop") %>% 
    dplyr::rename(x = group, y = n, type = type)
  order_bar <- names(color_bar)
  p_n <- plot_barplot(df_p, y_log = FALSE, order_bar, color_bar)
  p_log <- plot_barplot(df_p, y_log = TRUE, order_bar, color_bar)
  p <- patchwork::wrap_plots(p_n, p_log, ncol = 1)
  return(p)
}

get_xyrange <- function(x, y) {
  xmax <- max(x)
  xmin <- min(x)
  ymax <- max(y)
  ymin <- min(y)
  x0 <- mean(c(xmin, xmax))
  y0 <- mean(c(ymin, ymax))
  range <- max(c(xmax - xmin, ymax - ymin)) / 2
  xyrange <- list(x = c(x0 - range, x0 + range), y = c(y0 - range, y0 + range))
  return(xyrange)
}

format_pvalue <- function(pvalue) {
  dplyr::case_when(
    pvalue < 2.20e-16 ~ "p < 2.20e-16", 
    pvalue > 0.001 ~ sprintf("p = %.3f", round(pvalue, 3)),
    T ~ paste0("p = ", format(pvalue, scientific = TRUE, digits = 3))
  )
}

## Correlation dotplot #####
plot_cor <- function(df_data, target_type) {
  color_dot <- c("gene" = "#E41A1C", "negative" = "#377EB8", "SystemControl" =  "#4DAF4A")
  
  df_p <- df_data %>% 
    dplyr::filter(type == target_type) %>% 
    dplyr::group_by(gene, group) %>% 
    dplyr::summarise(n = log1p(sum(n)), .groups = "drop") %>% 
    tidyr::pivot_wider(id_cols = gene, names_from = group, values_from = n, values_fill = 0) %>%
    dplyr::filter(only + post != 0)
  cor_res <- cor.test(df_p$only, df_p$post, method = "pearson")
  pvalue <- format_pvalue(cor_res$p.value)
  corvalue <- sprintf("pearson = %.3f", round(cor_res$estimate, 3))
  subtitle <- paste(corvalue, pvalue, sep = ", ")
  
  xyrange <- get_xyrange(df_p$only, df_p$post)
  p <- df_p %>% 
    dplyr::mutate(type = target_type) %>% 
    ggplot(aes(x = only, y = post)) +
    geom_point(aes(color = type), size = 0.6) +
    geom_smooth(method = "lm", color = "black", linewidth = .8) +
    labs(subtitle = subtitle, y = "post") + 
    scale_color_manual(values = color_dot) +
    theme_bw() +
    theme(
      legend.position = "none", 
      plot.title = element_text(face = "bold", size = 11), 
      plot.subtitle = element_text(size = 11)
    ) +
    coord_equal(xlim = xyrange$x, ylim = xyrange$y)
  return(p)
}

plot_cor_roi <- function(df_data, target_type) {
  roi_label <- unique(df_data$roi)
  res_p <- vector("list", length = length(roi_label))
  for (i in seq_along(res_p)) {
    res_p[[i]] <- df_data %>% 
      dplyr::filter(roi == roi_label[i]) %>% 
      plot_cor(target_type = target_type) + 
      labs(title = paste0("ROI = ", roi_label[i]))
  }
  x_range <- res_p %>% purrr::map(~ .x$data$only) %>% unlist()
  y_range <- res_p %>% purrr::map(~ .x$data$post) %>% unlist()
  xyrange <- get_xyrange(x_range, y_range)
  p <- res_p %>% 
    purrr::map(~ .x + coord_equal(xlim = xyrange$x, ylim = xyrange$y)) %>% 
    patchwork::wrap_plots()
  return(p)
}

plot_cor_geomx_negative <- function(df_data) {
  color_dot <- c("gene" = "#E41A1C", "negative" = "#377EB8", "SystemControl" =  "#4DAF4A")
  
  df_p <- df_data %>% 
    dplyr::filter(type == "negative") %>% 
    dplyr::group_by(roi, group) %>% 
    dplyr::summarise(n = log1p(sum(n)), .groups = "drop") %>% 
    tidyr::pivot_wider(id_cols = roi, names_from = group, values_from = n, values_fill = 0) %>%
    dplyr::filter(only + post != 0)
  cor_res <- cor.test(df_p$only, df_p$post, method = "pearson")
  pvalue <- format_pvalue(cor_res$p.value)
  corvalue <- sprintf("pearson = %.3f", round(cor_res$estimate, 3))
  subtitle <- paste(corvalue, pvalue, sep = ", ")
  
  xyrange <- get_xyrange(df_p$only, df_p$post)
  p <- df_p %>% 
    dplyr::mutate(type = "negative") %>% 
    ggplot(aes(x = only, y = post)) +
    geom_point(aes(color = type), size = 0.6) +
    geom_smooth(method = "lm", color = "black", linewidth = .8) +
    labs(subtitle = subtitle, y = "post") + 
    scale_color_manual(values = color_dot) +
    theme_bw() +
    theme(
      legend.position = "none", 
      plot.title = element_text(face = "bold", size = 11), 
      plot.subtitle = element_text(size = 11)
    ) +
    coord_equal(xlim = xyrange$x, ylim = xyrange$y)
  return(p)
}

## Corelation Heatmap #####
calculate_roi_pearson_cor <- function(tag) {
  df_data <- readr::read_csv(glue::glue("data/output/{tag}.csv"), show_col_types = FALSE)
  df_p <- df_data %>% 
    dplyr::filter(type == "gene") %>% 
    dplyr::group_by(gene, group, roi) %>% 
    dplyr::summarise(n = log1p(sum(n)), .groups = "drop") %>% 
    tidyr::pivot_wider(id_cols = c(roi, gene), names_from = group, values_from = n, values_fill = 0) %>%
    dplyr::filter(only + post != 0)
  df_p <- df_p %>% 
    dplyr::group_split(roi) %>% 
    purrr::map(~ {
      cor_res <- cor.test(.x$only, .x$post, method = "pearson")
      tibble::tibble(
        roi = unique(.x$roi), 
        pvalue = cor_res$p.value, 
        pvalue_tag = format_pvalue(pvalue), 
        corvalue = cor_res$estimate,
        corvalue_tag = sprintf("%.3f", round(corvalue, 3))
      )
    }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(
      tag = tag, 
      roi = as.factor(as.numeric(as.factor(roi)))
    )
  return(df_p)
}

heatmap_roi_pearson_cor <- function(tags, tags_label) {
  df_p <- tags %>% 
    purrr::map(calculate_roi_pearson_cor, .progress = TRUE) %>% 
    dplyr::bind_rows()
  p <- df_p %>% 
    dplyr::mutate(
      corvalue_tag = sprintf("%.2f", round(corvalue, 2)),
      dplyr::across(tag, ~ factor(.x, level = tags, labels = tags_label))
    ) %>% 
    ggplot(aes(x = roi, y = tag, fill = corvalue)) + 
    geom_tile(height = 1, width = 1, color = "black", linewidth = .75) +
    geom_text(aes(label = corvalue_tag)) +
    coord_fixed(ratio = 2/3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1)) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), axis.text.x = element_blank()) +
    labs(x = NULL, y = NULL, fill = "pearson")
  return(p)
}


# Output #####

## CosMx #####
output_figure_cosmx <- function(tag) {
  path_data <- paste0("data/output/", tag, ".csv")
  dir_output <- paste0("output/figure/figure_1/", tag)
  fs::dir_create(dir_output)
  df_data <- readr::read_csv(path_data)
  
  p_1 <- plot_total_count(df_data)
  ggsave(paste0(dir_output, "/01_total_count.pdf"), p_1, height = 7, width = 4)
  ggsave(paste0(dir_output, "/01_total_count.png"), p_1, height = 7, width = 4)
  
  p_2 <- plot_cor(df_data, target_type = "gene")
  ggsave(paste0(dir_output, "/02_cor_slide_gene.pdf"), p_2, height = 4, width = 4)
  ggsave(paste0(dir_output, "/02_cor_slide_gene.png"), p_2, height = 4, width = 4)
  
  p_3 <- plot_cor(df_data, target_type = "negative")
  ggsave(paste0(dir_output, "/03_cor_slide_negative.pdf"), p_3, height = 4, width = 4)
  ggsave(paste0(dir_output, "/03_cor_slide_negative.png"), p_3, height = 4, width = 4)
  
  p_4 <- plot_cor_roi(df_data, target_type = "gene")
  ggsave(paste0(dir_output, "/04_cor_roi_gene.pdf"), p_4, height = max(str(p_4)$t) * 4, width = max(str(p_4)$l) * 4, bg = "white", limitsize = FALSE)
  ggsave(paste0(dir_output, "/04_cor_roi_gene.png"), p_4, height = max(str(p_4)$t) * 4, width = max(str(p_4)$l) * 4, bg = "white", limitsize = FALSE)
  
  p_5 <- plot_cor_roi(df_data, target_type = "negative")
  ggsave(paste0(dir_output, "/05_cor_roi_negative.pdf"), p_5, height = max(str(p_5)$t) * 4, width = max(str(p_5)$l) * 4, bg = "white", limitsize = FALSE)
  ggsave(paste0(dir_output, "/05_cor_roi_negative.png"), p_5, height = max(str(p_5)$t) * 4, width = max(str(p_5)$l) * 4, bg = "white", limitsize = FALSE)
}

## GeoMx #####
output_figure_geomx <- function(tag) {
  path_data <- paste0("data/output/", tag, ".csv")
  dir_output <- paste0("output/figure/figure_1/", tag)
  fs::dir_create(dir_output)
  df_data <- readr::read_csv(path_data)
  
  p_1 <- plot_total_count(df_data)
  ggsave(paste0(dir_output, "/01_total_count.pdf"), p_1, height = 7, width = 4)
  ggsave(paste0(dir_output, "/01_total_count.png"), p_1, height = 7, width = 4)
  
  p_2 <- plot_cor(df_data, target_type = "gene")
  ggsave(paste0(dir_output, "/02_cor_slide_gene.pdf"), p_2, height = 4, width = 4)
  ggsave(paste0(dir_output, "/02_cor_slide_gene.png"), p_2, height = 4, width = 4)
  
  p_3 <- plot_cor_geomx_negative(df_data)
  ggsave(paste0(dir_output, "/03_cor_negative.pdf"), p_3, height = 4, width = 4)
  ggsave(paste0(dir_output, "/03_cor_negative.png"), p_3, height = 4, width = 4)
  
  p_4 <- plot_cor_roi(df_data, target_type = "gene")
  ggsave(paste0(dir_output, "/04_cor_roi_gene.pdf"), p_4, height = max(str(p_4)$t) * 4, width = max(str(p_4)$l) * 4, bg = "white", limitsize = FALSE)
  ggsave(paste0(dir_output, "/04_cor_roi_gene.png"), p_4, height = max(str(p_4)$t) * 4, width = max(str(p_4)$l) * 4, bg = "white", limitsize = FALSE)
}

## VisiumHD #####
output_figure_visiumhd <- function(tag, y_step = NULL) {
  path_data <- paste0("data/output/", tag, ".csv")
  dir_output <- paste0("output/figure/figure_1/", tag)
  fs::dir_create(dir_output)
  df_data <- readr::read_csv(path_data)
  
  p_1 <- plot_total_count(df_data)
  ggsave(paste0(dir_output, "/01_total_count.pdf"), p_1, height = 7, width = 4)
  ggsave(paste0(dir_output, "/01_total_count.png"), p_1, height = 7, width = 4)
  
  p_2 <- plot_cor(df_data, target_type = "gene")
  if (!is.null(y_step)) {
    p_2 <- p_2 + ggplot2::scale_y_continuous(breaks = scales::breaks_width(4))
  }
  ggsave(paste0(dir_output, "/02_cor_slide_gene.pdf"), p_2, height = 4, width = 4)
  ggsave(paste0(dir_output, "/02_cor_slide_gene.png"), p_2, height = 4, width = 4)
  
  p_4 <- plot_cor_roi(df_data, target_type = "gene")
  # if (!is.null(y_step)) {
  #   p_4 <- p_4 & ggplot2::scale_y_continuous(breaks = scales::breaks_width(4))
  # }
  ggsave(paste0(dir_output, "/04_cor_roi_gene.pdf"), p_4, height = max(str(p_4)$t) * 4, width = max(str(p_4)$l) * 4, bg = "white", limitsize = FALSE)
  ggsave(paste0(dir_output, "/04_cor_roi_gene.png"), p_4, height = max(str(p_4)$t) * 4, width = max(str(p_4)$l) * 4, bg = "white", limitsize = FALSE)
}


## General #####
output_figure_cor_slide_gene <- function() {
  modify_xyrange <- function(all_path) {
    list_p <- all_path %>% purrr::map(~ readr::read_csv(.x) %>% plot_cor(target_type = "gene"))
    xyrange <- list_p %>% purrr::map(~ c(.x$data$only, .x$data$post)) %>% unlist() %>% range()
    list_p <- list_p %>% purrr::map(~ .x + coord_equal(xlim = xyrange, ylim = xyrange))
    return(list_p)
  }
  
  path_cosmx <- fs::dir_ls("data/output/", regexp = "cosmx")
  p_cosmx <- modify_xyrange(path_cosmx)
  
  path_geomx <- fs::dir_ls("data/output/", regexp = "geomx")
  p_geomx <- modify_xyrange(path_geomx)
  
  p_list <- c(p_cosmx, p_geomx)
  names(p_list) %>% 
    purrr::map(~ {
      tag <- .x %>% fs::path_file() %>% fs::path_ext_remove()
      p <- p_list[[.x]]
      dir_output <- "output/figure/figure_1/cor_slide_gene_xyrange"
      fs::dir_create(dir_output)
      ggsave(fs::path(dir_output, tag, ext = "pdf"), p, height = 4, width = 4)
      ggsave(fs::path(dir_output, tag, ext = "png"), p, height = 4, width = 4)
    })
}

output_heatmap_roi_pearson_cor <- function() {
  tags <- c(
    "cosmx_1k", "cosmx_6k", "geomx_DLBCL", "geomx_LN", "geomx_RCC", 
    "geomx_USC", "geomx_cycling", "geomx_fusion", "visiumhd_only=jiang_post=postcodex_cpb10000"
  )
  tags_label <- c(
    "cosmx_1k", "cosmx_6k", "geomx_DLBCL", "geomx_LN", "geomx_RCC", 
    "geomx_USC", "geomx_cycling", "geomx_fusion", "visiumhd"
  )
  p <- heatmap_roi_pearson_cor(tags, tags_label)
  ggsave(
    "output/figure/cor_heatmap/heatmap_roi_pearson_cor.pdf", 
    p, 
    height = length(unique(p$data$tag)) / 2, 
    width = length(unique(p$data$roi)) / 2
  )
  ggsave(
    "output/figure/cor_heatmap/heatmap_roi_pearson_cor.png", 
    p, 
    height = length(unique(p$data$tag)) / 2, 
    width = length(unique(p$data$roi)) / 2
  )
}
