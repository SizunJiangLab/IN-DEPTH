# Load required libraries
library(ggplot2)
library(dplyr)
library(qs)
library(ggpubr)
# set working directory
wdpath <- c("./data/SGCC_data/")
radius.range <- seq(from = 2.5,to= 20, by = 2.5)
plot.list <- list()
for (i in 1:length(radius.range)){
  # Set dimensions
  x_border <- 60
  
  # Create a base scatter plot
  base_data <- expand.grid(x = 1:x_border, y = 1:x_border)
  
  # Define circle properties
  radius <- radius.range[i]
  center_x <- x_border / 2
  center_y <- x_border / 2
  
  # Generate movement data for the outer ring to close in
  movements <- data.frame(
    x_shift = seq(from = 40, to = 1.2*radius, length.out = 10)
  )
  
  # Add circle data to base_data for each movement
  plot_data <- lapply(1:nrow(movements), function(i) {
    r_outer <- movements$x_shift[i]
    base_data %>%
      mutate(
        inside_circle1 = sqrt((x - center_x)^2 + (y - center_y)^2) <= radius,
        inside_circle2_outer = sqrt((x - center_x)^2 + (y - center_y)^2) <= r_outer,
        inside_ring = inside_circle2_outer & !inside_circle1,
        movement = paste("Movement", i)
      )
  }) %>%
    bind_rows()
  
  # Create a new variable for coloring
  plot_data <- plot_data %>%
    mutate(circle_overlap = case_when(
      inside_circle1 ~ "Solid Circle",        # In the solid circle
      inside_ring ~ "Concentric Ring",        # In the concentric ring
      TRUE ~ "Outside"                        # Outside both
    ))
  # save name
  save.name <- paste0("border","-",x_border,"_",
                      "radius","-",radius)
  # Check the first few rows of data and the structure for the first movement
  # plot_data[1:5,]
  list_of_movements <- split(plot_data, plot_data$movement)
  qsave(list_of_movements,file =  file.path(wdpath,"Simulation_Data","RingPattern",paste0("data_",save.name,".qs")))
  # Plotting
  plot_data$movement <- factor(plot_data$movement, levels = paste0("Movement ", 1:10))
  p <- ggplot(plot_data, aes(x = x, y = y, fill = circle_overlap)) +
    geom_tile(alpha =0.95) +
    facet_wrap(~ movement, nrow = 1) +
    scale_fill_manual(values = c("Solid Circle" = "#16964a", "Concentric Ring" = "#2958a8", "Outside" = "grey")) +
    #scale_fill_manual(values = c("Solid Circle" = "#a9669b", "Concentric Ring" = "#6da061", "Outside" = "grey")) +
    theme_void() +
    theme(panel.margin = unit(-0.5, "lines")) +
    theme(plot.margin = unit(c(0,0,0,0), "lines")) +
    theme(strip.background = element_blank()) +
    theme(plot.background = element_blank()) +
    theme(strip.text = element_blank()) +
    theme(axis.ticks.margin = unit(0, "lines"))+
    theme(legend.position = "none")
  plot.list[[i]] <- p
}




combined_plot <- ggarrange(plot.list[[1]],
                           NULL,
                           plot.list[[2]],
                           NULL,
                           plot.list[[3]],
                           NULL,
                           plot.list[[4]],
                           NULL,
                           plot.list[[5]],
                           NULL,
                           plot.list[[6]],
                           NULL,
                           plot.list[[7]],
                           NULL,
                           plot.list[[8]],
                           nrow = 15,
                           ncol = 1,
                           heights = rep(c(1,-0.06),times =8)[-16])
combined_plot   
ggsave(plot = combined_plot,filename = "Figure 3C.svg",
       device = "svg",width = 15,height = 12)  
     