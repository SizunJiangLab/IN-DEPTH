# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(qs)
# setworking directory
wdpath <- c("./data/SGCC_data/")
from.1 = 6
to.1 = 14
by.1 = 1
radius.range <- seq(from = from.1,to= to.1, by = by.1 )
plot.list <- list()
for (i in 1:length(radius.range)){
  # Set dimensions
  x_border <- 60
  
  # Create a base scatter plot
  base_data <- expand.grid(x = 1:x_border, y = 1:x_border)
  
  # Define circle properties
  radius <- radius.range[i]
  center_distance <- 30  # Distance of each circle's center from the middle line
  
  # Generate movement data
  movements <- data.frame(
    x_shift = seq(from = center_distance, to = 0, length.out = 10),
    circle1_x = center_distance - seq(from = center_distance/2, to = 0, length.out = 10),
    circle2_x = center_distance + seq(from = center_distance/2, to = 0, length.out = 10)
  )
  
  # Add circle data to base_data for each movement
  plot_data <- lapply(1:nrow(movements), function(i) {
    movement <- movements[i,]
    base_data %>%
      mutate(
        inside_circle1 = sqrt((x - movement$circle1_x)^2 + (y - x_border/2)^2) <= radius,
        inside_circle2 = sqrt((x - movement$circle2_x)^2 + (y - x_border/2)^2) <= 1.5*radius,
        movement = paste("Movement", i)
      )
  }) %>%
    bind_rows()
  
  # Create a new variable for coloring
  plot_data <- plot_data %>%
    mutate(circle_overlap = case_when(
      inside_circle1 & !inside_circle2 ~ "Red",    # Only in the first circle
      !inside_circle1 & inside_circle2 ~ "Blue",   # Only in the second circle
      inside_circle1 & inside_circle2 ~ "Purple",  # In both circles (if possible)
      TRUE ~ "Black"                               # In neither circle
    ))
  plot_data$inside_circle2[plot_data$inside_circle1 & plot_data$inside_circle2 == T] <- FALSE
  # Plotting
  plot_data$movement <- factor(plot_data$movement,levels = paste0("Movement ",1:10))
  p <- ggplot(plot_data, aes(x = x, y = y, color = circle_overlap)) +
    geom_point(size = 1.5) +
    facet_wrap(~ movement, nrow = 1) +
    scale_color_manual(values = c("Red" = "#16964a", "Blue" = "#2958a8", "Purple" = "#16964a", "Black" = "grey")) +
    theme_void() +
    theme(panel.margin = unit(-0.5, "lines")) +
    theme(plot.margin = unit(c(0,0,0,0), "lines")) +
    theme(strip.background = element_blank()) +
    theme(plot.background = element_blank()) +
    theme(strip.text = element_blank()) +
    theme(axis.ticks.margin = unit(0, "lines"))+
    theme(legend.position = "none")
  plot.list[[i]] <- p
  # Split the dataframe into a list based on the 'movement' column
  list_of_movements <- split(plot_data, plot_data$movement)
  save.name <- paste0("border","-",x_border,"_",
                      "radius","-",radius)
  qsave(list_of_movements,file = file.path(wdpath,"Simulation_Data","MovingPattern",paste0("data_",save.name,".qs")))
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
ggsave(plot = combined_plot,filename = paste0("Supplementary Figure 3C.svg"),
       device = "svg",width = 20,height = 18)  


