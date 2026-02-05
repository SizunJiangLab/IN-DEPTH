library(dplyr)
library(igraph)
library(deldir)
library(tidyr)
library(tidyverse)
library(ggplot2)

generate_edge_data <- function(annotation, cell_types = c("CD8 T", "M2"), pixel_to_um = 0.5) {
   # Prepare cell data
  cell_data <- annotation %>%
    select(cellLabel, X_cent_fusion, Y_cent_fusion, Annotation) %>%
    dplyr::rename(x = X_cent_fusion, y = Y_cent_fusion, domain = Annotation)
  
  # Ensure the specified cell types are present in the data
  if (!all(cell_types %in% unique(cell_data$domain))) {
    stop("One or both specified cell types are not present in the data.")
  }
  
  # Perform Delaunay triangulation
  delaunay_result <- deldir(cell_data$x, cell_data$y)
  
  # Extract the Delaunay edges (neighbor relationships) without filtering
  edges <- delaunay_result$delsgs %>%
    mutate(
      distance_um = sqrt((x2 - x1)^2 + (y2 - y1)^2) * pixel_to_um,
      ind1 = as.character(cell_data$cellLabel[ind1]),
      ind2 = as.character(cell_data$cellLabel[ind2])
    )
  
  # Create a unique identifier for each edge (smallest to largest to account for undirected nature)
  edges <- edges %>%
    mutate(
      edge_id = paste(pmin(ind1, ind2), pmax(ind1, ind2), sep = "_")
    )
  
  # Create an igraph object using all nodes and edges
  graph <- graph_from_data_frame(d = edges[, c("ind1", "ind2")], directed = FALSE, vertices = cell_data)
  
  # Add domain information to the graph
  V(graph)$domain <- cell_data$domain[match(V(graph)$name, cell_data$cellLabel)]
  
  # Classify the edges based on the two specified cell types
  E(graph)$edge_type <- apply(as.data.frame(ends(graph, E(graph))), 1, function(e) {
    domain1 <- V(graph)$domain[match(e[1], V(graph)$name)]
    domain2 <- V(graph)$domain[match(e[2], V(graph)$name)]
    
    # Check if both ends belong to cell type 1
    if (domain1 == domain2 && domain1 == cell_types[1]) {
      return("within_type1")  # Edges within domain 1
    }
    # Check if both ends belong to cell type 2
    else if (domain1 == domain2 && domain1 == cell_types[2]) {
      return("within_type2")  # Edges within domain 2
    }
    # Check if ends belong to different cell types
    else if (domain1 != domain2 && (domain1 %in% cell_types && domain2 %in% cell_types)) {
      return("across")  # Edges between two domains
    } else {
      return(NA)  # Ignore edges outside the specified domains
    }
  })
  
  # Filter out edges that are not classified (NA)
  valid_edges <- which(!is.na(E(graph)$edge_type))
  graph <- subgraph.edges(graph, valid_edges, delete.vertices = FALSE)
  
  # Create a unique identifier for edges in the graph (to match with the original edges)
  graph_edges <- data.frame(
    edge_id = paste(
      pmin(ends(graph, E(graph))[, 1], ends(graph, E(graph))[, 2]),
      pmax(ends(graph, E(graph))[, 1], ends(graph, E(graph))[, 2]),
      sep = "_"
    )
  )
  
  # Match distances from the original edges
  E(graph)$distance_um <- edges$distance_um[match(graph_edges$edge_id, edges$edge_id)]
  
  # Apply the distance filter to retain edges shorter than 50 µm
  graph <- subgraph.edges(graph, which(E(graph)$distance_um < 50), delete.vertices = FALSE)
  
  # Prepare edge data for plotting
  edge_data <- data.frame(
    x1 = cell_data$x[match(ends(graph, E(graph))[, 1], cell_data$cellLabel)],
    y1 = cell_data$y[match(ends(graph, E(graph))[, 1], cell_data$cellLabel)],
    x2 = cell_data$x[match(ends(graph, E(graph))[, 2], cell_data$cellLabel)],
    y2 = cell_data$y[match(ends(graph, E(graph))[, 2], cell_data$cellLabel)],
    edge_type = E(graph)$edge_type,
    distance = E(graph)$distance_um
  )
  
  # Return the edge data and cell data
  list(edge_data = edge_data, cell_data = cell_data)
  
}

filter_samples_summary <- function(annotation, 
                                   max_diff = 0.8, 
                                   cell_types = c("CD8T", "M1"), 
                                   tumor_annotation = "Tumor",  # Specify the annotation used to identify tumor cells
                                   tumor_upper = 0.5,
                                   tumor_lower = 0.2,
                                   logfc=1.25) {

  # Check if the specified cell types are present in the data
  if (!all(cell_types %in% unique(annotation$Annotation))) {
    stop("One or both specified cell types are not present in the data.")
  }
  
  # Calculate the total number of cells per sample
  total_counts <- annotation %>%
    group_by(coreName) %>%
    summarise(total_cells = n(), .groups = "drop")
  
  # Calculate the number of cells of each specified type per sample
  cell_type_counts <- annotation %>%
    filter(Annotation %in% cell_types) %>%
    group_by(coreName, Annotation) %>%
    summarise(cell_type_count = n(), .groups = "drop")
  
  # Merge total cell counts to compute relative proportions
  cell_type_proportions <- cell_type_counts %>%
    left_join(total_counts, by = "coreName") %>%
    mutate(relative_proportion = cell_type_count / total_cells)
  
  # Reshape data to wide format for easier difference calculation
  CT1_count_name <- paste0("cell_type_count_", cell_types[1])
  CT2_count_name <- paste0("cell_type_count_", cell_types[2])
  CT1_prop_name <- paste0("relative_proportion_", cell_types[1])
  CT2_prop_name <- paste0("relative_proportion_", cell_types[2])
  cell_type_wide <- cell_type_proportions %>%
    pivot_wider(
      names_from = Annotation,
      values_from = c(cell_type_count, relative_proportion),
      names_sep = "_"
    ) %>%
    dplyr::rename(
      CT1_count = CT1_count_name,
      CT2_count = CT2_count_name,
      CT1_proportion = CT1_prop_name,
      CT2_proportion = CT2_prop_name
    ) %>%
    select(coreName, total_cells, CT1_count, CT2_count, CT1_proportion, CT2_proportion)
  
  
  # Compute the absolute difference in relative proportions for each sample between the two cell types
  cell_type_wide <- cell_type_wide %>%
    mutate(diff = abs(CT1_proportion - CT2_proportion),
           log_ratio = log(CT1_proportion/CT2_proportion),
           )
  
  # Calculate the number of tumor cells in each sample
  tumor_counts <- annotation %>%
    filter(Annotation == tumor_annotation) %>%
    group_by(coreName) %>%
    summarise(tumor_cells = n(), .groups = "drop")
  
  # Merge with total cell counts to calculate the proportion of tumor cells
  tumor_proportions <- tumor_counts %>%
    left_join(total_counts, by = "coreName") %>%
    mutate(tumor_proportion = tumor_cells / total_cells)
  
  # Combine data for final table
  final_summary <- total_counts %>%
    left_join(cell_type_wide, by = "coreName") %>%
    left_join(tumor_proportions, by = "coreName") %>%
    mutate(filtered_out = ifelse(diff > max_diff | 
                                   tumor_proportion > tumor_upper |
                                   tumor_proportion < tumor_lower |
                                   abs(log_ratio) > logfc, TRUE, FALSE))
  
  # Return the final summary table
  return(final_summary)
}





# plot function
plot_cell_data_with_colors <- function(edge_data = NULL, cell_data = NULL, 
                                       cell_type_colors = cell_type_colors,
                                       cell_types = c("CD8 T", "CD4 T"),
                                       size_domain = 6, size_other = 3,corename = NULL,...) {
 
  
  # Check if cell types in cell_data match the provided cell_types
  if (!all(cell_types %in% unique(cell_data$domain))) {
    stop("One or both specified cell types are not present in the provided cell data.")
  }
  
  # Assign size to cells based on domain
  cell_data$size <- ifelse(cell_data$domain %in% cell_types, size_domain, size_other)
  
  # Assign colors to cells: use color table for domain 1 and 2, grey for others
  cell_data$fill_color <- sapply(cell_data$domain, function(domain) {
    if (domain %in% cell_types) {
      return(cell_type_colors[domain])  # Use color for domain 1 and 2
    } else {
      return("grey")  # Grey color for other cell types
    }
  })
  
  # Assign colors to edges: use the color of the domain for within-domain edges, and black for across-domain edges
  edge_data$edge_color <- sapply(edge_data$edge_type, function(edge_type) {
    if (edge_type == "across") {
      return("black")  # Across domain edges are black
    } else if (edge_type == "within_type1") {
      return(cell_type_colors[cell_types[1]])  # Color for the first specified cell type
    } else if (edge_type == "within_type2") {
      return(cell_type_colors[cell_types[2]])  # Color for the second specified cell type
    } else {
      return("grey")  # Default color if not classified correctly
    }
  })
  
  # Plotting using ggplot2
  p <- ggplot() +
    # Plot edges with assigned colors
    geom_segment(data = edge_data, 
                 aes(x = x1, y = -y1, xend = x2, yend = -y2, color = edge_color), 
                 size = 1) +
    # Plot cells with customized sizes and colors
    geom_point(data = cell_data, 
               aes(x = x, y = -y, fill = fill_color, size = size), 
               shape = 21, color = "white") +
    # Customize colors for edges
    scale_color_identity() +  # Use identity scale to keep assigned colors
    # Use the extracted color mapping for cell types, including grey for others
    scale_fill_identity() +   # Use identity scale to keep assigned fill colors
    # Control the scale of the size of the points
    scale_size_identity() #+
    # labs(title = corename,
    #      x = "X Coordinate", y = "Y Coordinate", color = "Edge Type")
  
  # Print the plot
  return(p)
}
# calculate odds and proportion
calculate_odds <- function(cell_data, edge_data, ct1_name, ct2_name) {
  # Check if specified cell types are present in the dataset
  if (!(ct1_name %in% cell_data$domain)) {
    stop(paste("The specified CT1 type:", ct1_name, "is not present in the cell data."))
  }
  if (!(ct2_name %in% cell_data$domain)) {
    stop(paste("The specified CT2 type:", ct2_name, "is not present in the cell data."))
  }
  
  # Extract the number of each cell type
  num_ct1 <- nrow(cell_data[cell_data$domain == ct1_name,])
  num_ct2 <- nrow(cell_data[cell_data$domain == ct2_name,])
  
  # Calculate total possible connections
  total_across <- num_ct1 * num_ct2
  total_within_ct1 <- num_ct1 * (num_ct1 - 1) / 2
  total_within_ct2 <- num_ct2 * (num_ct2 - 1) / 2
  
  # Count the actual connections from edge_data
  actual_across <- nrow(edge_data[edge_data$edge_type == "across",])
  mean_distance_actual_across <- mean(edge_data$distance[edge_data$edge_type == "across"])
  sd_distance_actual_across <- sd(edge_data$distance[edge_data$edge_type == "across"])
  actual_within_ct1 <- nrow(edge_data[edge_data$edge_type == "within_type1", ])
  actual_within_ct2 <- nrow(edge_data[edge_data$edge_type == "within_type2", ])
  mean_distance_actual_within_ct1 <- mean(edge_data$distance[edge_data$edge_type == "within_type1"])
  sd_distance_actual_within_ct1 <- sd(edge_data$distance[edge_data$edge_type == "within_type1"])
  mean_distance_actual_within_ct2  <- mean(edge_data$distance[edge_data$edge_type == "within_type2"])
  sd_distance_actual_within_ct2  <- sd(edge_data$distance[edge_data$edge_type == "within_type2"])
  
  # Calculate odds
  odds_across <- actual_across / (total_across - actual_across)
  odds_within_ct1 <- actual_within_ct1 / (total_within_ct1 - actual_within_ct1)
  odds_within_ct2 <- actual_within_ct2 / (total_within_ct2 - actual_within_ct2)
  #
  total_edges <- actual_across + actual_within_ct1 + actual_within_ct2
  # Calculate proportion
  pp_across <- actual_across / total_edges
  pp_within_ct1 <- actual_within_ct1 / total_edges
  pp_within_ct2 <- actual_within_ct2 / total_edges
  # Return results as a list
  list(
    # CT1 = ct1_name,
    # CT2 = ct2_name,
    CT1_n = num_ct1,
    CT2_n = num_ct2,
    total_n = nrow(cell_data),
    edge_across = actual_across,
    mean_distance_edge_across = mean_distance_actual_across,
    sd_distance_edge_across = sd_distance_actual_across,
    mean_distance_actual_within_ct1 = mean_distance_actual_within_ct1,
    sd_distance_actual_within_ct1 = sd_distance_actual_within_ct1,
    mean_distance_actual_within_ct2 = mean_distance_actual_within_ct2,
    sd_distance_actual_within_ct2 = sd_distance_actual_within_ct2,
    edge_within_ct1 = actual_within_ct1,
    edge_within_ct2 = actual_within_ct2,
    total_edges =  total_edges,
    odds_across = odds_across,
    odds_within_ct1 = odds_within_ct1,
    odds_within_ct2 = odds_within_ct2,
    pp_across = pp_across,
    pp_within_ct1 = pp_within_ct1,
    pp_within_ct2 = pp_within_ct2
  )
}
# plot tri cell types function
plot_tri_cell_data <- function(cell_data = NULL, 
                               cell_type_colors = cell_type_colors,
                               cell_types = c("CD4T", "Tumor","Macro"),
                               size_domain = 6, size_other = 3,corename = NULL,...) {
  # Check if cell types in cell_data match the provided cell_types
  if (!all(cell_types %in% unique(cell_data$domain))) {
    stop("One or both specified cell types are not present in the provided cell data.")
  }
  
  # Assign size to cells based on domain
  cell_data$size <- ifelse(cell_data$domain %in% cell_types, size_domain, size_other)
  
  # Assign colors to cells: use color table for domain 1 and 2, grey for others
  cell_data$fill_color <- sapply(cell_data$domain, function(domain) {
    if (domain %in% cell_types) {
      return(cell_type_colors[domain])  # Use color for domain 1 and 2
    } else {
      return("grey")  # Grey color for other cell types
    }
  })

  
  # Plotting using ggplot2
  p <- ggplot() +
    # Plot cells with customized sizes and colors
    geom_point(data = cell_data, 
               aes(x = x, y = -y, fill = fill_color, size = size), 
               shape = 21, color = "white") +
    # Customize colors for edges
    scale_color_identity() +  # Use identity scale to keep assigned colors
    # Use the extracted color mapping for cell types, including grey for others
    scale_fill_identity() +   # Use identity scale to keep assigned fill colors
    # Control the scale of the size of the points
    scale_size_identity() #+
  # labs(title = corename,
  #      x = "X Coordinate", y = "Y Coordinate", color = "Edge Type")
  
  # Print the plot
  return(p)
}
