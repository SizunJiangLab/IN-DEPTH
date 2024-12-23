# install packages
install_and_load <- function(packages) {
  # Ensure that the 'remotes' package is installed and loaded
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  library(remotes)
  
  # Loop through each package and install from the appropriate source
  for (pkg in names(packages)) {
    package_name <- pkg
    source_url <- packages[pkg]
    
    # Check if the package is already installed
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
      # Install from GitHub if URL is provided
      if (grepl("github", source_url)) {
        cat(sprintf("Installing %s from GitHub repository %s\n", package_name, source_url))
        remotes::install_github(source_url)
      } else {  # Install from CRAN
        cat(sprintf("Installing package: %s from CRAN\n", package_name))
        install.packages(package_name)
      }
      
      # Load the package
      library(package_name, character.only = TRUE)
    } else {
      cat(sprintf("Package already installed: %s\n", package_name))
    }
  }
}

# List of packages with their installation sources
packages <- c(
  spdep = "CRAN", 
  spatialEco = "CRAN", 
  sp = "CRAN", 
  RANN = "CRAN", 
  igraph = "CRAN", 
  RSpectra = "CRAN", 
  Matrix = "CRAN", 
  patchwork = "CRAN", 
  ggpubr = "CRAN", 
  egg = "CRAN",
  ImpulseDE2 = "YosefLab/ImpulseDE2", 
  kneedle = "etam4260/kneedle", 
  gasper = "fabnavarro/gasper"
)

# Run the installation and loading function
install_and_load(packages)

# 
library(spdep)
library(spatialEco)
library(sp)
library(RANN)
library(igraph)
library(RSpectra)
library(kneedle)
library(gasper)
library(Matrix)
library(patchwork)
library(ImpulseDE2)
library(ggpubr)
library(egg)
# calculate laplacian
cal_laplacian <- function(W, type = "unnormalized") {
  if (is(W, 'sparseMatrix')) {
    D <- Diagonal(nrow(W), rowSums(W))
    D_half_inv <- Diagonal(nrow(W), 1/sqrt(rowSums(W)))
    D_inv <- Diagonal(nrow(W), 1/rowSums(W))
  } else {
    D <- diag(rowSums(W))
    D_half_inv <- diag(1/sqrt(rowSums(W)))
    D_inv <- diag(1/rowSums(W))
  }
  if (type == "unnormalized") {
    L <- D - W
  } else if (type == "normalized") {
    L <- diag(nrow(W)) - D_half_inv %*% W %*% D_half_inv
  } else if (type == "randomwalk") {
    L <- diag(nrow(W)) - D_inv %*% W
  } else {
    stop("Invalid type provided.
         Choose from 'unnormalized',
         'normalized', or 'randomwalk'.")
  }
  return(L)
}

# fast eigen calculation
FastDecompositionLap <- function(laplacianMat = NULL, k_fold = 1.5,which = "LM", sigma = NULL, opts = list(),
                                 lower = TRUE, ...) {
  res_decom <- eigs_sym(laplacianMat, k = k_fold* sqrt(ncol(laplacianMat)),which = which, sigma = sigma, opts = opts,
                        lower = lower)
  return(list(evalues = rev(res_decom$values),
              evectors = res_decom$vectors[, ncol(res_decom$vectors):1]))
}

# Graph Fourier Transform
gft <- function(signal, U) {
  # Convert the signal to a matrix if it is a dense vector
  if (is.vector(signal)) {
    signal <- matrix(signal, ncol = 1)  # Convert vector to a column matrix
  }
  
  # Ensure that U is a compatible dense matrix or sparse matrix
  U <- as(U, "CsparseMatrix")  # Convert U to sparse matrix format if necessary
  
  # Compute the Graph Fourier Transform
  result <- t(U) %*% signal
  
  return(result)
}


# plot FM
plot_FM <- function(input = NULL,FM_idx =c(1:20),ncol = 5){
  
   # Create a vector of frequency column names
  freq_columns <- paste0("Freq", FM_idx)
  # Create a list to store the plots
  plot_list <- list()
  
  # Loop through each frequency column and generate plots
  for (freq in freq_columns) {
    # Generate the scatter plot
    p <- ggplot(df_hex_combine, aes_string(x = "x", y = "y")) +
      geom_point(aes_string(color = freq), size =1) +
      # labs(title = paste("L-",freq),
      #      x = "X",
      #      y = "Y",
      #      color = freq) 
      guides(color = "none")+
      theme_void()+scale_color_viridis_c(option = "magma")
    
    # Add the plot to the list
    plot_list[[freq]] <- p
  }
  
  # Combine all plots into a single layout
  combined_plot <- do.call(ggarrange, c(plot_list, list(ncol = ncol)))
  # Print the combined plot
  return(combined_plot)
}

# calculate Eign value
Cal_Eigen <- function(data.in = NULL, k=25, k_fold = 15,sensitivity = 2){
  data.in <- data.in
  k_fold <- k_fold # how many eigen values and vectors need from low frequency
  sensitivity <- sensitivity
  k <- k  # Number of nearest neighbors
  nn <- nn2(data.in[, c("x", "y")], k = k + 1)  # Include the point itself in neighbors
  
  adj_list <- lapply(1:nrow(data.in), function(i) setdiff(nn$nn.idx[i, ], i))  # Remove self-loops
  # Convert adjacency list to edge list
  edges <- do.call(rbind, lapply(1:length(adj_list), function(i) cbind(i, adj_list[[i]])))
  edges <- unique(t(apply(edges, 1, sort)))  # Remove duplicate edges
  # Create the graph
  g <- graph_from_edgelist(edges, directed = FALSE)
  # Ensure vertex names are correctly set
  V(g)$name <- 1:vcount(g)
  
  A <- as_adjacency_matrix(g, sparse = TRUE)
  #
  # Create the graph
  # Create a Graph object
  L <- cal_laplacian(A,"normalized") 
  #
  L_decompose_res <-  FastDecompositionLap(L, k_fold = k_fold,which = "SM")
  # calculate eigen value and vector
  eigenvalue <- L_decompose_res$evalues
  eigenvector <- L_decompose_res$evectors
  # cut low frequency
  knee <- kneedle(x = seq_along(eigenvalue),y = eigenvalue,sensitivity = sensitivity)
  print(paste0("cutoff of low-frequency FM: ",knee[1]))
  plot(eigenvalue, type = "o", col = "blue", xlab = "Index", ylab = "Value", main = "Line Plot of a Numeric Vector")
  abline(v = knee[1],col = "red",lwd =2, lty = 2)
  return(list(knee, eigenvector,eigenvalue))
}

Cal_GCC <- function(data.in = NULL, knee = NULL,signal1 = NULL, signal2 = NULL, eigenvector = NULL){
  data.in <- as.data.frame(data.in)
  knee <- knee
  eigenvector <- eigenvector
  #
  # evectors_named <- eigenvector
  #colnames(eigenvector) <- paste0("Freq",1:ncol(eigenvector))
  #df_hex_combine <- cbind(data.in,eigenvector)
  # plot_FM(df_hex_combine,FM_idx = 1:49,ncol = 7)
  # Apply GFT to both signals
  gft_signal1 <- gft(as.numeric(data.in[,signal1]), eigenvector)
  gft_signal2 <- gft(as.numeric(data.in[,signal2]), eigenvector)
  consin.similiary <- as.numeric(lsa::cosine((gft_signal1[2:knee[1]]),(gft_signal2[2:knee[1]])))
  return(consin.similiary)
}

