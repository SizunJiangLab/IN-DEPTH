## set working directory
loadingfunction_wdpath <- c("./src/")
source(file = file.path(loadingfunction_wdpath,"SGCC_code","1-SpaGFT_preload_function.R"))


# plot eigen value and vector
x_border <- 60
# Create a base scatter plot
base_data <- expand.grid(x = 1:x_border, y = 1:x_border)
k_range <- seq(from = 50, to =600, by =50)
eigne.value.list <- c()
for (i in 1:length(k_range)){
  tmp.name <- k_range[i]
  output.res <- Cal_Eigen(data.in = base_data,k = k_range[i],k_fold = 10)
  eigne.value.list.tmp <- data.frame(eignvalue = output.res[[3]], name = paste0("k_",tmp.name),order = 1:length(output.res[[3]]))
  eigne.value.list <- rbind.data.frame(eigne.value.list, eigne.value.list.tmp)
  print(i)
}
p.eigenvalue <- ggplot(eigne.value.list,aes(x = order , y = eignvalue , color = name ))+
  geom_line() + theme_classic()
ggsave(plot = p.eigenvalue,filename = paste0("Supplementary Figure 3A.svg"),
       device = "svg",width = 12,height = 12)


#-----------------------------------------------
# revision-calculate the Pairwise Correlations
library(dplyr)
library(tidyr)

eigne.value.list <- eigne.value.list %>%
  mutate(k_num = as.numeric(sub("k_", "", name)))

eigvals_wide <- eigne.value.list %>%
  arrange(k_num, order) %>%
  pivot_wider(
    id_cols = name,
    names_from = order,
    values_from = eignvalue
  )

# Save the eigenvalues for further analysis
saveRDS(eigvals_wide, "/Users/rongting/Documents/ResearchProject/IN-DEPTH/data/eigvals_wide.rds")
# write.csv(eigvals_wide, "/Users/rongting/Documents/ResearchProject/IN-DEPTH/data/eigvals_wide.csv", row.names = FALSE)
# Load later
# eigvals_wide <- readRDS("eigvals_wide.rds")
# eigvals_wide <- read.csv("eigvals_wide.csv")

# source("03_SuppFig3_revision.R")
#-----------------------------------------------


# 
output.res <- Cal_Eigen(data.in = base_data,k = 50,k_fold = 10)
eigenvector <- output.res[[2]]
colnames(eigenvector) <- paste0("Freq",1:ncol(eigenvector))
df_hex_combine <- data.frame(x= base_data$x, y = base_data$y, eigenvector)
p.eigvector <- plot_FM(input = df_hex_combine,FM_idx = 2:71,ncol = 10)
p.eigvector
ggsave(plot = p.eigvector,filename = paste0("Supplementary Figure 3B left.tiff"),
       device = "tiff",dpi = 300,width = 20,height = 15)  

#
output.res <- Cal_Eigen(data.in = base_data,k = 400,k_fold = 10)
eigenvector <- output.res[[2]]
colnames(eigenvector) <- paste0("Freq",1:ncol(eigenvector))
df_hex_combine <- data.frame(x= base_data$x, y = base_data$y, eigenvector)
p.eigvector <- plot_FM(input = df_hex_combine,FM_idx = 2:71,ncol = 10)
p.eigvector
ggsave(plot = p.eigvector,filename = paste0("Supplementary Figure 3B right.tiff"),
       device = "tiff",dpi = 300,width = 20,height = 15)  