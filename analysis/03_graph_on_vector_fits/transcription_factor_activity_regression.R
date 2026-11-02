

#-----------------------------------
# LUAD TF
#-----------------------------------

rm(list=ls())
set.seed(2020)

pacman::p_load(spdep,fpc,GWmodel,sp,ppcor,spatstat,reshape2, dplyr,
               tmap,ggplot2,tidyverse, Matrix, igraph, ggraph, matrixcalc,
               mvtnorm,MCMCpack)


finalPATIDs <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\luad_HR_PATIDs.rds")


source("C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/H&E Lung TCGA Graph Dataset/LSGPR/lsGPR_VS.R")

#LUAD 
activity<- readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\top_tf_df_highrisk_luad.rds')
grouped <- activity %>% group_by(source) %>%
  group_split()
length(grouped)
for (j in 1:length(grouped)) {
  # j=12
  act <- grouped[[j]]$source[1]
  print(paste("Group", j, ":",act )) 
  
  currentPATIDs<-grouped[[j]]$condition 
  table(finalPATIDs %in% currentPATIDs)
  #315
  
  Y <- as.matrix(grouped[[j]]$score, ncol=1)
  rownames(Y) <- currentPATIDs
  
  Y <- Y[match(finalPATIDs, rownames(Y)), , drop=FALSE]
  
  finalmya <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\luad_HR_mya.rds")
  finalmyU <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\luad_HR_myU.rds")
  
  test_indices <- 1:10  # Adjust as needed
  a.te <- finalmya[test_indices]
  Y.te <- Y[test_indices]
  U.te <- finalmyU[test_indices]
  
  a <- finalmya[-test_indices]
  U <- finalmyU[-test_indices]
  Y <- Y[-test_indices]
  
  result <- lsGPR_VS(a,a.te,U,U.te,z=NULL,z.te=NULL,Y,burn=500,nrun=1000)
  burn<- result$burn
  nrun <- result$nrun
  
  GAMMA <- result$GAMMA
  # After MCMC
  post_inclusion = rowMeans(GAMMA[, (burn+1):nrun]);post_inclusion
  
  #important_nodes = which(post_inclusion > 0.5);important_nodes
  
  # These nodes have network topology features that predict pathway activity
  
  BETA <- result$BETA
  post_beta = rowMeans(BETA[, (burn+1):nrun]);post_beta
  
  saveRDS(result, paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\luad\\HR_luad_TF_lsGPRVS_result1000_",act,"_activity.rds"))
  
}

#---------------------------------------------
# LUAD LR TF
#---------------------------------------------

rm(list=ls())
set.seed(2020)

pacman::p_load(spdep,fpc,GWmodel,sp,ppcor,spatstat,reshape2, dplyr,
               tmap,ggplot2,tidyverse, Matrix, igraph, ggraph, matrixcalc,
               mvtnorm,MCMCpack)


finalPATIDs <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\luad_LR_PATIDs.rds")


source("C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/H&E Lung TCGA Graph Dataset/LSGPR/lsGPR_VS.R")

#LUAD 
activity<- readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\top_tf_df_lowrisk_luad.rds')
grouped <- activity %>% group_by(source) %>%
  group_split()
length(grouped)
for (j in 1:length(grouped)) {
  # j=12
  act <- grouped[[j]]$source[1]
  print(paste("Group", j, ":",act )) 
  
  currentPATIDs<-grouped[[j]]$condition 
  table(finalPATIDs %in% currentPATIDs)
  #315
  
  Y <- as.matrix(grouped[[j]]$score, ncol=1)
  rownames(Y) <- currentPATIDs
  
  Y <- Y[match(finalPATIDs, rownames(Y)), , drop=FALSE]
  
  finalmya <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\luad_LR_mya.rds")
  finalmyU <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\luad_LR_myU.rds")
  
  test_indices <- 1:10  # Adjust as needed
  a.te <- finalmya[test_indices]
  Y.te <- Y[test_indices]
  U.te <- finalmyU[test_indices]
  
  a <- finalmya[-test_indices]
  U <- finalmyU[-test_indices]
  Y <- Y[-test_indices]
  
  result <- lsGPR_VS(a,a.te,U,U.te,z=NULL,z.te=NULL,Y,burn=500,nrun=1000)
  burn<- result$burn
  nrun <- result$nrun
  
  GAMMA <- result$GAMMA
  # After MCMC
  post_inclusion = rowMeans(GAMMA[, (burn+1):nrun]);post_inclusion
  
  #important_nodes = which(post_inclusion > 0.5);important_nodes
  
  # These nodes have network topology features that predict pathway activity
  
  BETA <- result$BETA
  post_beta = rowMeans(BETA[, (burn+1):nrun]);post_beta
  
  saveRDS(result, paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\luad\\LR_luad_TF_lsGPRVS_result1000_",act,"_activity.rds"))
  
}



#-------------------------------------------------
# lsgpr_VS_plot_results
#-------------------------------------------------

rm(list=ls())
set.seed(2020)

pacman::p_load(spdep,fpc,GWmodel,sp,ppcor,spatstat,reshape2, dplyr,
               tmap,ggplot2,tidyverse, Matrix, igraph, ggraph, matrixcalc,
               mvtnorm,MCMCpack)
#----------------------------------------
# LUAD HR
#----------------------------------------

activity<- readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\top_tf_df_highrisk_luad.rds')
grouped <- activity %>% group_by(source) %>%
  group_split()
length(grouped)


beta<- gamma_inc <- matrix(0, nrow=20, ncol=length(grouped))
for (j in 1:length(grouped)) {
  # j=1
  act <- grouped[[j]]$source[1]
  print(paste("Group", j, ":",act )) 
  
  result <- readRDS(paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\luad\\HR_luad_TF_lsGPRVS_result1000_",act,"_activity.rds"))
  
  burn<- result$burn
  nrun <- result$nrun
  
  GAMMA <- result$GAMMA
  # After MCMC
  post_inclusion = rowMeans(GAMMA[, (burn+1):nrun]);post_inclusion
  post_inclusion = rowMeans(GAMMA);post_inclusion
  
  #important_nodes = which(post_inclusion > 0.5);important_nodes
  
  # These nodes have network topology features that predict pathway activity
  
  BETA <- result$BETA
  post_beta = rowMeans(BETA[, (burn+1):nrun]);post_beta
  post_beta = rowMeans(BETA);post_beta
  # Rank nodes by importance
  ranked_nodes = order(post_beta, decreasing=TRUE);ranked_nodes
  
  # PRED <- result$post.pred     # Predicted pathway activity for test patients
  # # Posterior mean prediction
  # pred_mean = rowMeans(PRED[, (burn+1):nrun]);pred_mean
  # 
  # # 95% credible interval
  # pred_lower = apply(PRED[, (burn+1):nrun], 1, quantile, 0.025);pred_lower
  # pred_upper = apply(PRED[, (burn+1):nrun], 1, quantile, 0.975);pred_upper
  
  #store
  gamma_inc[,j] <- post_inclusion
  beta[,j] <- post_beta
  
}
colnames(gamma_inc) <- colnames(beta) <- sapply(1:length(grouped), function(j) grouped[[j]]$source[1])
read_path<- "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\adj_matricesV3\\"

adj <- as.matrix(read.csv(paste0(read_path, 'TCGA-62-A46R-01Z-00-DX1.csv'))[,-1]) # Remove the first column
adj
rownames(beta) = rownames(gamma_inc) <- colnames(adj)

# Heatmap plotting
library(ggplot2)
library(reshape2)  # for melt() function
library(tidyverse)

# Create a sample matrix (replace this with your matrix)
set.seed(123)


#--------------------------Beta------------------------------------------ 

# Set your threshold value
threshold <- quantile(beta, probs = c(0.25, 0.5, 0.75))[2]
threshold <- round(threshold,5)

# Convert matrix to long format for ggplot
mat_melt <- melt(beta)
colnames(mat_melt) <- c("HPC", "Pathways", "Betas")

# Create heatmap with threshold
beta<- ggplot(mat_melt, aes(x = Pathways, y = HPC, fill = Betas)) +
  geom_tile() +
  # Color scale with threshold - diverging at the threshold value
  scale_fill_gradient2(
    low = "red",
    mid = "white", 
    high = "blue",
    midpoint = threshold,
    limits = c(min(mat_melt$Betas), max(mat_melt$Betas)),
    name = "Beta"
  ) +
  # Add text labels showing values
  #geom_text(aes(label = round(Value, 1)), size = 3) +
  # Theme customization
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Beta's for LUAD high risk",
    #subtitle = paste('median Beta ', threshold),
    x = "TF Activities",
    y = "HPCs"
  ) +
  coord_fixed()  # Keep aspect ratio square

beta

#--------------------------PIP------------------------------------------ 

# Set your threshold value
threshold <- quantile(gamma_inc, probs = c(0.25, 0.5, 0.75))[2]
threshold <- round(threshold,5)

threshold1 <- 0.5

mat_melt1 <- melt(gamma_inc)
colnames(mat_melt1) <- c("HPC", "Pathways", "PIP")

# Create heatmap with threshold
pip1<- ggplot(mat_melt1, aes(x = Pathways, y = HPC, fill = PIP)) +
  geom_tile() +
  # Color scale with threshold - diverging at the threshold value
  scale_fill_gradient2(
    low = "blue",
    mid = "white", 
    high = "red",
    midpoint = threshold,
    limits = c(min(mat_melt1$PIP), max(mat_melt1$PIP)),
    name = "PIP"
  ) +
  # Add asterisk where PIP > threshold
  geom_text(
    data = subset(mat_melt1, PIP > threshold1),
    aes(label = "*"),
    color = "black",
    size = 5
  ) +
  # Add text labels showing values
  #geom_text(aes(label = round(Value, 1)), size = 3) +
  # Theme customization
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Posterior Inclusion Probability for LUAD high risk",
    #subtitle = paste('median PIP ', threshold),
    x = "TF Activities",
    y = "HPCs"
  ) +
  coord_fixed()  # Keep aspect ratio square

pip1

# Save the plots

library(patchwork)

new<- beta + pip1 + plot_layout(ncol = 2)
new
ggsave("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\Plots\\luad_HR_TF_pip_heatmap.png", 
       plot = new, width = 10, height = 7, dpi = 300,bg = "white")


#----------------------------------------
# LUAD LR
#----------------------------------------

#LUAD or LUSC
activity<- readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\top_tf_df_lowrisk_luad.rds')
grouped <- activity %>% group_by(source) %>%
  group_split()
length(grouped)


beta<- gamma_inc <- matrix(0, nrow=20, ncol=length(grouped))
for (j in 1:length(grouped)) {
  # j=1
  act <- grouped[[j]]$source[1]
  print(paste("Group", j, ":",act )) 
  
  result <- readRDS(paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\luad\\LR_luad_TF_lsGPRVS_result1000_",act,"_activity.rds"))
  
  burn<- result$burn
  nrun <- result$nrun
  
  GAMMA <- result$GAMMA
  # After MCMC
  post_inclusion = rowMeans(GAMMA[, (burn+1):nrun]);post_inclusion
  post_inclusion = rowMeans(GAMMA);post_inclusion
  #important_nodes = which(post_inclusion > 0.5);important_nodes
  
  # These nodes have network topology features that predict pathway activity
  
  BETA <- result$BETA
  post_beta = rowMeans(BETA[, (burn+1):nrun]);post_beta
  post_beta = rowMeans(BETA);post_beta
  
  # Rank nodes by importance
  ranked_nodes = order(post_beta, decreasing=TRUE);ranked_nodes
  
  # PRED <- result$post.pred     # Predicted pathway activity for test patients
  # # Posterior mean prediction
  # pred_mean = rowMeans(PRED[, (burn+1):nrun]);pred_mean
  # 
  # # 95% credible interval
  # pred_lower = apply(PRED[, (burn+1):nrun], 1, quantile, 0.025);pred_lower
  # pred_upper = apply(PRED[, (burn+1):nrun], 1, quantile, 0.975);pred_upper
  
  #store
  gamma_inc[,j] <- post_inclusion
  beta[,j] <- post_beta
  
}
colnames(gamma_inc) <- colnames(beta) <- sapply(1:length(grouped), function(j) grouped[[j]]$source[1])
read_path<- "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\adj_matricesV3\\"

adj <- as.matrix(read.csv(paste0(read_path, 'TCGA-62-A46R-01Z-00-DX1.csv'))[,-1]) # Remove the first column
adj
rownames(beta) = rownames(gamma_inc) <- colnames(adj)

# Heatmap plotting
library(ggplot2)
library(reshape2)  # for melt() function
library(tidyverse)

# Create a sample matrix (replace this with your matrix)
set.seed(123)


#--------------------------Beta------------------------------------------ 

# Set your threshold value
threshold <- quantile(beta, probs = c(0.25, 0.5, 0.75))[2]
threshold <- round(threshold,5)

# Convert matrix to long format for ggplot
mat_melt <- melt(beta)
colnames(mat_melt) <- c("HPC", "Pathways", "Betas")

# Create heatmap with threshold
beta<- ggplot(mat_melt, aes(x = Pathways, y = HPC, fill = Betas)) +
  geom_tile() +
  # Color scale with threshold - diverging at the threshold value
  scale_fill_gradient2(
    low = "red",
    mid = "white", 
    high = "blue",
    midpoint = threshold,
    limits = c(min(mat_melt$Betas), max(mat_melt$Betas)),
    name = "Beta"
  ) +
  # Add text labels showing values
  #geom_text(aes(label = round(Value, 1)), size = 3) +
  # Theme customization
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Beta's for LUAD low risk",
    #subtitle = paste('median Beta ', threshold),
    x = "TF Activities",
    y = "HPCs"
  ) +
  coord_fixed()  # Keep aspect ratio square

beta

#--------------------------PIP------------------------------------------ 

# Set your threshold value
threshold <- quantile(gamma_inc, probs = c(0.25, 0.5, 0.75))[2]
threshold <- round(threshold,5)

threshold1<-0.5;

mat_melt2 <- melt(gamma_inc)
colnames(mat_melt2) <- c("HPC", "Pathways", "PIP")

# Create heatmap with threshold
pip2<- ggplot(mat_melt2, aes(x = Pathways, y = HPC, fill = PIP)) +
  geom_tile() +
  # Color scale with threshold - diverging at the threshold value
  scale_fill_gradient2(
    low = "blue",
    mid = "white", 
    high = "red",
    midpoint = threshold,
    limits = c(min(mat_melt2$PIP), max(mat_melt2$PIP)),
    name = "PIP"
  ) +
  # Add asterisk where PIP > threshold
  geom_text(
    data = subset(mat_melt2, PIP > threshold1),
    aes(label = "*"),
    color = "black",
    size = 5
  ) +
  # Add text labels showing values
  #geom_text(aes(label = round(Value, 1)), size = 3) +
  # Theme customization
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Posterior Inclusion Probability for LUAD low risk",
    #subtitle = paste('median PIP ', threshold),
    x = "TF Activities",
    y = "HPCs"
  ) +
  coord_fixed()  # Keep aspect ratio square

pip2

# Save the plots

library(patchwork)

new <- beta + pip2 + plot_layout(ncol = 2)
new

ggsave("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\Plots\\luad_LR_tf_pip_heatmap.png", 
       plot = new, width = 10, height = 7, dpi = 300,bg = "white")


new <- pip1 + pip2 + plot_layout(ncol = 2)
new

ggsave("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\Plots\\luad_HR-LR_TF_pipmap.png", 
       plot = new, width = 10, height = 7, dpi = 300,bg = "white")






# Or combine into long format
combined <- rbind(
  mat_melt1 %>% as.data.frame() %>% mutate(risk = "high"),
  mat_melt2 %>% as.data.frame() %>% mutate(risk = "low")
)

write.csv(combined, "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\luad\\luad_tf_lsGPRVS_pip_longformat.csv")

