#--------------------------------------------------------------
# Read files
#--------------------------------------------------------------
rm(list=ls())
set.seed(2020)

pacman::p_load(spdep,fpc,GWmodel,sp,ppcor,spatstat,reshape2, dplyr,
               tmap,ggplot2,tidyverse, Matrix, igraph, ggraph, matrixcalc)

#-----------------------------------------------------------------------------
# LSEM
source("C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/H&E Lung TCGA Graph Dataset/LSGPR/lsEM.R")
# constrained_adjacency_covariance_function
source("C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/H&E Lung TCGA Graph Dataset/codes/constrained_adjacency_covariance_function.R")
#-----------------------------------------------------------------------------

# number of latent channels
myd=5

# read_path<- "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\H&E Lung TCGA Graph Dataset\\complete_graphs\\adj_matrices\\"
read_path<- "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\adj_matricesV3\\"

setwd(read_path)
EllFiles<-Sys.glob("*")
num_files<-length(EllFiles)
length(unique(EllFiles))

log_cts<- readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\log_cts_high_risk_luad.rds')
#log_cts<- readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\log_cts_luad.rds')

log_cnts_PATID<-rownames(log_cts)
#log_cnts_PATID <- gsub("^(.{12}).*", "\\1", log_cnts_PATID) # Extract the first 12 characters


library(magrittr)  # or dplyr
PATID <- EllFiles %>%
  gsub("\\.csv$", "", .)
  # gsub("-", ".", .) %>%
  # gsub("^(.{12}).*", "\\1", .) # Extract the first 12 characters

table(PATID %in% log_cnts_PATID) # check if the PATID in the files are present in the log counts
# luAD - high risk
# FALSE  TRUE 
# 204   198 

# # Create PATID from file names (first 12 characters)
# PATID <- EllFiles %>%
#   gsub("\\.gml$", "", .) %>%
#   gsub("-", ".", .) %>%
#   gsub("^(.{12}).*", "\\1", .)

# Create logical index for matching files
matching_files_idx <- PATID %in% log_cnts_PATID
cat("Total files:", length(EllFiles), "\n")
cat("Matching files:", sum(matching_files_idx), "\n")

# Filter to only matching files
EllFiles_filtered <- EllFiles[matching_files_idx]
PATID_filtered <- PATID[matching_files_idx]

# Initialize storage for filtered files only
num_filtered_files <- length(EllFiles_filtered)
adj_matrix <- list()
myU <- list()
mya <- vector()
PATIDs <- vector()

#for 3
exclude_indices_luad_HR <- c(101,102,114,167)
#for 4
exclude_indices_luad_HR <- c(2, 101,102,114, 115,167)
#for 5
exclude_indices_luad_HR <- c(2, 101,102,114, 115,167,198)

exclude_indices <- exclude_indices_luad_HR # Change as needed
valid_indices <- setdiff(1:num_filtered_files, exclude_indices)

# 
# exclude_indices_luad_LR <- c(246)
# exclude_indices <- exclude_indices_luad_LR # Change as needed
# valid_indices <- setdiff(1:num_filtered_files, exclude_indices)


# Loop through filtered files only
for (i in valid_indices) {
  # i=1
  print(paste("Processing file", i, "of", num_filtered_files))
  files <- EllFiles_filtered[i]
  
   PATIDs[i] <- files %>%
    gsub("\\.csv$", "", .) 
   
  adj <- as.matrix(read.csv(paste0(read_path, files))[,-1]) # Remove the first column

  # Identify and replace
  zero_rows <- apply(adj == 0, 1, all)
  zero_cols <- apply(adj == 0, 2, all)
  
  adj_sp <- adj
  adj_sp[, zero_cols] <- NA
  adj_sp[zero_rows, ] <- NA
  #diag(adj) <- 1
  adj_sp
  
  library(softImpute)
  fit <- softImpute::softImpute(adj_sp, rank.max = nrow(adj)-1, lambda = 0.1, type = "svd") 
  A_completed <- softImpute::complete(adj_sp,fit)
  A_imputed <- adj
  A_imputed[adj == 0] <- A_completed[adj == 0]
  diag(A_imputed) <- 1
  A_imputed[A_imputed < 0] <- 0  # Non-negative
  # preserve edges that exist in only one direction
  adj <- ifelse((A_imputed > 0) & (t(A_imputed) > 0), sqrt(A_imputed * t(A_imputed)), pmax(A_imputed, t(A_imputed)))
  
  
  adj_matrix[[i]] <- adj
  
  
  # library(eimpute)
  # x_impute <- eimpute(adj, r=nrow(adj)-1)
  # adj <- x_impute$x.imp

  # # # Convert directed to undirected using geometric mean
  # directed_to_undirected_geomean <- function(adj_matrix) {
  #   sqrt(adj_matrix * t(adj_matrix))
  # }
  # adj <- directed_to_undirected_geomean(adj)
  # 
  
  result<-EM_toU(adj, d=myd, u1=0.5, sig.a=2, sig.u=2, tol=10^(-10), maxiter=3)
  
  mya[i]<-result$a
  myU[[i]] <-result$U
  
  #saveRDS(result, paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\lusc\\Ua_",PATIDs[i],".rds"))
}

finalPATIDs <- PATIDs[valid_indices]
finalmya <- mya[valid_indices]
finalmyU <- myU[valid_indices]

saveRDS(finalPATIDs, paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\luad_HR_", myd,"_PATIDs.rds"))
saveRDS(finalmya, paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\luad_HR_d_", myd,"_mya.rds"))
saveRDS(finalmyU, paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\luad_HR_d_", myd,"_myU.rds"))


#----------------------------------------------------
# LUAD LR
#----------------------------------------------------

rm(list=ls())
set.seed(2020)

pacman::p_load(spdep,fpc,GWmodel,sp,ppcor,spatstat,reshape2, dplyr,
               tmap,ggplot2,tidyverse, Matrix, igraph, ggraph, matrixcalc)

#-----------------------------------------------------------------------------
# LSEM
source("C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/H&E Lung TCGA Graph Dataset/LSGPR/lsEM.R")
# constrained_adjacency_covariance_function
source("C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/H&E Lung TCGA Graph Dataset/codes/constrained_adjacency_covariance_function.R")
#-----------------------------------------------------------------------------

myd=5

# read_path<- "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\H&E Lung TCGA Graph Dataset\\complete_graphs\\adj_matrices\\"
read_path<- "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\adj_matricesV3\\"

setwd(read_path)
EllFiles<-Sys.glob("*")
num_files<-length(EllFiles)
length(unique(EllFiles))

log_cts<- readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\log_cts_low_risk_luad.rds')

log_cnts_PATID<-rownames(log_cts)
#log_cnts_PATID <- gsub("^(.{12}).*", "\\1", log_cnts_PATID) # Extract the first 12 characters


library(magrittr)  # or dplyr
PATID <- EllFiles %>%
  gsub("\\.csv$", "", .)
# gsub("-", ".", .) %>%
# gsub("^(.{12}).*", "\\1", .) # Extract the first 12 characters

table(PATID %in% log_cnts_PATID) # check if the PATID in the files are present in the log counts
# luAD - high risk
# FALSE  TRUE 
# 201     201

# # Create PATID from file names (first 12 characters)
# PATID <- EllFiles %>%
#   gsub("\\.gml$", "", .) %>%
#   gsub("-", ".", .) %>%
#   gsub("^(.{12}).*", "\\1", .)

# Create logical index for matching files
matching_files_idx <- PATID %in% log_cnts_PATID
cat("Total files:", length(EllFiles), "\n")
cat("Matching files:", sum(matching_files_idx), "\n")

# Filter to only matching files
EllFiles_filtered <- EllFiles[matching_files_idx]
PATID_filtered <- PATID[matching_files_idx]

# Initialize storage for filtered files only
num_filtered_files <- length(EllFiles_filtered)
adj_matrix <- list()
myU <- list()
mya <- vector()
PATIDs <- vector()


#for 3
# exclude_indices_luad_LR <- c(11,23,61,77,82,101,112,121,138)
# for 4
# exclude_indices_luad_LR <- c(11,23,61,77,82,101,108,112,121,138, 167)
#for 5
# exclude_indices_luad_LR <- c(11,23,61,77,82,101,108, 112,121,138,167)


exclude_indices <- exclude_indices_luad_LR # Change as needed
valid_indices <- setdiff(1:num_filtered_files, exclude_indices)


# Loop through filtered files only
for (i in valid_indices) {
  # i=1
  print(paste("Processing file", i, "of", num_filtered_files))
  files <- EllFiles_filtered[i]
  
  PATIDs[i] <- files %>%
    gsub("\\.csv$", "", .) 
  
  adj <- as.matrix(read.csv(paste0(read_path, files))[,-1]) # Remove the first column
  
  # Identify and replace
  zero_rows <- apply(adj == 0, 1, all)
  zero_cols <- apply(adj == 0, 2, all)
  
  adj_sp <- adj
  adj_sp[, zero_cols] <- NA
  adj_sp[zero_rows, ] <- NA
  #diag(adj) <- 1
  adj_sp
  
  library(softImpute)
  fit <- softImpute::softImpute(adj_sp, rank.max = nrow(adj)-1, lambda = 0.1, type = "svd") 
  A_completed <- softImpute::complete(adj_sp,fit)
  A_imputed <- adj
  A_imputed[adj == 0] <- A_completed[adj == 0]
  diag(A_imputed) <- 1
  A_imputed[A_imputed < 0] <- 0  # Non-negative
  # preserve edges that exist in only one direction
  adj <- ifelse((A_imputed > 0) & (t(A_imputed) > 0), sqrt(A_imputed * t(A_imputed)), pmax(A_imputed, t(A_imputed)))
  
  
  adj_matrix[[i]] <- adj
  
  
  # library(eimpute)
  # x_impute <- eimpute(adj, r=nrow(adj)-1)
  # adj <- x_impute$x.imp
  
  # # # Convert directed to undirected using geometric mean
  # directed_to_undirected_geomean <- function(adj_matrix) {
  #   sqrt(adj_matrix * t(adj_matrix))
  # }
  # adj <- directed_to_undirected_geomean(adj)
  # 
  
  result<-EM_toU(adj, d=myd, u1=0.5, sig.a=2, sig.u=2, tol=10^(-10), maxiter=3)
  
  mya[i]<-result$a
  myU[[i]] <-result$U
  
  #saveRDS(result, paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\lusc\\Ua_",PATIDs[i],".rds"))
}

finalPATIDs <- PATIDs[valid_indices]
finalmya <- mya[valid_indices]
finalmyU <- myU[valid_indices]


saveRDS(finalPATIDs, paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\luad_LR_", myd,"_PATIDs.rds"))
saveRDS(finalmya, paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\luad_LR_d_",myd,"_mya.rds"))
saveRDS(finalmyU, paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\luad_LR_d_",myd,"_myU.rds"))
# adj_matrix
# table(hpc_names)
# table(lengths)

# saveRDS(adj_matrix, "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\lusc_adj_matrixs.rds")
# #saveRDS(adj_matrix, "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\luad_adj_matrixs.rds")
# 
# 
# patdata<-as.data.frame(cbind(PATIDs, lengths, hpc_names))
# saveRDS(patdata, "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\PATID_lengths_hpc_names.rds")
# 
# patdata %>% group_by(lengths,hpc_names) -> patdata2
# saveRDS(patdata2, "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\PATID_lengths_hpc_names_grouped.rds")
# 


# #---------------------------------------------------------------------------------------------------------
#  #  ls.mcmc2pred.gpa(a,a.te,U,U.te,Y)
# #---------------------------------------------------------------------------------------------------------
# rm(list=ls())
# set.seed(2020)
# 
# pacman::p_load(spdep,fpc,GWmodel,sp,ppcor,spatstat,reshape2, dplyr,
#                tmap,ggplot2,tidyverse, Matrix, igraph, ggraph, matrixcalc)
# 
# 
# #finalPATIDs <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\lusc\\lusc_PATIDs.rds")
# finalPATIDs <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\luad\\luad_PATIDs.rds")
# 
# 
# # genename<- 'DPM1' #Adjust as needed
# # Y <- log_cts[genename, match(finalPATIDs, log_cnts_PATID)]
# 
# source("C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/H&E Lung TCGA Graph Dataset/LSGPR/lsGPR.R")
# #LUAD or LUSC
# #activity<- readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\activity_scores_lusc.rds')
# activity<- readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\activity_scores_luad.rds')
# grouped <- activity %>% group_by(source) %>%
#   group_split()
# length(grouped)
# for (j in 1:length(grouped)) {
#   act <- grouped[[j]]$source[1]
#   print(paste("Group", j, ":",act )) 
#   
#   currentPATIDs<-grouped[[j]]$condition %>%  gsub("^(.{12}).*", "\\1", .) # Extract the first 12 characters
#   table(finalPATIDs %in% currentPATIDs)
#   #315
#   
#   Y <- as.matrix(grouped[[j]]$score, ncol=1)
#   rownames(Y) <- currentPATIDs
#   
#   Y <- Y[match(finalPATIDs, rownames(Y)), , drop=FALSE]
#   
#   finalmya <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\luad\\luad_mya.rds")
#   finalmyU <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\luad\\luad_myU.rds")
#   
#   # finalmya <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\lusc\\lusc_mya.rds")
#   # finalmyU <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\lusc\\lusc_myU.rds")
#   
#   test_indices <- 1:10  # Adjust as needed
#   a.te <- finalmya[test_indices]
#   Y.te <- Y[test_indices]
#   U.te <- finalmyU[test_indices]
#   
#   a <- finalmya[-test_indices]
#   U <- finalmyU[-test_indices]
#   Y <- Y[-test_indices]
#   
#   result <- ls.mcmc2pred.gpa(a, a.te, U, U.te, Y)
#   
#   # saveRDS(result, paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\lusc\\lusc_lsGPR_result_",act,"_activity.rds"))
#   saveRDS(result, paste0("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\luad\\luad_lsGPR_result_",act,"_activity.rds"))
#   
# }
# 
# plot(result$post.loglik, type='l')
# plot(result$PSI1, type='l')
# plot(result$TAU, type='l')



#---------------------------------------------------------------------------------------------------------
# lsGPR_VS(a,a.te,U,U.te,z=NULL,z.te=NULL,Y,init,burn,nrun,prior)
# with Spike n Slab
#---------------------------------------------------------------------------------------------------------


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ComplexHeatmap")
# install.packages("circlize")
# 
# library(ComplexHeatmap)
# library(circlize)
# 
# # Create matrices
# mat1 <- matrix(rnorm(100), nrow = 10)
# mat2 <- matrix(sample(letters[1:3], 100, replace = TRUE), nrow = 10)
# 
# # Create color mapping
# col_fun1 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
# 
# # Create the heatmap with cell annotations
# Heatmap(mat1, 
#         col = col_fun1,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(mat2[i, j], x, y, gp = gpar(fontsize = 8, col = "black"))
#         },
#         name = "Value1")











