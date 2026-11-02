#---------------------------------------------------
# Splitting rnaseq data into two groups based on a condition
#---------------------------------------------------

rm(list=ls())
set.seed(2025)

pacman::p_load(dplyr, tidyverse, ggplot2, cowplot, patchwork, ggridges, ggpubr, MASS)

rnaseq_data <- read.csv("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\raw_counts_luad_418.csv")
rownames(rnaseq_data) <- rnaseq_data[,1]
rnaseq_data <- rnaseq_data[,-1]
rnaseq_data <- as.data.frame(t(rnaseq_data))
rnaseq_data$slide_id <- rownames(rnaseq_data)


metadata<- read.csv("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\luad_slides_metadata.csv")
mapping <- read.csv('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\wsi_rna-seq_mapping_luad_418.csv', row.names = 1)
mapping$slide_id <- rownames(mapping)


#-------------------------------------------------------------------------------
# STRATEGY 1: KEEP FIRST SAMPLE (MOST COMMON)
#-------------------------------------------------------------------------------

# Extract first RNA-seq sample only
mapping$RNA_seq_first <- sapply(strsplit(mapping$RNA.seq_name, ","), 
                                function(x) trimws(x[1]))

# Merge with metadata
merged <- merge(metadata, 
                mapping[, c("slide_id", "RNA_seq_first")],
                by.x = "slides", 
                by.y = "slide_id",
                all.x = TRUE)
merged$RNA_seq_first[1:5]
merged$RNA_seq_first <- gsub("-", ".", merged$RNA_seq_first)
merged$RNA_seq_first[1:5]

# Check how many RNA-seq samples are found in rnaseq_data
table(merged$RNA_seq_first %in% rnaseq_data$slide_id)


# Merge with RNA-seq data

final_data <- merge(merged, 
                    rnaseq_data, 
                    by.x = "RNA_seq_first", 
                    by.y = "slide_id",
                    all.x = TRUE)
# remove NA entries
final_data <- final_data %>% filter(!is.na(RNA_seq_first))

saveRDS(final_data, file = "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\final_data_luad.rds")

str(final_data)

#-------------------------------------------------------------------------------

#subset final data with expression data only with risk groups
# stratify by risk group


high_risk <- final_data %>% filter(risk_group == 1) 
low_risk <- final_data %>% filter(risk_group == 0)

saveRDS(high_risk, file = "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\high_risk_luad.rds")
saveRDS(low_risk, file = "C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\low_risk_luad.rds")

#-------------------------------------------------------------------------------

rm(list=ls())
set.seed(2025)

# read back
pacman::p_load(dplyr, tidyverse, ggplot2, cowplot, patchwork, ggridges, ggpubr, MASS)

final_data <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\final_data_luad.rds")
high_risk <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\high_risk_luad.rds")
low_risk <- readRDS("C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\low_risk_luad.rds")



# Extract expression data only for high and low risk groups
expr_data_high <- high_risk[, c( (12+1) : ncol(final_data))]
rownames(expr_data_high) <- high_risk$slides
# Calculate z-scores
# Quick one-liner to remove zero-variance genes
# log_cts <- expr_data_high %>%
#   select(where(~var(log2(. + 1), na.rm = TRUE) > 0)) %>%
#   mutate(across(everything(), ~log2(. + 1))) %>%
#   scale() %>%
#   as.data.frame()

# Log transform
log_expr <- log2(expr_data_high + 1)

# Calculate variance for each gene (column)
gene_vars <- apply(log_expr, 2, var, na.rm = TRUE)
gene_means <- apply(log_expr, 2, mean, na.rm = TRUE)

# Filter
genes_to_keep <- gene_vars > 0.01 & gene_means > log2(2)

cat("Keeping", sum(genes_to_keep), "genes\n")

# Scale
log_cts_clean <- scale(log_expr[, genes_to_keep]) %>% as.data.frame()

# Verify
cat("Any NaNs?", any(is.nan(as.matrix(log_cts_clean))), "\n")

saveRDS(log_cts_clean, 'C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\log_cts_high_risk_luad.rds')

# Repeat the same process
# for the low risk group

expr_data_low <- low_risk[, c( (12+1) : ncol(final_data))]
rownames(expr_data_low) <- low_risk$slides

# Log transform
log_expr <- log2(expr_data_low + 1)

# Calculate variance for each gene (column)
gene_vars <- apply(log_expr, 2, var, na.rm = TRUE)
gene_means <- apply(log_expr, 2, mean, na.rm = TRUE)

# Filter
genes_to_keep <- gene_vars > 0.01 & gene_means > log2(2)

cat("Keeping", sum(genes_to_keep), "genes\n")

# Scale
log_cts_clean <- scale(log_expr[, genes_to_keep]) %>% as.data.frame()

# Verify
cat("Any NaNs?", any(is.nan(as.matrix(log_cts_clean))), "\n")

saveRDS(log_cts_clean, 'C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\log_cts_low_risk_luad.rds')



#--------------------------------------------------------------------------------------------
# decoupleR
#--------------------------------------------------------------------------------------------
rm(list=ls())
set.seed(2025)

# read back
.libPaths(c('C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/R_packages/',.libPaths()))
lib <- 'C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/R_packages/'


## We load the required packages
library(Seurat)
# install.packages("BiocManager")
# BiocManager::install("decoupleR")
# BiocManager::install("OmnipathR")

library(decoupleR)
library(OmnipathR)
# Only needed for data handling and plotting
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)


net <- decoupleR::get_progeny(organism = 'human', top = 500)
net

# Extract the normalized log-transformed counts
mat <- t(as.matrix(readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\log_cts_high_risk_luad.rds')))
dim(mat)


# Run mlm
acts <- decoupleR::run_mlm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor = 'weight', 
                           minsize = 5)
str(acts)
# tibble [2,772 × 5] (S3: tbl_df/tbl/data.frame)
# $ statistic: chr [1:2772] "mlm" "mlm" "mlm" "mlm" ...
# $ source   : chr [1:2772] "Androgen" "EGFR" "Estrogen" "Hypoxia" ...
# $ condition: chr [1:2772] "TCGA-05-4250-01Z-00-DX1" "TCGA-05-4250-01Z-00-DX1" "TCGA-05-4250-01Z-00-DX1" "TCGA-05-4250-01Z-00-DX1" ...
# $ score    : num [1:2772] 2.045 1.377 0.501 3.176 5.023 ...
# $ p_value  : num [1:2772] 4.08e-02 1.69e-01 6.17e-01 1.49e-03 5.13e-07 ...

table(acts$source)
# Androgen     EGFR Estrogen  Hypoxia JAK-STAT     MAPK     NFkB      p53     PI3K     TGFb     TNFa    Trail 
# 198      198      198      198      198      198      198      198      198      198      198      198 
# VEGF      WNT 
# 198      198 
saveRDS(acts,'C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\LUAD_high_risk_activity_scores.rds')

#------------------------------------------------------------

# Extract the normalized log-transformed counts
mat <- t(as.matrix(readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\log_cts_low_risk_luad.rds')))
dim(mat)


# Run mlm
acts <- decoupleR::run_mlm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor = 'weight', 
                           minsize = 5)
str(acts)
# tibble [2,814 × 5] (S3: tbl_df/tbl/data.frame)
# $ statistic: chr [1:2814] "mlm" "mlm" "mlm" "mlm" ...
# $ source   : chr [1:2814] "Androgen" "EGFR" "Estrogen" "Hypoxia" ...
# $ condition: chr [1:2814] "TCGA-05-4245-01Z-00-DX1" "TCGA-05-4245-01Z-00-DX1" "TCGA-05-4245-01Z-00-DX1" "TCGA-05-4245-01Z-00-DX1" ...
# $ score    : num [1:2814] -3.727 2.3 2.595 -0.825 -7.86 ...
# $ p_value  : num [1:2814] 1.94e-04 2.15e-02 9.48e-03 4.10e-01 4.03e-15 ...

table(acts$source)
# Androgen     EGFR Estrogen  Hypoxia JAK-STAT     MAPK     NFkB      p53     PI3K     TGFb     TNFa    Trail 
# 201      201      201      201      201      201      201      201      201      201      201      201 
# VEGF      WNT 
# 201      201 

saveRDS(acts,'C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\LUAD_low_risk_activity_scores.rds')


#------------------------------------------------------------------------------------------------------------
# Kinase activity inference
#------------------------------------------------------------------------------------------------------------
rm(list=ls())
set.seed(2025)

# Install (run once)
# if(!requireNamespace("BiocManager")) install.packages("BiocManager")
# BiocManager::install(c("decoupleR", "OmnipathR", "org.Hs.eg.db", "limma"))

.libPaths(c('C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/R_packages/',.libPaths()))
lib <- 'C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/R_packages/'

# Load
library(BiocManager)
library(decoupleR)
library(OmnipathR)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(limma)
library(tidyverse)

# Extract the normalized log-transformed counts

expr <- t(as.matrix(readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\log_cts_high_risk_luad.rds')))
head(expr[,1:5])
# Check for duplicate column names
sum(duplicated(colnames(expr)))
# FIX: Make column names unique (samples)
colnames(expr) <- make.unique(colnames(expr))
# Check for duplicate column names
sum(duplicated(colnames(expr)))
# Check for duplicate row names  
sum(duplicated(rownames(expr)))

# # New function, equivalent but modern
# interactions <- omnipath_interactions(resources = 'ARN')
# 
# # Check columns
# colnames(interactions)[1:10]


# Build kinase-substrate network
ks <- omnipath_interactions(datasets = c("kinaseextra", "omnipath")) %>%
  filter(!is.na(source_genesymbol) & !is.na(target_genesymbol)) %>%
  transmute(
    source = source_genesymbol,   # kinase
    target = target_genesymbol,   # downstream gene
    mor = 1                       # +1 (unknown directionality)
  ) %>%
  distinct()

# keep only targets present in your data
ks <- ks %>% filter(target %in% rownames(expr))%>%
  dplyr::group_by(source) %>%
  dplyr::filter(dplyr::n() >= 50) %>%   # only kinases with ≥10 targets
  dplyr::ungroup()

# Check if network is valid
cat("Network size:", nrow(ks), "interactions\n")
cat("Number of unique kinases:", length(unique(ks$source)), "\n")
cat("Number of unique targets:", length(unique(ks$target)), "\n")



library(dplyr)

# Step 1: Basic quality filters
ks_clean <- ks %>%
  filter(target %in% rownames(expr)) %>%
  distinct(source, target, .keep_all = TRUE)

# Step 2: Count substrates
network_size <- ks_clean %>%
  group_by(source) %>%
  summarise(n_substrates = n())

# Step 3: Filter by expression
expressed_kinases <- intersect(unique(ks_clean$source), rownames(expr))
kinase_expression <- rowMeans(expr[expressed_kinases, , drop = FALSE])
well_expressed <- names(kinase_expression[kinase_expression > median(kinase_expression)])

# Step 4: Apply combined filters
ks_filtered <- ks_clean %>%
  filter(
    source %in% well_expressed,           # Expressed kinases
    source %in% network_size$source[network_size$n_substrates >= 10 &  # 10-100 substrates
                                      network_size$n_substrates <= 100]
  )

cat("Original kinases:", length(unique(ks$source)), "\n")
cat("Filtered kinases:", length(unique(ks_filtered$source)), "\n")

library(decoupleR)

# expr <- expr[!duplicated(rownames(expr)), ]    # remove duplicate genes
# ks   <- ks  %>% distinct(source, target, .keep_all = TRUE)


kinase_scores <- run_ulm( 
  # or run_mlm --: error at least 135 sources in the network are colinear with other sources.
  # or run_ulm
  mat = expr,
  network = ks,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize  = 10
)


# Calculate variance of kinase activities
kinase_variance <- kinase_scores %>%
  group_by(source) %>%
  summarise(score_var = var(score, na.rm = TRUE)) %>%
  arrange(desc(score_var))

# Keep top 14 most variable kinases
top_kinases <- kinase_variance %>%
  top_n(14, score_var) %>%
  pull(source)

ks_filtered <- ks %>%
  filter(source %in% top_kinases)

# Extract scores for top 14 kinases
top_kinase_scores <- kinase_scores %>%
  filter(source %in% top_kinases)

grouped <- top_kinase_scores %>% group_by(source) %>%
  group_split()
length(grouped)

top_kinase_scores$source %>% unique()

# "AKT1"    "ATM"     "CDK1"    "CDK2"    "CSNK2A1" "EGF"     
# "EGFR"    "GSK3B"   "MAPK1"   "MAPK14"  "MAPK3"  
# "MAPK8"   "PRKCA"   "SRC"

# Summary stats for each kinase
kinase_summary <- top_kinase_scores %>%
  group_by(source) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    min_score = min(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    variance = var(score, na.rm = TRUE)
  ) %>%
  arrange(desc(variance))

print(kinase_summary)

library(dplyr)
# Long format (good for ggplot)
top_kinase_df <- kinase_scores %>%
  filter(source %in% top_kinases) %>%
  dplyr::select(source, condition, score)

head(top_kinase_df)
saveRDS(top_kinase_df, 'C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\top_kinase_df_highrisk_luad.rds')


# ggplot of top_kinase_df
library(ggplot2)
ggplot(top_kinase_df, aes(x = condition, y = score, color = source)) +
  geom_point(alpha = 0.6) +
  geom_line(aes(group = source), alpha = 0.3) +
  labs(title = "Kinase Activity Scores for Top 14 Variable Kinases",
       x = "Samples",
       y = "Activity Score") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~ source, ncol = 3)


library(ggplot2)
library(dplyr)

# Calculate mean activity for each kinase
total_kinase_activity <- kinase_scores %>%
  filter(source %in% top_kinases) %>%
  group_by(source) %>%
  summarise(
    mean_activity = mean(score, na.rm = TRUE),
    sd_activity = sd(score, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_activity))

top_kinase_df <- kinase_scores %>%
  filter(source %in% top_kinases)

library(plotly)

p <- ggplot(top_kinase_df, aes(x = reorder(source, score, FUN = median), y = score, fill = source)) +
  geom_boxplot() +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Kinase Activity Distribution",
       x = "Kinase",
       y = "Activity Score")

ggplotly(p)


#------------------------------------------
# LUAD low risk 
#------------------------------------------


expr <- t(as.matrix(readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\log_cts_low_risk_luad.rds')))
head(expr[,1:5])
# Check for duplicate column names
sum(duplicated(colnames(expr)))
# FIX: Make column names unique (samples)
colnames(expr) <- make.unique(colnames(expr))
# Check for duplicate column names
sum(duplicated(colnames(expr)))
# Check for duplicate row names  
sum(duplicated(rownames(expr)))

# # New function, equivalent but modern
# interactions <- omnipath_interactions(resources = 'ARN')
# 
# # Check columns
# colnames(interactions)[1:10]


# Build kinase-substrate network
ks <- omnipath_interactions(datasets = c("kinaseextra", "omnipath")) %>%
  filter(!is.na(source_genesymbol) & !is.na(target_genesymbol)) %>%
  transmute(
    source = source_genesymbol,   # kinase
    target = target_genesymbol,   # downstream gene
    mor = 1                       # +1 (unknown directionality)
  ) %>%
  distinct()

# keep only targets present in your data
ks <- ks %>% filter(target %in% rownames(expr))%>%
  dplyr::group_by(source) %>%
  dplyr::filter(dplyr::n() >= 50) %>%   # only kinases with ≥10 targets
  dplyr::ungroup()

# Check if network is valid
cat("Network size:", nrow(ks), "interactions\n")
cat("Number of unique kinases:", length(unique(ks$source)), "\n")
cat("Number of unique targets:", length(unique(ks$target)), "\n")



library(dplyr)

# Step 1: Basic quality filters
ks_clean <- ks %>%
  filter(target %in% rownames(expr)) %>%
  distinct(source, target, .keep_all = TRUE)

# Step 2: Count substrates
network_size <- ks_clean %>%
  group_by(source) %>%
  summarise(n_substrates = n())

# Step 3: Filter by expression
expressed_kinases <- intersect(unique(ks_clean$source), rownames(expr))
kinase_expression <- rowMeans(expr[expressed_kinases, , drop = FALSE])
well_expressed <- names(kinase_expression[kinase_expression > median(kinase_expression)])

# Step 4: Apply combined filters
ks_filtered <- ks_clean %>%
  filter(
    source %in% well_expressed,           # Expressed kinases
    source %in% network_size$source[network_size$n_substrates >= 10 &  # 10-100 substrates
                                      network_size$n_substrates <= 100]
  )

cat("Original kinases:", length(unique(ks$source)), "\n")
cat("Filtered kinases:", length(unique(ks_filtered$source)), "\n")

library(decoupleR)

# expr <- expr[!duplicated(rownames(expr)), ]    # remove duplicate genes
# ks   <- ks  %>% distinct(source, target, .keep_all = TRUE)


kinase_scores <- run_ulm( 
  # or run_mlm --: error at least 135 sources in the network are colinear with other sources.
  # or run_ulm
  mat = expr,
  network = ks,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize  = 10
)


# Calculate variance of kinase activities
kinase_variance <- kinase_scores %>%
  group_by(source) %>%
  summarise(score_var = var(score, na.rm = TRUE)) %>%
  arrange(desc(score_var))

# Keep top 14 most variable kinases
top_kinases <- kinase_variance %>%
  top_n(14, score_var) %>%
  pull(source)

ks_filtered <- ks %>%
  filter(source %in% top_kinases)

# Extract scores for top 14 kinases
top_kinase_scores <- kinase_scores %>%
  filter(source %in% top_kinases)

grouped <- top_kinase_scores %>% group_by(source) %>%
  group_split()
length(grouped)

top_kinase_scores$source %>% unique()

# "AKT1"    "ATM"     "CDK1"    "CDK2"    "CSNK2A1" "EGF"
# "GSK3B"   "MAPK1"   "MAPK14"  "MAPK3"  
# "MAPK8"   "PRKCA"   "RPS6KA3" "SRC"  

# Summary stats for each kinase
kinase_summary <- top_kinase_scores %>%
  group_by(source) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    min_score = min(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    variance = var(score, na.rm = TRUE)
  ) %>%
  arrange(desc(variance))

print(kinase_summary)

# Long format (good for ggplot)
top_kinase_df <- kinase_scores %>%
  filter(source %in% top_kinases) %>%
  select(source, condition, score)

head(top_kinase_df)

# ggplot of top_kinase_df
library(ggplot2)
ggplot(top_kinase_df, aes(x = condition, y = score, color = source)) +
  geom_point(alpha = 0.6) +
  geom_line(aes(group = source), alpha = 0.3) +
  labs(title = "Kinase Activity Scores for Top 14 Variable Kinases - LUSC",
       x = "Samples",
       y = "Activity Score") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~ source, ncol = 3)


library(ggplot2)
library(dplyr)

# Calculate mean activity for each kinase
total_kinase_activity <- kinase_scores %>%
  filter(source %in% top_kinases) %>%
  group_by(source) %>%
  summarise(
    mean_activity = mean(score, na.rm = TRUE),
    sd_activity = sd(score, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_activity))

top_kinase_df <- kinase_scores %>%
  filter(source %in% top_kinases)

saveRDS(top_kinase_df, 'C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\top_kinase_df_lowrisk_luad.rds')


library(plotly)

p <- ggplot(top_kinase_df, aes(x = reorder(source, score, FUN = median), y = score, fill = source)) +
  #order with median score
  #ggplot(top_kinase_df, aes(x = reorder(source, score, FUN = mean), y = score, fill = source))+
  geom_boxplot() +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Kinase Activity Distribution - LUSC",
       x = "Kinase",
       y = "Activity Score")

ggplotly(p)



#------------------------------------------------------------------------------------------------------------
# Transcriptome factor activity inference
#------------------------------------------------------------------------------------------------------------
rm(list=ls())
set.seed(2025)

# Install (run once)
# if(!requireNamespace("BiocManager")) install.packages("BiocManager")
# BiocManager::install(c("decoupleR", "OmnipathR", "org.Hs.eg.db", "limma"))

.libPaths(c('C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/R_packages/',.libPaths()))
lib <- 'C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/R_packages/'

# Load
library(BiocManager)
library(decoupleR)
library(OmnipathR)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(limma)

# Extract the normalized log-transformed counts

expr <- t(as.matrix(readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\log_cts_high_risk_luad.rds')))
head(expr[,1:5])
# Check for duplicate column names
sum(duplicated(colnames(expr)))
# FIX: Make column names unique (samples)
colnames(expr) <- make.unique(colnames(expr))
# Check for duplicate column names
sum(duplicated(colnames(expr)))
# Check for duplicate row names  
sum(duplicated(rownames(expr)))



# Get dorothea regulons (TF-target gene relationships)
# Dorothea has confidence levels: A (highest) to E (lowest)

net <- get_dorothea(organism = 'human', levels = c('A', 'B'))

# # 2. Run TF activity inference
# tf_acts <- run_mlm(
#   mat = expr,
#   network = net,
#   minsize = 5
# )

# Recommended: Use multiple methods and compare
tf_acts_mlm <- run_mlm(
  mat = expr,
  network = net,
  .source = 'source',
  .target = 'target',
  .mor = 'mor',
  minsize = 10  # minimum number of targets per TF
)
library(ggplot2)

# # Find top 10 most active TFs in each sample
# tf_acts_mlm %>%
#   group_by(condition) %>%
#   slice_max(order_by = score, n = 10) %>%
#   ggplot(aes(x = reorder(source, score), y = score, fill = score)) +
#   geom_bar(stat = 'identity') +
#   facet_wrap(~source, scales = 'free_y') +
#   coord_flip() +
#   theme_minimal() +
#   labs(x = 'Transcription Factor', y = 'Activity Score')


library(textshape)

# 3. Convert to matrix
tf_matrix <- tf_acts_mlm %>%
  pivot_wider(
    id_cols = 'source',
    names_from = 'condition',
    values_from = 'score'
  ) %>%
  column_to_rownames('source') %>%
  as.matrix()

# Now you have TFs as rows, samples as columns
dim(tf_matrix)

# 4. Plot heatmap of TF activities

# BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)

# Create color scale
col_fun <- colorRamp2(
  c(-3, 0, 3),
  c("blue", "white", "red")
)

# Plot heatmap
Heatmap(
  tf_matrix,
  name = "TF Activity",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 5)
)
#------------------------------------------


# Calculate variance of kinase activities
tf_variance <- tf_acts_mlm %>%
  group_by(source) %>%
  summarise(score_var = var(score, na.rm = TRUE)) %>%
  arrange(desc(score_var))

# Keep top 14 most variable kinases
top_tf <- tf_variance %>%
  top_n(14, score_var) %>%
  pull(source)

# Extract scores for top 14 kinases
top_tf_scores <- tf_acts_mlm %>%
  filter(source %in% top_tf)

grouped <- top_tf_scores %>% group_by(source) %>%
  group_split()
length(grouped)

top_tf_scores$source %>% unique()

# "MYC"    "TP53"   "STAT1"  "HNF1A"  "HNF4A"  "EGR1"
# "E2F4"   "SPI1"   "CDX2"   "RFX5"  
# "SOX2"   "KMT2A"  "PRDM14" "ZNF263"

# Summary stats for each kinase
tf_summary <- top_tf_scores %>%
  group_by(source) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    min_score = min(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    variance = var(score, na.rm = TRUE)
  ) %>%
  arrange(desc(variance))

print(tf_summary)

# Long format (good for ggplot)
top_tf_df <- top_tf_scores %>%
  filter(source %in% top_tf) %>%
  select(source, condition, score)

head(top_tf_df)

# ggplot of top_kinase_df
library(ggplot2)
ggplot(top_tf_df, aes(x = condition, y = score, color = source)) +
  geom_point(alpha = 0.6) +
  geom_line(aes(group = source), alpha = 0.3) +
  labs(title = "Kinase Activity Scores for Top 14 Variable Kinases - LUSC",
       x = "Samples",
       y = "Activity Score") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~ source, ncol = 3)


library(ggplot2)
library(dplyr)

# Calculate mean activity for each kinase
total_tf_activity <- top_tf_scores %>%
  filter(source %in% top_tf_scores) %>%
  group_by(source) %>%
  summarise(
    mean_activity = mean(score, na.rm = TRUE),
    sd_activity = sd(score, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_activity))

top_tf_df <- top_tf_scores %>%
  filter(source %in% top_tf)

saveRDS(top_tf_df, 'C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\top_tf_df_highrisk_luad.rds')

library(plotly)

p <- ggplot(top_tf_df, aes(x = reorder(source, score, FUN = median), y = score, fill = source)) +
  #order with median score
  #ggplot(top_kinase_df, aes(x = reorder(source, score, FUN = mean), y = score, fill = source))+
  geom_boxplot() +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "TF Activity Distribution - LUAD",
       x = "TF",
       y = "Activity Score")

ggplotly(p)

#---------------------------------------
# LUAD low risk
#---------------------------------------


expr <- t(as.matrix(readRDS('C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\log_cts_low_risk_luad.rds')))
head(expr[,1:5])
# Check for duplicate column names
sum(duplicated(colnames(expr)))
# FIX: Make column names unique (samples)
colnames(expr) <- make.unique(colnames(expr))
# Check for duplicate column names
sum(duplicated(colnames(expr)))
# Check for duplicate row names  
sum(duplicated(rownames(expr)))



# Get dorothea regulons (TF-target gene relationships)
# Dorothea has confidence levels: A (highest) to E (lowest)

net <- get_dorothea(organism = 'human', levels = c('A', 'B'))

# # 2. Run TF activity inference
# tf_acts <- run_mlm(
#   mat = expr,
#   network = net,
#   minsize = 5
# )

# Recommended: Use multiple methods and compare
tf_acts_mlm <- run_mlm(
  mat = expr,
  network = net,
  .source = 'source',
  .target = 'target',
  .mor = 'mor',
  minsize = 10  # minimum number of targets per TF
)
library(ggplot2)

# # Find top 10 most active TFs in each sample
# tf_acts_mlm %>%
#   group_by(condition) %>%
#   slice_max(order_by = score, n = 10) %>%
#   ggplot(aes(x = reorder(source, score), y = score, fill = score)) +
#   geom_bar(stat = 'identity') +
#   facet_wrap(~source, scales = 'free_y') +
#   coord_flip() +
#   theme_minimal() +
#   labs(x = 'Transcription Factor', y = 'Activity Score')


library(textshape)

# 3. Convert to matrix
tf_matrix <- tf_acts_mlm %>%
  pivot_wider(
    id_cols = 'source',
    names_from = 'condition',
    values_from = 'score'
  ) %>%
  column_to_rownames('source') %>%
  as.matrix()

# Now you have TFs as rows, samples as columns
dim(tf_matrix)

# 4. Plot heatmap of TF activities

# BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)

# Create color scale
col_fun <- colorRamp2(
  c(-3, 0, 3),
  c("blue", "white", "red")
)

# Plot heatmap
Heatmap(
  tf_matrix,
  name = "TF Activity",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 5)
)
#------------------------------------------


# Calculate variance of kinase activities
tf_variance <- tf_acts_mlm %>%
  group_by(source) %>%
  summarise(score_var = var(score, na.rm = TRUE)) %>%
  arrange(desc(score_var))

# Keep top 14 most variable kinases
top_tf <- tf_variance %>%
  top_n(14, score_var) %>%
  pull(source)

# Extract scores for top 14 kinases
top_tf_scores <- tf_acts_mlm %>%
  filter(source %in% top_tf)

grouped <- top_tf_scores %>% group_by(source) %>%
  group_split()
length(grouped)

top_tf_scores$source %>% unique()

# "TP53"   "STAT1"  "HNF1A"  "HNF4A"  "EGR1"
# "E2F4"   "SPI1"   "NFE2L2" "STAT2"  "RFX5"  
# "SOX2"   "KMT2A"  "PRDM14" "ZNF263"

# Summary stats for each kinase
tf_summary <- top_tf_scores %>%
  group_by(source) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    min_score = min(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    variance = var(score, na.rm = TRUE)
  ) %>%
  arrange(desc(variance))

print(tf_summary)

# Long format (good for ggplot)
top_tf_df <- top_tf_scores %>%
  filter(source %in% top_tf) %>%
  select(source, condition, score)

head(top_tf_df)

# ggplot of top_kinase_df
library(ggplot2)
ggplot(top_tf_df, aes(x = condition, y = score, color = source)) +
  geom_point(alpha = 0.6) +
  geom_line(aes(group = source), alpha = 0.3) +
  labs(title = "Kinase Activity Scores for Top 14 Variable Kinases - LUSC",
       x = "Samples",
       y = "Activity Score") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~ source, ncol = 3)


library(ggplot2)
library(dplyr)

# Calculate mean activity for each kinase
total_kinase_activity <- top_tf_scores %>%
  filter(source %in% top_tf_scores) %>%
  group_by(source) %>%
  summarise(
    mean_activity = mean(score, na.rm = TRUE),
    sd_activity = sd(score, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_activity))

top_tf_df <- top_tf_scores %>%
  filter(source %in% top_tf)

saveRDS(top_tf_df, 'C:\\Users\\bhadury\\University of Michigan Dropbox\\Sagnik Bhadury\\CODES\\SBTurboCodes\\H&E Lung TCGA Graph Dataset\\TCGA-LUAD-data\\top_tf_df_lowrisk_luad.rds')


library(plotly)

p <- ggplot(top_tf_df, aes(x = reorder(source, score, FUN = median), y = score, fill = source)) +
  #order with median score
  #ggplot(top_kinase_df, aes(x = reorder(source, score, FUN = mean), y = score, fill = source))+
  geom_boxplot() +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "TF Activity Distribution - LUAD",
       x = "TF",
       y = "Activity Score")

ggplotly(p)






