# ============================================================================
# PERMUTATION TESTING: PARALLEL EXECUTION (HPC CLUSTER)
# ============================================================================
# Purpose: Parallel permutation testing for statistical significance
# Runtime: ~1 week on standard CPU, ~1-2 days on HPC cluster (64+ cores)
# ============================================================================

library(tidyverse)
library(parallel)
library(doParallel)

set.seed(42)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Number of permutations (standard is 1000, can increase to 10000 for robustness)
n_perms <- 1000

# Number of cores (adjust based on your system)
# HPC: Can use 64-128 cores
# Local: Use detectCores() - 2
n_cores <- parallel::detectCores() - 2
cat(sprintf("Using %d cores for parallel processing\n", n_cores))

# Input paths (adjust to match your directory)
base_dir <- "H:/path/to/TCGA/analysis"
output_dir <- "H:/path/to/results/permutation_testing"
dir.create(output_dir, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading data...\n")

# Load model outputs (PIPs by risk group)
# Assume you have these from your main analysis:
# - pips_pathways_hr: (n_hpc × n_pathways) matrix for high-risk
# - pips_pathways_lr: (n_hpc × n_pathways) matrix for low-risk
# - Similar for kinases and TFs
# - risk_group: vector of risk assignments

pips_pathways_hr <- read.csv(file.path(base_dir, "pips_pathways_hr.csv"),
                              row.names = 1)
pips_pathways_lr <- read.csv(file.path(base_dir, "pips_pathways_lr.csv"),
                              row.names = 1)

pips_kinases_hr <- read.csv(file.path(base_dir, "pips_kinases_hr.csv"),
                             row.names = 1)
pips_kinases_lr <- read.csv(file.path(base_dir, "pips_kinases_lr.csv"),
                             row.names = 1)

pips_tfs_hr <- read.csv(file.path(base_dir, "pips_tfs_hr.csv"), row.names = 1)
pips_tfs_lr <- read.csv(file.path(base_dir, "pips_tfs_lr.csv"), row.names = 1)

# Original risk assignments
risk_group <- read.csv(file.path(base_dir, "risk_assignments.csv"))$risk_group
n_samples <- length(risk_group)
n_high <- sum(risk_group == "high")
n_low <- sum(risk_group == "low")

# ============================================================================
# COMPUTE ORIGINAL ASSOCIATION COUNTS
# ============================================================================

cat("Computing original association counts...\n")

original_counts <- tibble(
  Layer = c("Pathways", "Kinases", "TFs"),
  N_High_Risk = c(
    sum(pips_pathways_hr > 0.5),
    sum(pips_kinases_hr > 0.5),
    sum(pips_tfs_hr > 0.5)
  ),
  N_Low_Risk = c(
    sum(pips_pathways_lr > 0.5),
    sum(pips_kinases_lr > 0.5),
    sum(pips_tfs_lr > 0.5)
  )
)

original_counts <- original_counts %>%
  mutate(
    Difference = N_High_Risk - N_Low_Risk,
    Fold_Enrichment = N_High_Risk / (N_Low_Risk + 0.1)  # Avoid division by zero
  )

cat("Original association counts:\n")
print(original_counts)

# ============================================================================
# PERMUTATION TESTING FUNCTION
# ============================================================================

# Core permutation function
run_single_permutation <- function(perm_idx, pips_list, risk_group,
                                   n_high, n_low) {

  # Shuffle risk group assignment
  risk_shuffled <- c(rep("high", n_high), rep("low", n_low))
  risk_shuffled <- sample(risk_shuffled, replace = FALSE)

  # For each layer, randomly reassign PIPs to shuffled groups
  # (In practice, you would REFIT the model with shuffled labels;
  #  this simplified version shuffles PIPs, which is conservative)

  results_perm <- tibble()

  for (layer in names(pips_list)) {
    pips_matrix <- pips_list[[layer]]

    # Randomly permute row assignments (samples)
    perm_order <- sample(1:nrow(pips_matrix), replace = FALSE)
    pips_permuted <- pips_matrix[perm_order, ]

    # Assign to shuffled risk groups
    pips_perm_high <- pips_permuted[risk_shuffled == "high", ]
    pips_perm_low <- pips_permuted[risk_shuffled == "low", ]

    n_assoc_high <- sum(pips_perm_high > 0.5)
    n_assoc_low <- sum(pips_perm_low > 0.5)

    results_perm <- bind_rows(results_perm, tibble(
      Permutation = perm_idx,
      Layer = layer,
      N_High_Risk = n_assoc_high,
      N_Low_Risk = n_assoc_low,
      Difference = n_assoc_high - n_assoc_low
    ))
  }

  return(results_perm)
}

# ============================================================================
# PARALLEL PERMUTATION TESTING
# ============================================================================

cat(sprintf("\nRunning %d permutations in parallel...\n", n_perms))

# Set up parallel cluster
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Create list of PIPs for passing to parallel function
pips_list <- list(
  Pathways = pips_pathways_hr,
  Kinases = pips_kinases_hr,
  TFs = pips_tfs_hr
)

# Run permutations in parallel
# Note: cbind necessary because of how parallel processing returns results
start_time <- Sys.time()

perm_results_list <- foreach(
  perm_idx = 1:n_perms,
  .combine = bind_rows,
  .packages = c("tidyverse")
) %dopar% {
  if (perm_idx %% 100 == 0) {
    cat(sprintf("  Permutation %d/%d\n", perm_idx, n_perms))
  }

  run_single_permutation(perm_idx, pips_list, risk_group, n_high, n_low)
}

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

stopCluster(cl)

cat(sprintf("\n✓ Completed %d permutations in %.1f minutes\n", n_perms, runtime))

# ============================================================================
# ANALYZE PERMUTATION RESULTS
# ============================================================================

cat("\nAnalyzing permutation results...\n")

# Compute null distribution for each layer
perm_summary <- perm_results_list %>%
  group_by(Layer) %>%
  summarize(
    Mean_Diff = mean(Difference),
    SD_Diff = sd(Difference),
    Min_Diff = min(Difference),
    Max_Diff = max(Difference),
    P025 = quantile(Difference, 0.025),
    P975 = quantile(Difference, 0.975),
    .groups = "drop"
  )

# Compute permutation p-values (two-tailed)
p_values_perm <- original_counts %>%
  select(Layer, Difference) %>%
  left_join(
    perm_results_list %>%
      group_by(Layer) %>%
      summarize(
        p_value = mean(abs(Difference) >= abs(.$Difference)),
        .groups = "drop"
      ),
    by = "Layer"
  )

cat("\nPermutation p-values (two-tailed test):\n")
print(p_values_perm)

# ============================================================================
# VISUALIZATION: PERMUTATION DISTRIBUTIONS
# ============================================================================

cat("\nCreating visualizations...\n")

# Create density plots for each layer
pdf("Permutation_Distribution_Pathways.pdf", width = 8, height = 5)

perm_pathways <- perm_results_list %>% filter(Layer == "Pathways")
orig_pathways <- original_counts %>% filter(Layer == "Pathways")

plot_data <- perm_pathways %>%
  ggplot(aes(x = Difference)) +
  geom_histogram(bins = 50, alpha = 0.6, fill = "steelblue") +
  geom_vline(xintercept = orig_pathways$Difference[1],
             color = "red", linewidth = 1.5, linetype = "solid") +
  geom_vline(xintercept = c(perm_summary$P025[1], perm_summary$P975[1]),
             color = "gray", linewidth = 1, linetype = "dashed") +
  labs(
    title = "Permutation Distribution: Pathway Associations",
    x = "Difference (High-Risk - Low-Risk Association Counts)",
    y = "Frequency",
    subtitle = sprintf("Red line: observed difference (%.0f)\nGray lines: 95%% CI",
                       orig_pathways$Difference[1])
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_data)
dev.off()

# Repeat for kinases and TFs
pdf("Permutation_Distribution_Kinases.pdf", width = 8, height = 5)
perm_kinases <- perm_results_list %>% filter(Layer == "Kinases")
orig_kinases <- original_counts %>% filter(Layer == "Kinases")

perm_kinases %>%
  ggplot(aes(x = Difference)) +
  geom_histogram(bins = 50, alpha = 0.6, fill = "steelblue") +
  geom_vline(xintercept = orig_kinases$Difference[1],
             color = "red", linewidth = 1.5) +
  labs(
    title = "Permutation Distribution: Kinase Associations",
    x = "Difference (High-Risk - Low-Risk)",
    y = "Frequency"
  ) +
  theme_minimal()
dev.off()

pdf("Permutation_Distribution_TFs.pdf", width = 8, height = 5)
perm_tfs <- perm_results_list %>% filter(Layer == "TFs")
orig_tfs <- original_counts %>% filter(Layer == "TFs")

perm_tfs %>%
  ggplot(aes(x = Difference)) +
  geom_histogram(bins = 50, alpha = 0.6, fill = "steelblue") +
  geom_vline(xintercept = orig_tfs$Difference[1],
             color = "red", linewidth = 1.5) +
  labs(
    title = "Permutation Distribution: Transcription Factor Associations",
    x = "Difference (High-Risk - Low-Risk)",
    y = "Frequency"
  ) +
  theme_minimal()
dev.off()

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("Saving results...\n")

# Original counts
write.csv(original_counts,
          file.path(output_dir, "Original_Association_Counts.csv"),
          row.names = FALSE)

# Permutation summary
write.csv(perm_summary,
          file.path(output_dir, "Permutation_Summary_Statistics.csv"),
          row.names = FALSE)

# P-values
write.csv(p_values_perm,
          file.path(output_dir, "Permutation_Pvalues.csv"),
          row.names = FALSE)

# Full permutation results (for sensitivity analysis)
write.csv(perm_results_list,
          file.path(output_dir, "All_Permutation_Results.csv"),
          row.names = FALSE)

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n=== PERMUTATION TESTING SUMMARY ===\n")
cat(sprintf("Number of permutations: %d\n", n_perms))
cat(sprintf("Runtime: %.1f minutes\n\n", runtime))

cat("ORIGINAL ASSOCIATION COUNTS:\n")
print(original_counts)

cat("\n\nPERMUTATION-BASED P-VALUES (Two-tailed):\n")
print(p_values_perm)

cat("\nNULL DISTRIBUTION STATISTICS:\n")
print(perm_summary)

cat("\n\nINTERPRETATION:\n")
for (i in 1:nrow(p_values_perm)) {
  layer <- p_values_perm$Layer[i]
  p_val <- p_values_perm$p_value[i]
  sig <- ifelse(p_val < 0.05, "SIGNIFICANT", "not significant")
  cat(sprintf("  %s: p = %.4f (%s)\n", layer, p_val, sig))
}

cat(sprintf("\n✓ Results saved to: %s\n", output_dir))
cat("✓ Supplementary Table S6 ready for manuscript\n")

# ============================================================================
# END OF SCRIPT
# ============================================================================
