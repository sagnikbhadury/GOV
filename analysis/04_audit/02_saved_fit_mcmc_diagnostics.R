# Single-chain mixing diagnostics for the 84 saved TCGA-LUAD model objects.
#
# The saved objects contain one MCMC chain per fitted activity. Consequently,
# these summaries cannot establish convergence and cannot provide multi-chain
# R-hat. They quantify the information in the available post-burn-in traces by
# reporting lag-1 autocorrelation and an initial-positive-sequence estimate of
# effective sample size (ESS) for the five saved scalar chains.

options(stringsAsFactors = FALSE)

root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(root, "luad"))) {
  stop("Run this script from the repository root.")
}

out_dir <- file.path(root, "reviewer_reanalysis", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

result_files <- list.files(
  file.path(root, "luad"),
  pattern = "^(HR|LR)_luad_.*result1000_.*_activity\\.rds$",
  full.names = TRUE
)
stopifnot(length(result_files) == 84L)

parse_fit_name <- function(path) {
  nm <- basename(path)
  risk <- sub("^((HR|LR)).*$", "\\1", nm)
  if (grepl("_kinase_", nm, fixed = TRUE)) {
    layer <- "Kinase"
    feature <- sub("^(HR|LR)_luad_kinase_lsGPRVS_result1000_", "", nm)
  } else if (grepl("_TF_", nm, fixed = TRUE)) {
    layer <- "Transcription factor"
    feature <- sub("^(HR|LR)_luad_TF_lsGPRVS_result1000_", "", nm)
  } else {
    layer <- "Pathway"
    feature <- sub("^(HR|LR)_luad_lsGPRVS_result1000_", "", nm)
  }
  feature <- sub("_activity\\.rds$", "", feature)
  list(risk = risk, layer = layer, feature = feature)
}

lag1_acf <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 3L || !is.finite(sd(x)) || sd(x) == 0) return(NA_real_)
  unname(acf(x, plot = FALSE, lag.max = 1L, demean = TRUE)$acf[[2L]])
}

ess_ips <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 3L || !is.finite(sd(x)) || sd(x) == 0) return(NA_real_)
  rho <- as.numeric(acf(
    x, plot = FALSE, lag.max = min(n - 1L, floor(n / 2)),
    type = "correlation", demean = TRUE
  )$acf)[-1L]
  if (!length(rho)) return(as.numeric(n))

  # Geyer's initial positive sequence: retain adjacent autocorrelation pairs
  # until the first pair whose sum is non-positive.
  pair_count <- floor(length(rho) / 2L)
  if (pair_count == 0L) return(as.numeric(n))
  pair_sums <- rho[2L * seq_len(pair_count) - 1L] +
    rho[2L * seq_len(pair_count)]
  nonpositive <- which(!is.finite(pair_sums) | pair_sums <= 0)
  keep_pairs <- if (length(nonpositive)) nonpositive[[1L]] - 1L else pair_count
  tau <- if (keep_pairs > 0L) {
    1 + 2 * sum(rho[seq_len(2L * keep_pairs)])
  } else {
    1
  }
  min(as.numeric(n), max(1, n / tau))
}

parameters <- c(
  post_log_likelihood = "post.loglik",
  residual_precision = "TAU",
  kernel_variance = "PSI1",
  density_lengthscale = "PSIA",
  network_lengthscale = "PSIU"
)

rows <- list()
gamma_rows <- list()
k <- 0L
g <- 0L
for (path in result_files) {
  info <- parse_fit_name(path)
  fit <- readRDS(path)
  keep <- seq.int(fit$burn + 1L, fit$nrun)
  for (parameter in names(parameters)) {
    k <- k + 1L
    chain <- fit[[parameters[[parameter]]]][keep]
    rows[[k]] <- data.frame(
      risk_group = info$risk,
      molecular_layer = info$layer,
      activity = info$feature,
      parameter = parameter,
      post_burnin_draws = length(chain),
      lag1_autocorrelation = lag1_acf(chain),
      effective_sample_size = ess_ips(chain)
    )
  }
  for (node in seq_len(nrow(fit$GAMMA))) {
    g <- g + 1L
    chain <- fit$GAMMA[node, keep]
    pip <- mean(chain)
    ess <- ess_ips(chain)
    gamma_rows[[g]] <- data.frame(
      risk_group = info$risk,
      molecular_layer = info$layer,
      activity = info$feature,
      hpc_index = node,
      pip = pip,
      constant_chain = length(unique(chain)) == 1L,
      lag1_autocorrelation = lag1_acf(chain),
      effective_sample_size = ess,
      pip_monte_carlo_se = if (is.finite(ess)) {
        sqrt(pip * (1 - pip) / ess)
      } else {
        NA_real_
      }
    )
  }
}

diagnostics <- do.call(rbind, rows)
write.csv(
  diagnostics,
  file.path(out_dir, "saved_fit_mcmc_diagnostics_by_activity.csv"),
  row.names = FALSE
)

gamma_diagnostics <- do.call(rbind, gamma_rows)
write.csv(
  gamma_diagnostics,
  file.path(out_dir, "saved_fit_pip_mcmc_diagnostics.csv"),
  row.names = FALSE
)

groups <- split(
  diagnostics,
  interaction(
    diagnostics$risk_group,
    diagnostics$molecular_layer,
    diagnostics$parameter,
    drop = TRUE
  )
)
summary_rows <- lapply(groups, function(dat) data.frame(
  risk_group = dat$risk_group[[1L]],
  molecular_layer = dat$molecular_layer[[1L]],
  parameter = dat$parameter[[1L]],
  models = nrow(dat),
  median_lag1_autocorrelation = median(
    dat$lag1_autocorrelation, na.rm = TRUE
  ),
  median_effective_sample_size = median(
    dat$effective_sample_size, na.rm = TRUE
  ),
  minimum_effective_sample_size = min(
    dat$effective_sample_size, na.rm = TRUE
  ),
  models_with_ess_at_least_100 = sum(
    dat$effective_sample_size >= 100, na.rm = TRUE
  )
))
diagnostic_summary <- do.call(rbind, summary_rows)
rownames(diagnostic_summary) <- NULL
diagnostic_summary <- diagnostic_summary[order(
  diagnostic_summary$molecular_layer,
  diagnostic_summary$risk_group,
  diagnostic_summary$parameter
), ]

write.csv(
  diagnostic_summary,
  file.path(out_dir, "saved_fit_mcmc_diagnostics_summary.csv"),
  row.names = FALSE
)

gamma_groups <- split(
  gamma_diagnostics,
  interaction(
    gamma_diagnostics$risk_group,
    gamma_diagnostics$molecular_layer,
    drop = TRUE
  )
)
gamma_summary <- do.call(rbind, lapply(gamma_groups, function(dat) {
  moving <- dat[!dat$constant_chain, , drop = FALSE]
  data.frame(
    risk_group = dat$risk_group[[1L]],
    molecular_layer = dat$molecular_layer[[1L]],
    indicator_chains = nrow(dat),
    constant_indicator_chains = sum(dat$constant_chain),
    nonconstant_indicator_chains = nrow(moving),
    median_ess_nonconstant = if (nrow(moving)) {
      median(moving$effective_sample_size, na.rm = TRUE)
    } else NA_real_,
    minimum_ess_nonconstant = if (nrow(moving)) {
      min(moving$effective_sample_size, na.rm = TRUE)
    } else NA_real_,
    nonconstant_chains_with_ess_at_least_100 = if (nrow(moving)) {
      sum(moving$effective_sample_size >= 100, na.rm = TRUE)
    } else 0L,
    median_pip_monte_carlo_se = if (nrow(moving)) {
      median(moving$pip_monte_carlo_se, na.rm = TRUE)
    } else NA_real_,
    maximum_pip_monte_carlo_se = if (nrow(moving)) {
      max(moving$pip_monte_carlo_se, na.rm = TRUE)
    } else NA_real_
  )
}))
rownames(gamma_summary) <- NULL
gamma_summary <- gamma_summary[order(
  gamma_summary$molecular_layer, gamma_summary$risk_group
), ]
write.csv(
  gamma_summary,
  file.path(out_dir, "saved_fit_pip_mcmc_diagnostics_summary.csv"),
  row.names = FALSE
)

cat("Saved-fit single-chain diagnostics complete.\n")
cat("These diagnostics do not establish convergence; multi-chain R-hat is unavailable.\n\n")
print(diagnostic_summary, row.names = FALSE)
cat("\nInclusion-indicator diagnostics (constant chains have undefined ESS):\n")
print(gamma_summary, row.names = FALSE)
