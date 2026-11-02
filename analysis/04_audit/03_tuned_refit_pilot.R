# Pilot proposal tuning for selected saved-latent-space TCGA-LUAD fits.
# Run from the repository root. Optional command-line arguments:
#   Rscript reviewer_reanalysis/03_tuned_refit_pilot.R grid
#   Rscript reviewer_reanalysis/03_tuned_refit_pilot.R final

options(stringsAsFactors = FALSE)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(root, "TCGA-LUAD-data"))) {
  stop("Run this script from the repository root.")
}

personal_library <- file.path(
  "C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/R_packages"
)
if (dir.exists(personal_library)) .libPaths(c(personal_library, .libPaths()))

suppressPackageStartupMessages({
  library(matrixcalc)
  library(mvtnorm)
  library(MCMCpack)
})
source(file.path(root, "LSGPR", "lsGPR_VS.R"))

out_dir <- file.path(root, "reviewer_reanalysis", "tuned_pilot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ess_ips <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 3L || !is.finite(sd(x)) || sd(x) == 0) return(NA_real_)
  rho <- as.numeric(acf(
    x, plot = FALSE, lag.max = min(n - 1L, floor(n / 2)),
    type = "correlation", demean = TRUE
  )$acf)[-1L]
  pair_count <- floor(length(rho) / 2L)
  if (!pair_count) return(as.numeric(n))
  pair_sums <- rho[2L * seq_len(pair_count) - 1L] +
    rho[2L * seq_len(pair_count)]
  stop_at <- which(!is.finite(pair_sums) | pair_sums <= 0)
  keep_pairs <- if (length(stop_at)) stop_at[[1L]] - 1L else pair_count
  tau <- if (keep_pairs > 0L) {
    1 + 2 * sum(rho[seq_len(2L * keep_pairs)])
  } else 1
  min(as.numeric(n), max(1, n / tau))
}

load_case <- function(risk, layer, activity_name) {
  risk_word <- if (risk == "HR") "high" else "low"
  response_path <- switch(
    layer,
    Pathway = file.path(
      root, "TCGA-LUAD-data",
      sprintf("LUAD_%s_risk_activity_scores.rds", risk_word)
    ),
    Kinase = file.path(
      root, "TCGA-LUAD-data",
      sprintf("top_kinase_df_%srisk_luad.rds", risk_word)
    ),
    TF = file.path(
      root, "TCGA-LUAD-data", sprintf("top_tf_df_%srisk_luad.rds", risk_word)
    )
  )
  response <- readRDS(response_path)
  response <- response[response$source == activity_name, , drop = FALSE]
  ids <- readRDS(file.path(
    root, "TCGA-LUAD-data", sprintf("luad_%s_PATIDs.rds", risk)
  ))
  index <- match(ids, response$condition)
  if (anyNA(index)) stop("Response alignment failed for ", activity_name)
  y <- as.numeric(response$score[index])
  a_all <- readRDS(file.path(
    root, "TCGA-LUAD-data", sprintf("luad_%s_mya.rds", risk)
  ))
  u_all <- readRDS(file.path(
    root, "TCGA-LUAD-data", sprintf("luad_%s_myU.rds", risk)
  ))
  test <- seq_len(10L)
  list(
    a = a_all[-test], a_test = a_all[test],
    U = u_all[-test], U_test = u_all[test],
    Y = y[-test], Y_test = y[test]
  )
}

fit_one <- function(spec, step_a, step_u, burn, nrun, seed) {
  dat <- load_case(spec$risk, spec$layer, spec$activity)
  set.seed(seed)
  elapsed <- system.time({
    fit <- lsGPR_VS(
      dat$a, dat$a_test, dat$U, dat$U_test,
      z = NULL, z.te = NULL, Y = dat$Y,
      burn = burn, nrun = nrun,
      prior = list(
        a.t = 15, b.t = 1, a.p1 = 2, b.p1 = 50,
        step = step_a, a.p = 1, b.p = 10
      ),
      control = list(
        step.a = step_a, step.u = step_u,
        # Keep the submitted node-update kernel unchanged; this pilot tunes
        # only the length-scale proposal and therefore does not alter the
        # inferential target or replace the primary analysis.
        correct.within = FALSE, show.progress = FALSE
      )
    )
  })[["elapsed"]]
  keep <- seq.int(burn + 1L, nrun)
  old_fit <- readRDS(file.path(
    root, "luad",
    if (spec$layer == "Pathway") {
      sprintf(
        "%s_luad_lsGPRVS_result1000_%s_activity.rds",
        spec$risk, spec$activity
      )
    } else if (spec$layer == "Kinase") {
      sprintf(
        "%s_luad_kinase_lsGPRVS_result1000_%s_activity.rds",
        spec$risk, spec$activity
      )
    } else {
      sprintf(
        "%s_luad_TF_lsGPRVS_result1000_%s_activity.rds",
        spec$risk, spec$activity
      )
    }
  ))
  old_keep <- seq.int(old_fit$burn + 1L, old_fit$nrun)
  old_pip <- rowMeans(old_fit$GAMMA[, old_keep, drop = FALSE])
  new_pip <- rowMeans(fit$GAMMA[, keep, drop = FALSE])
  pred <- rowMeans(fit$post.pred[, keep, drop = FALSE])
  baseline <- rep(mean(dat$Y), length(dat$Y_test))

  summary <- data.frame(
    risk_group = spec$risk,
    molecular_layer = spec$layer,
    activity = spec$activity,
    step_a = step_a,
    step_u = step_u,
    burn = burn,
    nrun = nrun,
    elapsed_seconds = elapsed,
    acceptance_psia = mean(fit$ACCEPT_PSIA[keep]),
    acceptance_psiu = mean(fit$ACCEPT_PSIU[keep]),
    ess_psia = ess_ips(fit$PSIA[keep]),
    ess_psiu = ess_ips(fit$PSIU[keep]),
    ess_log_likelihood = ess_ips(fit$post.loglik[keep]),
    pip_correlation_with_original = cor(old_pip, new_pip),
    selections_original = sum(old_pip > 0.5),
    selections_tuned = sum(new_pip > 0.5),
    median_absolute_pip_change = median(abs(new_pip - old_pip)),
    maximum_absolute_pip_change = max(abs(new_pip - old_pip)),
    heldout_rmse = sqrt(mean((dat$Y_test - pred)^2)),
    baseline_rmse = sqrt(mean((dat$Y_test - baseline)^2))
  )
  list(fit = fit, summary = summary)
}

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "grid"

if (mode == "grid") {
  specs <- list(list(risk = "HR", layer = "Pathway", activity = "NFkB"))
  proposal_pairs <- list(c(0.10, 0.10), c(0.30, 0.30), c(1.00, 1.00))
  burn <- 150L
  nrun <- 300L
} else if (mode == "final") {
  specs <- list(list(risk = "HR", layer = "Pathway", activity = "NFkB"))
  # The short screen showed that psi_a requires a larger move than psi_u.
  # These separate values interpolate/extrapolate to the usual 20--60%
  # acceptance window without altering any other update in the submitted fit.
  proposal_pairs <- list(c(1.50, 0.50))
  burn <- 500L
  nrun <- 1500L
} else {
  stop("Mode must be 'grid' or 'final'.")
}

summaries <- list()
k <- 0L
for (spec in specs) {
  for (proposal_pair in proposal_pairs) {
    k <- k + 1L
    step_a <- proposal_pair[[1L]]
    step_u <- proposal_pair[[2L]]
    cat(
      "Fitting", spec$risk, spec$layer, spec$activity,
      "steps", step_a, step_u, "\n"
    )
    ans <- fit_one(spec, step_a, step_u, burn, nrun, seed = 2020L + k)
    fit_name <- sprintf(
      "%s_%s_%s_stepa%s_stepu%s_%s.rds",
      spec$risk, spec$layer, spec$activity,
      format(step_a, trim = TRUE), format(step_u, trim = TRUE), mode
    )
    saveRDS(ans$fit, file.path(out_dir, fit_name))
    summaries[[k]] <- ans$summary
    print(ans$summary)
  }
}

summary_table <- do.call(rbind, summaries)
write.csv(
  summary_table,
  file.path(out_dir, sprintf("tuned_refit_%s_summary.csv", mode)),
  row.names = FALSE
)
