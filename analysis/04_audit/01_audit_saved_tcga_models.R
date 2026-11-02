# Reproducible audit of the saved TCGA-LUAD graph-on-vector model objects.
#
# This script deliberately uses only the original TCGA inputs and saved MCMC
# fits.  It does not use the exploratory HEST-1k files.  It reconstructs PIPs
# from post-burn-in draws, evaluates PIP-threshold sensitivity, reports the
# posterior expected false-discovery proportion (rather than treating 1-PIP as
# a frequentist p-value), and evaluates the ten observations held out by each
# original fitting script.

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

adjacency_example <- read.csv(
  list.files(file.path(root, "adj_matricesV3"), full.names = TRUE)[[1]],
  check.names = FALSE
)
hpc_names <- colnames(adjacency_example)[-1L]
stopifnot(length(hpc_names) == 20L)

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

load_response_table <- local({
  cache <- new.env(parent = emptyenv())
  function(risk, layer) {
    key <- paste(risk, layer, sep = "|")
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    risk_name <- if (risk == "HR") "high" else "low"
    path <- switch(
      layer,
      "Pathway" = file.path(
        root, "TCGA-LUAD-data",
        sprintf("LUAD_%s_risk_activity_scores.rds", risk_name)
      ),
      "Kinase" = file.path(
        root, "TCGA-LUAD-data",
        sprintf("top_kinase_df_%srisk_luad.rds", risk_name)
      ),
      "Transcription factor" = file.path(
        root, "TCGA-LUAD-data",
        sprintf("top_tf_df_%srisk_luad.rds", risk_name)
      )
    )
    ans <- readRDS(path)
    assign(key, ans, envir = cache)
    ans
  }
})

load_patient_ids <- function(risk) {
  readRDS(file.path(
    root, "TCGA-LUAD-data", sprintf("luad_%s_PATIDs.rds", risk)
  ))
}

pip_rows <- list()
prediction_rows <- list()

for (i in seq_along(result_files)) {
  path <- result_files[[i]]
  info <- parse_fit_name(path)
  fit <- readRDS(path)
  keep <- seq.int(fit$burn + 1L, fit$nrun)

  pips <- rowMeans(fit$GAMMA[, keep, drop = FALSE])
  pip_rows[[i]] <- data.frame(
    risk_group = info$risk,
    molecular_layer = info$layer,
    activity = info$feature,
    hpc_index = seq_along(pips),
    hpc = hpc_names,
    pip = pips
  )

  response <- load_response_table(info$risk, info$layer)
  activity_response <- response[response$source == info$feature, , drop = FALSE]
  patient_ids <- load_patient_ids(info$risk)
  score_index <- match(patient_ids, activity_response$condition)
  if (anyNA(score_index)) {
    stop("Could not align response for ", basename(path))
  }
  y <- activity_response$score[score_index]
  y_test <- y[seq_len(10L)]
  y_train <- y[-seq_len(10L)]

  pred_draws <- fit$post.pred[, keep, drop = FALSE]
  pred_mean <- rowMeans(pred_draws)
  pred_lower <- apply(pred_draws, 1L, quantile, probs = 0.025, names = FALSE)
  pred_upper <- apply(pred_draws, 1L, quantile, probs = 0.975, names = FALSE)
  baseline <- rep(mean(y_train), length(y_test))
  sse <- sum((y_test - pred_mean)^2)
  sst <- sum((y_test - mean(y_test))^2)

  prediction_rows[[i]] <- data.frame(
    risk_group = info$risk,
    molecular_layer = info$layer,
    activity = info$feature,
    n_train = length(y_train),
    n_test = length(y_test),
    rmse = sqrt(mean((y_test - pred_mean)^2)),
    mae = mean(abs(y_test - pred_mean)),
    test_r2 = if (sst > 0) 1 - sse / sst else NA_real_,
    correlation = if (sd(y_test) > 0 && sd(pred_mean) > 0) {
      cor(y_test, pred_mean)
    } else {
      NA_real_
    },
    interval_coverage_95 = mean(y_test >= pred_lower & y_test <= pred_upper),
    baseline_rmse = sqrt(mean((y_test - baseline)^2)),
    rmse_improvement_percent = 100 * (
      sqrt(mean((y_test - baseline)^2)) - sqrt(mean((y_test - pred_mean)^2))
    ) / sqrt(mean((y_test - baseline)^2))
  )
}

pips <- do.call(rbind, pip_rows)
predictions <- do.call(rbind, prediction_rows)

write.csv(pips, file.path(out_dir, "post_burnin_pips.csv"), row.names = FALSE)
write.csv(
  predictions,
  file.path(out_dir, "heldout_prediction_metrics_by_activity.csv"),
  row.names = FALSE
)

thresholds <- c(0.3, 0.5, 0.7, 0.9)
threshold_summary <- do.call(rbind, lapply(
  split(pips, interaction(pips$risk_group, pips$molecular_layer, drop = TRUE)),
  function(dat) do.call(rbind, lapply(thresholds, function(threshold) {
    data.frame(
      risk_group = dat$risk_group[[1]],
      molecular_layer = dat$molecular_layer[[1]],
      threshold = threshold,
      selected = sum(dat$pip > threshold),
      total_tests = nrow(dat)
    )
  }))
))
rownames(threshold_summary) <- NULL
write.csv(
  threshold_summary,
  file.path(out_dir, "pip_threshold_sensitivity.csv"),
  row.names = FALSE
)

# Preserve a separate sensitivity summary for the PIPs already exported and
# used in the submitted figures/manuscript.  This lets the revision answer the
# threshold question without changing the manuscript's core TCGA estimates.
published_files <- c(
  Pathway = "luad_Activity_lsGPRVS_pip_longformat.csv",
  Kinase = "luad_kinase_lsGPRVS_pip_longformat.csv",
  `Transcription factor` = "luad_tf_lsGPRVS_pip_longformat.csv"
)
published_threshold_rows <- list()
k <- 0L
for (layer in names(published_files)) {
  dat <- read.csv(
    file.path(root, "luad", published_files[[layer]]),
    check.names = FALSE
  )
  for (risk in c("high", "low")) {
    for (threshold in thresholds) {
      k <- k + 1L
      published_threshold_rows[[k]] <- data.frame(
        risk_group = risk,
        molecular_layer = layer,
        threshold = threshold,
        selected = sum(dat$PIP[dat$risk == risk] > threshold),
        total_tests = sum(dat$risk == risk)
      )
    }
  }
}
published_threshold_summary <- do.call(rbind, published_threshold_rows)
write.csv(
  published_threshold_summary,
  file.path(out_dir, "published_pip_threshold_sensitivity.csv"),
  row.names = FALSE
)

# For a selected set S, mean(1-PIP) is the posterior expected false-discovery
# proportion.  Within each molecular layer, choose the largest prefix of PIPs
# sorted from high to low whose expected FDP is <= 0.05.
posterior_fdr_rows <- lapply(split(pips, pips$molecular_layer), function(dat) {
  dat <- dat[order(dat$pip, decreasing = TRUE), , drop = FALSE]
  dat$posterior_expected_fdp <- cumsum(1 - dat$pip) / seq_len(nrow(dat))
  admissible <- which(dat$posterior_expected_fdp <= 0.05)
  n_selected <- if (length(admissible)) max(admissible) else 0L
  dat$selected_at_5pct <- seq_len(nrow(dat)) <= n_selected
  dat
})
posterior_fdr <- do.call(rbind, posterior_fdr_rows)
rownames(posterior_fdr) <- NULL
write.csv(
  posterior_fdr,
  file.path(out_dir, "posterior_expected_fdr.csv"),
  row.names = FALSE
)

prediction_summary <- do.call(rbind, lapply(
  split(predictions, interaction(
    predictions$risk_group, predictions$molecular_layer, drop = TRUE
  )),
  function(dat) data.frame(
    risk_group = dat$risk_group[[1]],
    molecular_layer = dat$molecular_layer[[1]],
    models = nrow(dat),
    median_rmse = median(dat$rmse),
    median_test_r2 = median(dat$test_r2, na.rm = TRUE),
    median_correlation = median(dat$correlation, na.rm = TRUE),
    median_coverage_95 = median(dat$interval_coverage_95),
    models_beating_mean_baseline = sum(dat$rmse < dat$baseline_rmse),
    median_rmse_improvement_percent = median(dat$rmse_improvement_percent)
  )
))
rownames(prediction_summary) <- NULL
write.csv(
  prediction_summary,
  file.path(out_dir, "heldout_prediction_summary.csv"),
  row.names = FALSE
)

cat("Saved-model audit complete.\n")
cat("Post-burn-in PIP counts at 0.5:\n")
print(threshold_summary[threshold_summary$threshold == 0.5, ])
cat("\nSelections at posterior expected FDP <= 0.05:\n")
print(with(posterior_fdr, table(molecular_layer, selected_at_5pct)))
cat("\nHeld-out prediction summary:\n")
print(prediction_summary)
cat("\nOutputs written to:", out_dir, "\n")
