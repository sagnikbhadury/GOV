# Create a reviewer-facing trace-plot comparison for the NFkB pilot.

options(stringsAsFactors = FALSE)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
pilot_dir <- file.path(root, "reviewer_reanalysis", "tuned_pilot")
summary_path <- file.path(pilot_dir, "tuned_refit_final_summary.csv")
if (!file.exists(summary_path)) stop("Run the tuning grid and final pilot first.")

tab <- read.csv(summary_path, check.names = FALSE)
step_a <- tab$step_a[[1L]]
step_u <- tab$step_u[[1L]]

tuned_path <- file.path(
  pilot_dir,
  sprintf(
    "HR_Pathway_NFkB_stepa%s_stepu%s_final.rds",
    format(step_a, trim = TRUE), format(step_u, trim = TRUE)
  )
)
original_path <- file.path(
  root, "luad", "HR_luad_lsGPRVS_result1000_NFkB_activity.rds"
)
original <- readRDS(original_path)
tuned <- readRDS(tuned_path)

post_draws <- function(fit) {
  keep <- seq.int(fit$burn + 1L, fit$nrun)
  list(
    iteration = keep,
    psia = fit$PSIA[keep],
    psiu = fit$PSIU[keep],
    loglik = fit$post.loglik[keep],
    selected = colSums(fit$GAMMA[, keep, drop = FALSE])
  )
}
old <- post_draws(original)
new <- post_draws(tuned)

draw_panel <- function(draws, key, label, color) {
  plot(
    draws$iteration, draws[[key]], type = "l", col = color, lwd = 0.7,
    xlab = "Post-burn-in iteration", ylab = label
  )
}

draw_figure <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(
    mfrow = c(4, 2), mar = c(4.4, 4.2, 2.3, 0.8),
    oma = c(0.3, 0, 2, 0), mgp = c(2.5, 0.75, 0)
  )
  keys <- c("psia", "psiu", "loglik", "selected")
  labels <- c(
    expression(psi[a]), expression(psi[u]),
    "Log likelihood", "Selected nodes"
  )
  for (i in seq_along(keys)) {
    draw_panel(old, keys[[i]], labels[[i]], "#666666")
    if (i == 1L) title("(A) Original: step 0.01")
    draw_panel(new, keys[[i]], labels[[i]], "#0072B2")
    if (i == 1L) title(sprintf("(B) Tuned: steps %.2f / %.2f", step_a, step_u))
  }
  mtext(
    "High-risk NFkB saved-latent-space fit: original and tuned traces",
    outer = TRUE, cex = 1.05
  )
}

pdf_path <- file.path(root, "manuscript", "supplementary_traceplot_NFkB.pdf")
pdf(pdf_path, width = 8.5, height = 11)
draw_figure()
dev.off()

png_path <- file.path(root, "manuscript", "supplementary_traceplot_NFkB.png")
png(png_path, width = 1700, height = 2200, res = 200)
draw_figure()
dev.off()

cat("Selected steps (psi_a, psi_u):", step_a, step_u, "\n")
cat("Wrote:", pdf_path, "\n")
cat("Wrote:", png_path, "\n")
