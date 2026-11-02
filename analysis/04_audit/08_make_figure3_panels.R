# Build Figure 3 for submission: original vs tuned NFkB trace plots with panel
# letters a)-h) and no in-image figure title (npj Artificial Intelligence
# technical checks: figure titles/labels must not appear in the image file).

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
orig <- post_draws(original)
new <- post_draws(tuned)

draw_panel <- function(draws, key, label, color, letter, header) {
  plot(
    draws$iteration, draws[[key]], type = "l", col = color, lwd = 0.7,
    xlab = "Post-burn-in iteration", ylab = label, main = ""
  )
  mtext(letter, side = 3, line = 0.4, adj = 0, cex = 0.95, font = 2)
  mtext(header, side = 3, line = 0.4, adj = 0.5, cex = 0.72)
}

draw_figure <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(
    mfrow = c(4, 2), mar = c(4.4, 4.2, 2.6, 0.8),
    oma = c(0.3, 0, 0.6, 0), mgp = c(2.5, 0.75, 0)
  )
  keys <- c("psia", "psiu", "loglik", "selected")
  labels <- c(
    expression(psi[a]), expression(psi[u]),
    "Log likelihood", "Selected nodes"
  )
  head_orig <- "Original, proposal SD 0.01"
  head_tuned <- sprintf("Tuned, proposal SD %.2f / %.2f", step_a, step_u)
  letters8 <- letters[1:8]
  k <- 1L
  for (i in seq_along(keys)) {
    draw_panel(orig, keys[[i]], labels[[i]], "#666666", letters8[[k]], head_orig)
    k <- k + 1L
    draw_panel(new, keys[[i]], labels[[i]], "#0072B2", letters8[[k]], head_tuned)
    k <- k + 1L
  }
}

pdf_path <- file.path(root, "manuscript", "Figure3.pdf")
pdf(pdf_path, width = 8.5, height = 11)
draw_figure()
dev.off()

png_path <- file.path(root, "manuscript", "Figure3.png")
png(png_path, width = 1700, height = 2200, res = 200)
draw_figure()
dev.off()

cat("Steps (psi_a, psi_u):", step_a, step_u, "\n")
cat("Wrote:", pdf_path, "and", png_path, "\n")
