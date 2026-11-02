# Regenerate the quantitative PIP heatmaps from iterations 501--1000.

options(stringsAsFactors = FALSE)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
dat <- read.csv(
  file.path(root, "reviewer_reanalysis", "results", "post_burnin_pips.csv"),
  check.names = FALSE
)

layers <- c("Pathway", "Kinase", "Transcription factor")
risks <- c("HR", "LR")

make_matrix <- function(layer, risk) {
  x <- dat[dat$molecular_layer == layer & dat$risk_group == risk, ]
  hpcs <- unique(x$hpc)
  activities <- unique(x$activity)
  mat <- matrix(
    NA_real_, nrow = length(hpcs), ncol = length(activities),
    dimnames = list(hpcs, activities)
  )
  for (i in seq_len(nrow(x))) mat[x$hpc[[i]], x$activity[[i]]] <- x$pip[[i]]
  mat
}

draw_heatmap <- function(mat, main) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  image(
    x = seq_len(nc), y = seq_len(nr), z = t(mat),
    zlim = c(0, 1), col = hcl.colors(100, "Blue-Red 3"),
    axes = FALSE, xlab = "", ylab = "", main = main, cex.main = 1.15
  )
  axis(1, at = seq_len(nc), labels = colnames(mat), las = 2, cex.axis = 1.15)
  axis(2, at = seq_len(nr), labels = rownames(mat), las = 2, cex.axis = 1.15)
  selected <- which(mat > 0.5, arr.ind = TRUE)
  if (nrow(selected)) {
    text(selected[, "col"], selected[, "row"], labels = "*", cex = 0.55)
  }
  box()
}

draw_figure <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  layout(matrix(seq_len(6), nrow = 3, byrow = TRUE))
  par(mar = c(8.3, 6.0, 2.3, 1.0), oma = c(0, 0, 2.2, 0))
  for (layer in layers) {
    for (risk in risks) {
      draw_heatmap(
        make_matrix(layer, risk),
        paste(layer, if (risk == "HR") "high-risk" else "low-risk")
      )
    }
  }
  mtext(
    "Post-burn-in posterior inclusion probabilities (iterations 501--1000)",
    outer = TRUE, cex = 1.1
  )
}

pdf_path <- file.path(root, "manuscript", "postburnin_pip_heatmaps.pdf")
pdf(pdf_path, width = 13, height = 13)
draw_figure()
dev.off()

png_path <- file.path(root, "manuscript", "postburnin_pip_heatmaps.png")
png(png_path, width = 2600, height = 2600, res = 200)
draw_figure()
dev.off()

cat("Wrote", pdf_path, "and", png_path, "\n")
