# Build Figure 2 for submission: post-burn-in PIP heatmaps with panel letters
# a)-f) and no in-image figure title (npj Artificial Intelligence technical
# checks: figure titles/labels must not appear in the image file).

options(stringsAsFactors = FALSE)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
dat <- read.csv(
  file.path(root, "reviewer_reanalysis", "results", "post_burnin_pips.csv"),
  check.names = FALSE
)

layers <- c("Pathway", "Kinase", "Transcription factor")
risks <- c("HR", "LR")
panel_letters <- letters[1:6]

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

draw_heatmap <- function(mat, letter, main) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  image(
    x = seq_len(nc), y = seq_len(nr), z = t(mat),
    zlim = c(0, 1), col = hcl.colors(100, "Blue-Red 3"),
    axes = FALSE, xlab = "", ylab = "", main = "", cex.main = 1.15
  )
  axis(1, at = seq_len(nc), labels = colnames(mat), las = 2, cex.axis = 1.15)
  axis(2, at = seq_len(nr), labels = rownames(mat), las = 2, cex.axis = 1.15)
  selected <- which(mat > 0.5, arr.ind = TRUE)
  if (nrow(selected)) {
    text(selected[, "col"], selected[, "row"], labels = "*", cex = 0.55)
  }
  box()
  mtext(letter, side = 3, line = 0.45, adj = 0, cex = 1.15, font = 2)
  mtext(main, side = 3, line = 0.45, adj = 0.5, cex = 0.95)
}

draw_figure <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  layout(matrix(seq_len(6), nrow = 3, byrow = TRUE))
  par(mar = c(8.3, 6.0, 2.6, 1.0), oma = c(0, 0, 0.6, 0))
  k <- 1L
  for (layer in layers) {
    for (risk in risks) {
      draw_heatmap(
        make_matrix(layer, risk),
        panel_letters[[k]],
        paste(layer, if (risk == "HR") "high-risk" else "low-risk")
      )
      k <- k + 1L
    }
  }
}

pdf_path <- file.path(root, "manuscript", "Figure2.pdf")
pdf(pdf_path, width = 13, height = 13)
draw_figure()
dev.off()

png_path <- file.path(root, "manuscript", "Figure2.png")
png(png_path, width = 2600, height = 2600, res = 200)
draw_figure()
dev.off()

cat("Wrote", pdf_path, "and", png_path, "\n")
