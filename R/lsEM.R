## Stage 1 calculations for an EM algorithm to estimate latent positions in a
## network model (latent space model with a global shift + low-rank latent
## coordinates U). Includes utilities to convert between vector and matrix
## representations of lower-triangular matrices.

#' Extract the strict lower triangle of a matrix as a vector
#'
#' @param mat A square matrix.
#' @return A numeric vector of the strict lower-triangular entries of
#'   \code{mat}, in column-major order.
#' @export
func.mtv <- function(mat) {
  lower <- mat
  lower[upper.tri(lower, diag = TRUE)] <- NA
  c.lower <- c(lower)
  c.lower[!is.na(c.lower)]
}

#' Reconstruct a (optionally symmetric) matrix from its lower-triangular vector
#'
#' More flexible version of \code{\link{func.vtm}}: allows specifying the
#' target matrix size explicitly, custom fill values for the upper triangle
#' and diagonal, and optional input validation.
#'
#' @param vec Vector of lower-triangular entries.
#' @param matrix_size Optional target matrix dimension \code{n}; if omitted,
#'   inferred from \code{length(vec)}.
#' @param fill_upper Value to initially fill the upper triangle with, before
#'   symmetrizing (ignored if \code{symmetric = TRUE} other than as a base
#'   value). Default 0.
#' @param fill_diag Value to fill the diagonal with. Default 0.
#' @param symmetric If \code{TRUE} (default), mirror the lower triangle into
#'   the upper triangle.
#' @param validate If \code{TRUE} (default), check that \code{vec} has a
#'   length consistent with a valid triangular number (or with
#'   \code{matrix_size} if supplied), padding/truncating with a warning if not.
#'
#' @return An \code{n x n} matrix.
#' @export
vector_to_lower_tri <- function(vec, matrix_size = NULL,
                                 fill_upper = 0, fill_diag = 0,
                                 symmetric = TRUE, validate = TRUE) {
  k <- length(vec)

  if (!is.null(matrix_size)) {
    # Use provided matrix size
    n <- matrix_size
    expected_k <- n * (n - 1) / 2

    if (k != expected_k && validate) {
      stop(paste("Vector length", k, "doesn't match expected", expected_k,
                 "for", n, "x", n, "matrix"))
    }

    # Pad or truncate vector if needed
    if (k < expected_k) {
      vec <- c(vec, rep(0, expected_k - k))
      warning(paste("Vector padded with", expected_k - k, "zeros"))
    } else if (k > expected_k) {
      vec <- vec[1:expected_k]
      warning(paste("Vector truncated to", expected_k, "elements"))
    }

  } else {
    # Calculate matrix size from vector length
    n <- (1 + sqrt(1 + 8 * k)) / 2

    if (abs(n - round(n)) > 1e-10 && validate) {
      valid_lengths <- sapply(2:20, function(x) x * (x - 1) / 2)
      stop(paste("Invalid vector length:", k,
                 "\nValid lengths:", paste(valid_lengths, collapse = ", ")))
    }

    n <- round(n)
  }

  # Create matrix
  result_matrix <- matrix(fill_upper, n, n)
  diag(result_matrix) <- fill_diag

  # Fill lower triangle
  result_matrix[lower.tri(result_matrix)] <- vec

  # Make symmetric if requested
  if (symmetric) {
    result_matrix[upper.tri(result_matrix)] <- t(result_matrix)[upper.tri(result_matrix)]
  }

  result_matrix
}

#' Reconstruct a symmetric matrix from its lower-triangular vector
#'
#' Inverse of \code{\link{func.mtv}}: given the strict lower-triangular
#' entries of a symmetric matrix (diagonal assumed zero), reconstructs the
#' full symmetric matrix.
#'
#' @param vec Vector of lower-triangular entries, of length \code{n*(n-1)/2}
#'   for some integer \code{n}.
#' @return The reconstructed \code{n x n} symmetric matrix (zero diagonal).
#' @export
func.vtm <- function(vec) {
  v <- length(c(vec))
  n <- (1 + sqrt(8 * v + 1)) / 2
  lower <- matrix(rep(0, n^2), n)
  lower[lower.tri(lower, diag = FALSE)] <- vec
  lower + t(lower)
}

#' EM algorithm for a latent space network model (regularized, with intercept)
#'
#' Estimates a global intercept \code{a} and low-rank latent position matrix
#' \code{U} (with a fixed leading column, and \code{d-1} free coordinates per
#' node) for a network model where edge probabilities depend on the inner
#' product of nodes' latent positions, via an EM algorithm with Gaussian
#' priors on \code{a} and the latent coordinates. Latent coordinates are
#' initialized from classical scaling of graph-distance.
#'
#' @param G A symmetric adjacency matrix (0/1), or the vector of its strict
#'   lower triangle (as produced by \code{\link{func.mtv}}).
#' @param d Latent space dimension (including the fixed leading column).
#'   Default 3, giving 2 free coordinates per node.
#' @param u1 Fixed value of the leading latent coordinate for every node.
#'   Default 0.5.
#' @param sig.a Prior variance for the intercept \code{a}. Default 2.
#' @param sig.u Prior variance for the free latent coordinates. Default 2.
#' @param tol Convergence tolerance on the max absolute parameter change.
#'   Default 1e-10.
#' @param maxiter Maximum number of EM iterations. Default 200.
#'
#' @return A list with \code{U} (the \code{n x d} latent position matrix,
#'   first column fixed at \code{u1}) and \code{a} (the estimated intercept).
#' @export
EM_toU <- function(G, d = 3, u1 = 0.5, sig.a = 2, sig.u = 2, tol = 10^(-10), maxiter = 200) {
  G <- c(G[lower.tri(G, diag = FALSE)]); v <- length(G); p <- (1 + sqrt(8 * v + 1)) / 2
  sG <- func.vtm(vec = G)
  g <- igraph::graph_from_adjacency_matrix(sG)
  D <- igraph::distances(g, mode = 'all')
  D[D == Inf] <- max(D[D != Inf]) + igraph::mean_distance(g)
  z0 <- D[, 1:(d - 1), drop = FALSE]
  # Normalize
  z0 <- scale(z0)

  u1 <- 0.5
  U <- cbind(u1, matrix(z0, ncol = (d - 1)))
  a <- 0
  step <- 0; diff <- 1
  while (step <= maxiter & diff > tol) {
    # E step
    dotU <- func.mtv(U %*% t(U))
    phi <- a + dotU
    W <- (exp(phi) - 1) / (2 * phi * (1 + exp(phi)))
    mW <- func.vtm(vec = W)

    # M step
    # max a
    old.a <- a
    a <- sum(G - W * dotU - 0.5) / (1 / sig.a + sum(W))

    old.U <- U

    # max U[k,-1]
    for (k in 1:p) {
      tmp1 <- t(U[-k, -1]) %*% diag(mW[-k, k]) %*% U[-k, -1] + (1 / sig.u) * diag(rep(1, (d - 1)))
      tmp2 <- matrix(sG[-k, k] - mW[-k, k] * (a + u1^2) - 0.5, nrow = 1) %*% U[-k, -1]
      U[k, -1] <- tmp2 %*% MASS::ginv(tmp1)
    }

    diff.U <- max(abs(c(old.U - U)))
    diff.a <- abs(old.a - a)
    diff <- max(diff.U, diff.a)
    step <- step + 1
    message("iteration: ", step)
  }
  list(U = U, a = a)
}

#' EM algorithm for a latent space network model (unregularized baseline)
#'
#' Earlier, unregularized version of \code{\link{EM_toU}}: no priors on
#' \code{a} or the latent coordinates, and latent coordinates are
#' initialized from classical multidimensional scaling
#' (\code{\link[stats]{cmdscale}}) of graph distance rather than the raw
#' distance columns.
#'
#' @param G A symmetric adjacency matrix (0/1).
#' @param d Latent space dimension (including the fixed leading column).
#'   Default 3.
#' @param tol Convergence tolerance on the max absolute parameter change.
#'   Default 1e-10.
#' @param maxiter Maximum number of EM iterations. Default 200.
#'
#' @return A list with \code{U} (the \code{n x d} latent position matrix)
#'   and \code{a} (the estimated intercept).
#' @export
EM_toU_ver0 <- function(G, d = 3, tol = 10^(-10), maxiter = 200) {
  G <- c(G); v <- length(G); p <- (1 + sqrt(8 * v + 1)) / 2
  sG <- func.vtm(G)
  g <- igraph::graph_from_adjacency_matrix(sG)
  D <- igraph::distances(g, mode = 'all')
  D[D == Inf] <- max(D[D != Inf]) + 1
  z0 <- stats::cmdscale(D, (d - 1))

  u1 <- 0.5
  U <- cbind(u1, matrix(z0, ncol = (d - 1)))
  a <- -log(1 / mean(G) - 1)
  step <- 0; diff <- 1
  while (step <= maxiter & diff > tol) {
    # E step
    dotU <- func.mtv(U %*% t(U))
    phi <- a + dotU
    W <- (exp(phi) - 1) / (2 * phi * (1 + exp(phi)))
    mW <- func.vtm(W)

    # M step
    # max a
    old.a <- a
    a <- sum(G - W * dotU - 0.5) / sum(W)

    old.U <- U

    # max U[k,-1]
    for (k in 1:p) {
      tmp1 <- t(U[-k, -1]) %*% diag(mW[-k, k]) %*% U[-k, -1] + 0.5 * diag(rep(1, (d - 1)))
      tmp2 <- matrix(sG[-k, k] - mW[-k, k] * (a + u1^2) - 0.5, nrow = 1) %*% U[-k, -1]
      U[k, -1] <- tmp2 %*% solve(tmp1)
    }

    diff.U <- norm(old.U - U, "F")
    diff.a <- abs(old.a - a)
    diff <- max(diff.U, diff.a)
    step <- step + 1
  }
  list(U = U, a = a)
}
