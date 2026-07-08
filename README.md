# GOV: Graph on Vector regression with latent scale GP

An R package for regressing a scalar response on a network-valued
predictor:

* `EM_toU()` / `EM_toU_ver0()` -- fit a latent space network model to an
  observed adjacency matrix via EM, recovering low-rank latent node
  positions.
* `lsGPR_VS()` -- regress a scalar response on a Gaussian process kernel
  built from those latent positions (and optionally other covariates),
  with a spike-and-slab prior that selects which nodes' positions actually
  drive the kernel.

## Installation

```r
# install.packages("remotes")
remotes::install_github("sagnikbhadury/GOV")
```

## Usage

```r
library(GOV)

fit_latent <- EM_toU(G = adjacency_matrix, d = 3)

fit <- lsGPR_VS(a = a_train, a.te = a_test,
                U = U_train_list, U.te = U_test_list,
                Y = y_train,
                burn = 5000, nrun = 10000)

# posterior mean predictions on held-out units
pred_mean <- rowMeans(fit$post.pred[, (fit$burn + 1):fit$nrun])
```

See `?EM_toU` and `?lsGPR_VS` for full documentation.
