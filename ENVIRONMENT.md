# Software environment

Two R sessions produced the published results. Versions below are taken from
the OmnipathR session log written during the activity-inference run
(`analysis/01_data_preparation/omnipathr-20251115-1340_activity_inference_session.log`,
15 Nov 2025) and from the library that ran the graph-on-vector fits and the
audit. Where a version could not be read from a log or a saved object it is
marked as such.

## Activity inference (15 Nov 2025)

| Component | Version |
|---|---|
| R | 4.5.1 (2025-06-13 ucrt), Windows 11 x64 |
| decoupleR | 2.14.0 (2025-04-15), Bioconductor 3.21 |
| OmnipathR | 3.16.2 (2025-09-03), Bioconductor 3.21 |
| Seurat / SeuratObject | 5.3.1 / 5.2.0 |
| dplyr / tidyr / tibble | 1.1.4 / 1.3.1 / 3.3.0 |
| ggplot2 | 4.0.0 |
| igraph | 2.2.1 |

Prior-knowledge networks, all retrieved through OmnipathR on 15 Nov 2025
(the log records the interaction and annotation tables loaded from the
OmniPath web service cache on that date):

| Layer | Network | Call in `activities_splitting_data_into_two_groups.R` | Method |
|---|---|---|---|
| Pathway | PROGENy, top 500 genes per pathway (14 pathways) | `decoupleR::get_progeny(organism = "human", top = 500)` | `run_mlm` |
| Kinase | OmniPath kinase-substrate interactions, datasets `kinaseextra` + `omnipath`, kinases with >= 50 targets, expressed above the median | `OmnipathR::omnipath_interactions(datasets = c("kinaseextra", "omnipath"))` | `run_ulm`; 14 most variable kinases kept per risk group |
| Transcription factor | DoRothEA regulons, confidence levels A and B | `decoupleR::get_dorothea(organism = "human", levels = c("A", "B"))` | `run_mlm`; 14 most variable TFs kept per risk group |

OmniPath is a continuously updated web resource. The exact interaction tables
are therefore the ones served on 15 Nov 2025; the log records the cache record
identifiers and row counts (101,454 and 278,677 interactions; 699,831 annotation
records).

## Graph construction, model fits and audit (Nov 2025 to Aug 2026)

| Component | Version |
|---|---|
| R | 4.4.x (Windows); the graph-on-vector fits and the audit ran under R 4.4 with the packages below |
| MASS | 7.3-64 / 7.3-65 |
| matrixcalc | 1.0-6 |
| mvtnorm | 1.3-x |
| invgamma | 1.1 / 1.2 |
| MCMCpack | 1.7-1 |
| igraph | 2.1.4 to 2.3.0 |
| softImpute | version not recorded in a log or saved object |
| coda | 0.19-4.1 |
| ggplot2 | 3.5.1 / 4.0.3 |

The 84 saved fit objects record `burn = 500` and `nrun = 1000` but do not store
`sessionInfo()`, so the patch-level versions of these packages at fit time
(Nov 2025) cannot be confirmed; the ranges above span the R 4.4 libraries on the
analysis machine. None of the model code depends on package behaviour that
changed across these ranges (base linear algebra, `mvtnorm::rmvnorm`,
`MCMCpack::rinvgamma`, `MASS::ginv`).

## Prior hyperparameters of the deposited model code

`R/lsGPR_VS.R` (and `analysis/03_graph_on_vector_fits/lsGPR_VS_as_run.R`) use
`a.t = 15, b.t = 1` (tau), `a.p1 = 2, b.p1 = 50` (psi1), `a.p = 1, b.p = 10`
(prior inclusion probability), proposal step 0.01, and initial values
`c(1, 1, 1, 1)`. The saved fit objects do not store the prior list, so the
manuscript reports these deposited-code values. The two files differ only by a
`control` list in the as-run copy (per-parameter proposal steps, acceptance
tracking, progress switch) that was added for the NFkB proposal-tuned refit.
