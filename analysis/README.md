# Analysis pipeline

Scripts are deposited as run. File paths inside them are absolute paths on the
analysis machine (a Dropbox folder named `H&E Lung TCGA Graph Dataset`); the
table below gives the inputs and outputs each script expects so the paths can be
remapped. Exploratory sections that did not feed the paper (LUSC comparisons,
alternative plots) are left in place rather than edited out.

## 01_data_preparation

`activities_splitting_data_into_two_groups.R`

1. Reads TCGA-LUAD raw counts (418 cases), slide metadata and the
   whole-slide-image to RNA-seq mapping; keeps cases with both data types.
2. Splits cases at the median survival rank into high-risk and low-risk groups
   (402 cases: 200 high-risk, 202 low-risk) and writes `high_risk_luad.rds`,
   `low_risk_luad.rds`.
3. Log-normalises counts per group (`log_cts_high_risk_luad.rds`,
   `log_cts_low_risk_luad.rds`).
4. Infers activities with decoupleR (see `ENVIRONMENT.md` for networks and
   versions): PROGENy pathways with `run_mlm`
   (`LUAD_{high,low}_risk_activity_scores.rds`); OmniPath kinase-substrate
   network with `run_ulm`, 14 most variable kinases kept
   (`top_kinase_df_{high,low}risk_luad.rds`); DoRothEA A/B regulons with
   `run_mlm`, 14 most variable TFs kept (`top_tf_df_{high,low}risk_luad.rds`).

`omnipathr-20251115-1340_activity_inference_session.log` is the OmnipathR log
of the session that produced these files.

## 02_graph_construction

`run_on_files_obtain_U_and_a.R`

For each slide's 20 x 20 HPC adjacency matrix (G-cross statistic at radius
R = 5 between HPC point patterns; one CSV per slide, not distributed here):
rows and columns that are entirely zero are set to missing and completed with
Soft-Impute (`softImpute(rank.max = 19, lambda = 0.1, type = "svd")`), the
diagonal is set to 1 and negatives to 0, the matrix is symmetrised by the
geometric mean of positive bilateral entries with positive unilateral entries
kept through the pairwise maximum, and `EM_toU(adj, d = 3, u1 = 0.5, sig.a = 2,
sig.u = 2, tol = 1e-10, maxiter = 3)` returns the density intercept `a` and the
latent positions `U`. Outputs per risk group: `luad_{HR,LR}_PATIDs.rds`,
`luad_{HR,LR}_mya.rds`, `luad_{HR,LR}_myU.rds` (194 high-risk and 190 low-risk
aligned cases).

`constrained_adjacency_covariance_function.R` is sourced by that script
(graph-kernel utilities used for visualisation).

## 03_graph_on_vector_fits

`pathway_activity_regression.R`, `kinase_activity_regression.R`,
`transcription_factor_activity_regression.R`

For each risk group and each of the 14 activities in a layer: the first 10
aligned cases are held out, and
`lsGPR_VS(a, a.te, U, U.te, z = NULL, z.te = NULL, Y, burn = 500, nrun = 1000)`
is fitted once (single chain). 2 groups x 3 layers x 14 activities = 84 saved
fits (`luad/{HR,LR}_luad_[kinase_|TF_]lsGPRVS_result1000_<activity>_activity.rds`).
The remainder of each script produces the original PIP heatmaps; the published
Figure 2 uses the post-burn-in recomputation in `04_audit/07_make_figure2_panels.R`.

`lsGPR_VS_as_run.R` is the copy of the model function that these scripts
`source()`; it matches `R/lsGPR_VS.R` except for the `control` list described in
`ENVIRONMENT.md`.

## 04_audit

Run from the project root (the scripts locate `luad/`, `TCGA-LUAD-data/` and
`adj_matricesV3/` relative to the working directory).

| Script | Reads | Writes |
|---|---|---|
| `01_audit_saved_tcga_models.R` | the 84 saved fits | `results/audit/post_burnin_pips.csv`, `pip_threshold_sensitivity.csv`, `posterior_expected_fdr.csv`, `heldout_prediction_metrics_by_activity.csv`, `heldout_prediction_summary.csv` (Tables 1 and 2) |
| `02_saved_fit_mcmc_diagnostics.R` | the 84 saved fits | `saved_fit_mcmc_diagnostics_by_activity.csv`, `saved_fit_mcmc_diagnostics_summary.csv`, `saved_fit_pip_mcmc_diagnostics.csv`, `saved_fit_pip_mcmc_diagnostics_summary.csv` (Table 3) |
| `03_tuned_refit_pilot.R` | saved latent representations and the high-risk NFkB response | `results/nfkb_tuned_refit/tuned_refit_grid_summary.csv`, `tuned_refit_final_summary.csv` (proposal SD grid, then the final 1,500-iteration refit with steps 1.5/0.5) |
| `04_make_traceplot.R` | original and tuned NFkB fits | supplementary trace plot |
| `06_generate_postburnin_heatmaps.R` | `post_burnin_pips.csv` | supplementary PIP heatmaps |
| `07_make_figure2_panels.R` | `post_burnin_pips.csv` | Figure 2 |
| `08_make_figure3_panels.R` | original and tuned NFkB fits | Figure 3 |

Figure 1 is a schematic assembled in Adobe Illustrator from the survival curve,
example slides and the Figure 2 heatmaps; it has no generating script.

Script `05` of the original numbering was an intermediate that was superseded
by `06`; the numbering is kept so that it matches the analysis record.
