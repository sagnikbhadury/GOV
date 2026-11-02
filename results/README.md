# Processed outputs

All files are CSV. `risk_group` is `HR` (high-risk) or `LR` (low-risk);
`molecular_layer` is `Pathway`, `Kinase` or `Transcription factor`; `hpc` is one
of the 20 survival-associated HPCs (`hpc_0` ... `hpc_44`).

## audit/

| File | Rows | Content |
|---|---|---|
| `post_burnin_pips.csv` | 1,680 | Posterior inclusion probability of each HPC in each of the 84 fits, from iterations 501 to 1,000. At PIP > 0.5 this table gives the counts reported in the Abstract: pathway 97 vs 35, kinase 18 vs 2, transcription factor 49 vs 43 (high-risk vs low-risk). |
| `pip_threshold_sensitivity.csv` | 24 | Number of selected HPC-activity pairs at PIP thresholds 0.3, 0.5, 0.7 and 0.9 per group and layer (Table 1). |
| `published_pip_threshold_sensitivity.csv` | 24 | The same quantity computed from the full-chain PIPs of the first submission, kept for comparison. |
| `posterior_expected_fdr.csv` | 1,680 | Posterior expected false-discovery proportion of the selected set when HPCs are ranked by PIP, and the selection at a 5% posterior expected FDP. |
| `heldout_prediction_metrics_by_activity.csv` | 84 | RMSE, MAE, test R^2, correlation, 95% interval coverage and the improvement over the training-mean baseline for the 10 held-out cases of each fit. |
| `heldout_prediction_summary.csv` | 6 | Medians of the above per group and layer (Table 2). |
| `saved_fit_mcmc_diagnostics_by_activity.csv` | 420 | Lag-1 autocorrelation and effective sample size of the log-likelihood and the inverse length-scale chains in each fit. |
| `saved_fit_mcmc_diagnostics_summary.csv` | 30 | Medians and minima per group, layer and parameter (Table 3). |
| `saved_fit_pip_mcmc_diagnostics.csv` | 1,680 | ESS and Monte Carlo standard error of each inclusion-indicator chain. |
| `saved_fit_pip_mcmc_diagnostics_summary.csv` | 6 | Summary of the indicator-chain diagnostics per group and layer. |

These are single-chain diagnostics from 500 post-burn-in draws. They describe
mixing; they cannot establish convergence.

## nfkb_tuned_refit/

| File | Content |
|---|---|
| `tuned_refit_grid_summary.csv` | Pilot grid over a common proposal step (0.03, 0.1, 0.3, 1) for the high-risk NFkB pathway model: acceptance rates, ESS, correlation with the original PIPs, selected-node counts. |
| `tuned_refit_final_summary.csv` | The final refit with steps 1.5 (density scale) and 0.5 (network scale), 1,500 iterations and 500 burn-in: acceptance 0.367 / 0.459, ESS 95.1 / 18.9, PIP correlation 0.477 with the original chain, 7 nodes selected versus 19. |

The refit is a sensitivity analysis. It does not replace the primary PIPs in
`audit/post_burnin_pips.csv`.
