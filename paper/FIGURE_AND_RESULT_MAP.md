# Figure and numerical-result provenance

The current source for reproducible calculations is in `paper/` on the public repository, and the earlier exploratory scripts remain accessible under `archive/previous_extensions/`. The current article and supplement are distributed separately in an Overleaf ZIP, outside GitHub. The microscopy input tables are the previously segmented derivatives in `empirical/observables/`, whose source is Lee et al.'s Figshare image deposit. No segmentation of raw images was repeated for the corrected occupancy estimates.

| Location in article | Quantitative source | Generation or check |
|---|---|---|
| Main Fig. 1, directly observed sustained fivefold establishment in paired flow cells and FC3 | `paper/results/observed_establishment.csv` | `paper/analysis/observed_establishment.py` followed by `paper/analysis/make_observed_figure.py` |
| Supplement section 2, archival estimator check and secondary exponential growth-lag fits | `paper/results/shared_calibration.csv`, `paper/results/edge_corrected_summary.json` | `paper/analysis/shared_calibration.py` followed by `paper/analysis/make_lag_figure.py` |
| Main Fig. 2, all-time lag-carrier equivalence | Exact selection/Hill identity in Methods and Supplement section 3 | `archive/previous_extensions/v11/analysis/make_figures.py`; mathematical comparator in `archive/previous_extensions/v9/analysis/discriminating_test_v9.py` |
| Main Fig. 3, bounded attrition-power frontier | `paper/results/fine_loss_detection.csv` | `paper/analysis/fine_recovery.py` and `paper/analysis/model_core.py` |
| Main Fig. 4, repeated-contact information | `paper/results/repeated_contact_information.csv` | `paper/analysis/fine_recovery.py` and `paper/analysis/model_core.py` |
| Supplement Table 1, directly observed fold crossings | `paper/results/observed_establishment.csv` | `paper/analysis/observed_establishment.py` |
| Supplement Table 2, retrospective fold threshold sensitivity | `paper/results/observed_establishment_sensitivity.csv` | `paper/analysis/observed_establishment.py` |
| Supplement Table 3, occupancy lags and technical ranges | `paper/results/shared_calibration.csv`, `paper/results/technical_fit_grid.csv` | `paper/analysis/shared_calibration.py`; also checks FC3, common area, and end-window weights |
| Supplement section 1, log-convexity and decay constraint | Published fitted lag point estimates from Lee et al. (2018) | Analytic inequalities in Supplement section 1; original numerical scenario in `archive/previous_extensions/v8/analysis/lag_robustness_v8.py` |
| Supplement Table 2 and Fig. S1, two-arm recovery controls | `paper/results/two_arm_loss_grid.csv` | `paper/analysis/model_core.py`; the original figure generator is in `archive/previous_extensions/v12/analysis/critical_design_v12.py` |
| Supplement Fig. S2, calibrated missingness | `paper/results/assisted_calibration.csv` | `paper/analysis/fine_recovery.py`, with null false-positive stress from `paper/analysis/calibration_stress.py` |
| Supplement Fig. S3, repeated-contact oracle power | `archive/previous_extensions/v12/results/v12_repeated_contact_oracle_mc.csv` | `archive/previous_extensions/v12/analysis/critical_design_v12.py` |
| Supplement Fig. S4, localization selectivity | Assumed five-carrier matrix and calibrated assay coefficients in supplement | `archive/previous_extensions/v11/analysis/make_figures.py`, with stochastic coverage from `archive/previous_extensions/v10/analysis/generality_attrition_v10.py` |
| Fine-grid Monte Carlo percentages | `paper/results/fine_loss_detection.csv` | `paper/analysis/fine_recovery.py` |
| 32,000-study seed replication | `paper/results/cloglog_seed_replication.json` | `paper/analysis/repeat_cloglog_check.py` |
| Calibration-failure false-attribution rates | `paper/results/null_asymmetric_calibration_stress.csv` | `paper/analysis/calibration_stress.py` |

The old DOI `10.5281/zenodo.23990601` archives an earlier release and does not archive the current article's code. The exact code revision is the Git commit linked in the article's data statement. A current version-specific archival DOI remains pending, as requested by the author.
