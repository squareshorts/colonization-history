# Computational analyses for bacterial colonization-history discrimination

This directory contains the current article's numerical analysis and tests. Simulated attachment and recovery models are mathematical demonstrations and have not been calibrated to a bacterial population. The empirical observations are independent image-derived `counts.csv` and `meta.json` files created earlier from the Lee et al. (2018) Figshare binary microscopy stacks. They are retained in the public repository's `empirical/observables/` directories. The current edge-correction analysis does not repeat raw-image segmentation.

## Reproduction

In a checkout of the public repository `squareshorts/colonization-history`, use Python 3.10+.

```sh
python -m pip install -r paper/requirements.txt
python -m unittest discover -s paper/tests -p 'test_*.py' -v
python paper/analysis/shared_calibration.py --repo . --out paper/results
python paper/analysis/make_lag_figure.py --summary paper/results/shared_calibration.csv --output paper/figures/matched_flowcell_lags.pdf
```

The seed-replication and calibration-stress jobs are optional longer runs and overwrite their respective derived output files.

```sh
PYTHONPATH=paper/analysis python paper/analysis/repeat_cloglog_check.py
python paper/analysis/calibration_stress.py
```

`FIGURE_AND_RESULT_MAP.md` maps every main and supplementary figure and the central numerical claims to the exact scripts and tables. The earlier exploratory analyses remain archived under `archive/previous_extensions/`, and the main article LaTeX sources are supplied privately as an Overleaf ZIP. The README and CITATION metadata identify the current source revision separately from the earlier archived DOI.

## Directly observed surface occupancy and observation-rule sensitivity

Archived occupied-area records have substantially different early cell-area references in FC1 and FC2. The current article compares the original retrospective fivefold increase against common absolute thresholds of 20, 30 and 50 cell-area equivalents, applying a shared 262-pixel calibration. At 20 equivalents, FC2 series 1 is already above threshold in the initial frame, so its new crossing time is left-censored. At 50 equivalents, FC1 series 2 remains below threshold to the end of the recording, so the matched difference is a bound. The sign of the transfer-associated difference is negative under all tested thresholds, smoothing widths and persistence requirements. The inferred normalized holding-time advantage depends on the observation rule, with 37-hour values ranging from approximately 0.30 to 0.70.

To reproduce the archival contrasts:

```sh
python paper/analysis/observed_establishment.py --repo . --out paper/results/observed_establishment.csv --sensitivity paper/results/observed_establishment_sensitivity.csv
python paper/analysis/endpoint_observation_comparison.py --repo . --records paper/results/endpoint_observation_comparison.csv --summary paper/results/endpoint_observation_summary.json --sensitivity paper/results/absolute_threshold_sensitivity.csv
python paper/analysis/make_endpoint_comparison_figure.py --csv paper/results/endpoint_observation_comparison.csv --out paper/figures/endpoint_comparison.pdf
```

The threshold results summarize occupancy progression and cannot attribute its biological carrier. No new experiment or longer recording of FC1 series 2 was performed. The original exponential-lag calculation remains secondary and extrapolation-dependent.
