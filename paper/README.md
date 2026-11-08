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

## Directly observed surface-occupancy expansion

The current manuscript's main empirical endpoint is time to a sustained fivefold increase in recorded occupied area, measured within the archived flow-cell time series. The 26.75-hour FC1 series-2 recording reaches the threshold and contributes to the three paired comparisons without extrapolating its narrowly supported exponential lag. Recompute with `python paper/analysis/observed_establishment.py --repo . --out paper/results/observed_establishment.csv --sensitivity paper/results/observed_establishment_sensitivity.csv`. The fitted exponential lags remain secondary sensitivity results. The observed threshold combines initial surface delivery with later attachment and population expansion; it is a retrospective operational endpoint, not a biological memory marker.
