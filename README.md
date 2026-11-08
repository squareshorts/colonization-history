# Bacterial colonization history: mechanisms and discriminating measurements

This repository provides the computational source, processed microscopy observations, numerical results and tests associated with a theoretical analysis of history-dependent bacterial surface colonization. The central question is which observations can distinguish a transient state acquired by a cell from selection of pre-existing adhesive phenotypes, cell-associated extracellular material, or a history-dependent environment.

Distinct intracellular-state and differential-growth-selection models can produce exactly the same fitted colonization-lag trajectory. Randomized, founder-indexed contact assays and additional carrier-selective measurements can separate restricted alternatives under explicit assumptions. Missing attachment outcomes, especially among sham cells, can severely reduce power and lead to false attribution if calibration assumptions are violated.

## Current analyses

- [`paper/`](paper/) contains the current numerical source, results and tests. [`paper/FIGURE_AND_RESULT_MAP.md`](paper/FIGURE_AND_RESULT_MAP.md) maps every article figure and main numerical result to its generating program or the traceable archived source.
- [`empirical/observables/`](empirical/observables/) contains previously derived occupancy and tracking measurements from binary microscopy deposited by Lee et al. (2018). [`empirical/extract.py`](empirical/extract.py) documents the original segmentation process; [`empirical/occupancy.py`](empirical/occupancy.py) provides the corrected boundary-normalized smoothing estimator.
- [`archive/previous_extensions/`](archive/previous_extensions/) preserves previous exploratory analyses, agent ablations and alternative numerical experiments. These materials remain traceable for research history and are not the primary evidence for the current article.
- The original prototype simulator and its exploratory scripts remain in the repository history and at the root for backward compatibility.

## Reproduce the reported results

Python 3.10 or later is required. From the repository root:

```sh
python -m pip install -r paper/requirements.txt
python -m unittest discover -s paper/tests -p 'test_*.py' -v
python paper/analysis/observed_establishment.py --repo . --out paper/results/observed_establishment.csv --sensitivity paper/results/observed_establishment_sensitivity.csv
python paper/analysis/endpoint_observation_comparison.py --repo . --records paper/results/endpoint_observation_comparison.csv --summary paper/results/endpoint_observation_summary.json --sensitivity paper/results/absolute_threshold_sensitivity.csv
python paper/analysis/make_endpoint_comparison_figure.py --csv paper/results/endpoint_observation_comparison.csv --out paper/figures/endpoint_comparison.pdf
python paper/analysis/shared_calibration.py --repo . --out paper/results
python paper/analysis/make_lag_figure.py --summary paper/results/shared_calibration.csv --output paper/figures/matched_flowcell_lags.pdf
```

The tests check the analytic bounds, saved Monte Carlo results, smoothing edge behavior, the matched comparisons, FC3 inclusion, and censoring-aware common occupancy thresholds. Relative and common-area endpoints are retrospective measures that also reflect initial cell delivery. Computational runs with alternative seeds and larger replication counts are described in [`paper/README.md`](paper/README.md). The article and electronic supplement are maintained outside the GitHub repository.

## Data provenance

The primary binary microscopy images are available at [Figshare DOI 10.6084/m9.figshare.5969335](https://doi.org/10.6084/m9.figshare.5969335). The occupancy recordings in this repository were extracted from the deposited images using the repository's segmentation software. They were not supplied as counts by the original depositors. Derived occupancy lags use an operational estimator distinct from the growth-lag procedure reported by Lee et al. Biological replication for the individual 14-hour and 37-hour holding records is unavailable in the deposited set.

## Code version and citation

The current source version is **1.2.1**. The exact, immutable Git commit for this code is referenced in the corresponding article's Data accessibility statement. `CITATION.cff` describes this software version without attaching the DOI of a different release.

The existing [Zenodo DOI 10.5281/zenodo.23990601](https://doi.org/10.5281/zenodo.23990601) is for the **earlier v1.1.0** release. It does not archive the current results. A new archival DOI has not been minted.

Software is MIT-licensed (`LICENSE`); derived observations and results are governed by `LICENSE-DATA.md`. Cite the source article and the Figshare microscopy deposit when using the image-derived data.
