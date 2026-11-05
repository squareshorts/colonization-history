# Colonization history

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23990601.svg)](https://doi.org/10.5281/zenodo.23990601)

Code, simulation outputs and derived observables for attributing bacterial surface
colonization history to cells or surfaces. One stochastic attachment–division–detachment
simulator implements six hypotheses, calibrated to a common fresh-cell baseline:

| Hypothesis | Transfer effect arises from |
|---|---|
| Heterogeneous passive | nothing (null with heterogeneous patch and cell rates) |
| Residual-cell carryover | residents left on a visited surface below detection |
| Adhesive selection | exposed cells enriched for heritably adhesive variants |
| Cell history | a state acquired on the surface, inherited, decaying in suspension |
| Surface history | a state deposited by contacts and residents, decaying on the surface |
| Cell + surface history | both history states |

The repository contains three analyses.

1. **Calibrated factorial comparison** (`calibrate.py`, `run_experiments.py`). History
   hypotheses receive baseline rates fitted so that their fresh-cell / fresh-surface
   trajectory matches the passive null; confounder strengths reproduce the zero-gap
   history effects.
2. **Mechanism recovery** (`recovery.py`). Random-forest model choice on 3000 synthetic
   experiments per hypothesis, baseline-matched, across designs (replicates, gaps, arms)
   and observation tiers (counts against lineage tracking).
3. **Re-analysis of Lee et al. (2018)** (`empirical/`). Observable extraction from the
   deposited image series (doi:10.6084/m9.figshare.5969335), dose and lag validation,
   residence-time fits that compare exposure as an adhesiveness offset with exposure as a
   head start on surface-acquired history, time-matched and occupancy-adjusted fits, a
   pre-specified cAMP test, a posterior predictive check and prospective power.

Model quantities are dimensionless. No simulator parameter is fitted to the image data;
the empirical analysis uses statistical families with the simulator's residence law.

## Requirements

Python 3.10 or later with the packages in `requirements.txt`:

    pip install -r requirements.txt

## Reproduce

    python -m unittest test_model -v          # 18 tests, about 2 s
    python calibrate.py                        # about 15 min
    python run_experiments.py                  # about 2 min
    python recovery.py --n 3000 --workers 2    # about 25 min on 2 cores
    python make_figures.py && python make_tables.py

Empirical analyses from the observables shipped in `empirical/observables`:

    python empirical/analyze.py 60             # first-3-h windows, about 30 min
    python empirical/occupancy.py              # lags, dose and their correlations
    empirical/run_v2_parts.sh 200              # full-record, time-matched and cAMP analyses, about 30 min on 2 cores
    python empirical/design_power.py           # prospective power, about 10 min

The extended analyses follow `empirical/ANALYSIS_PLAN_v2.md`, written before they were
run; the same file logs results, deviations and the analyses added afterwards.

To re-extract observables from the raw deposit (2.5 GB), download it from figshare,
unzip the outer archive into a directory and run

    empirical/run_all.sh <directory with the 16 dataset zips>   # windowed tracking, about 3 min per dataset
    empirical/run_v2.sh  <directory with the 16 dataset zips>   # continuous tracking with fluorescence readings

`Colonization_History.ipynb` (executed) runs the tests, a reduced factorial, a reduced
recovery study and the empirical fits from the shipped observables in one session.
`make_tables.py` writes LaTeX tables to `manuscript/generated/`, which is not tracked.

## Layout

- `model.py` simulator with per-replicate parameters and exact cell accounting
- `calibrate.py` baseline and confounder calibration
- `run_experiments.py` factorial comparison, step-size and effect-strength checks
- `recovery.py` synthetic experiments, baseline matching, classifier evaluation
- `empirical/extract.py`, `extract_v2.py` windowed and continuous tracking, fluorescence registration
- `empirical/analyze.py`, `analyze_v2.py`, `run_v2_final.py` residence-time families, bootstrap,
  time-matched fits, cAMP test
- `empirical/occupancy.py`, `design_power.py` lag validation and prospective power
- `make_figures.py`, `make_tables.py`, `make_notebook.py` figures, tables and notebook
- `results/` tables, logs, figures and fits; `results/recovery_pool.npz` holds the
  synthetic experiments used by the classifier
- `SHA256SUMS.txt` checksums of every tracked file

## Citation and licence

Cite the archived release, doi:10.5281/zenodo.23990601; `CITATION.cff` gives the metadata.
Code is MIT-licensed (`LICENSE`). Results and derived observables are CC BY 4.0
(`LICENSE-DATA.md`); the observables derive from Lee *et al.*, figshare
doi:10.6084/m9.figshare.5969335 and PNAS 115:4471 (2018), doi:10.1073/pnas.1720071115,
which should be cited alongside this repository.

## Version history

- **1.1.0** Continuous tracking of every arrival (22,888 attached cells), time-matched and
  occupancy-adjusted fits, pre-specified cAMP test (uninformative: the constitutive channel
  is at the noise floor), posterior predictive check and prospective power.
- **1.0.0** Heritable cell heterogeneity in every hypothesis, carryover and selection
  confounders, calibrated baselines, recovery study and empirical re-analysis.
- **0.1.0** Four mechanisms at shared rates (superseded).
