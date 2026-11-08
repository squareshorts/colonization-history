# Bacterial colonization-history inference: independent computational analyses

The package reproduces an exploratory, mathematical and simulation-based investigation of mechanistic identifiability in bacterial colonization history.

## Completed computations

1. `analysis/archival_lag_audit.py` summarizes per-independent-field lag estimates from derived occupancy records supplied in the public repository. The 14-hour and 37-hour holding conditions each have one deposited run; their between-replicate variance cannot be estimated. Those lags use a different estimator from the original published exponential fit.
2. `analysis/recompute_occupancy_lags.py` recomputes window/smoothing sensitivity from the public `empirical/observables/` folder of the source repository. The source-Git blob SHAs are recorded in `results/technical_fit_summary.csv`. The script needs a checkout of `squareshorts/colonization-history` containing the original derived CSVs. Raw microscope images were not downloaded or reprocessed.
3. `analysis/attrition_power_and_observables.py` simulates 96 attrition/effect/kinetics experimental designs (1,000 independent study replicates each, 32 clonal blocks with 48 originally assigned founders per arm) and generates the paired-contact identifiability benchmark. The programs require NumPy, SciPy, pandas and Matplotlib.
4. `analysis/make_figures.py` and `analysis/make_supplement_figures.py` regenerate the main and supplementary graphics.
5. `tests/test_design.py` verifies exact model equivalence, positive-mixture inequality, sharp missing-data limits, expected power/bias behavior, archive unit count, and manuscript preflight assertions.

From this directory run:

    python -m pip install -r requirements.txt
    python analysis/archival_lag_audit.py
    python analysis/attrition_power_and_observables.py --studies 1000 --seed 20261107
    python analysis/make_figures.py
    python analysis/make_supplement_figures.py
    python -m pytest -q tests

The synthetic model parameters were selected to illustrate inferential constraints. The generated power numbers are not calibrated to an actual bacterial assay. Simulated repeatedly measured attachment outcomes assume stable cell traits and conditionally independent contacts. The mathematical result is input-output non-identifiability for specified observation classes, not a claim of a novel theorem or a universal biological mechanism.

## Data provenance

Archival occupancy data are derived from Lee et al. (2018), DOI 10.1073/pnas.1720071115, original dataset DOI 10.6084/m9.figshare.5969335. The derived records were retrieved from https://github.com/squareshorts/colonization-history. The snapshots and blob provenance are recorded in `results/archive_provenance.md`. The 0-, 14- and 37-hour published fitted lags are taken from the original article.

The code and manuscript were prepared with substantial generative AI assistance. Scientific interpretations and analytical choices require independent author review before submission, particularly if the target publisher restricts generative AI in substantive research tasks.
