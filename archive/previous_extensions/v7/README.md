# Colonization-history analysis patch, v7

This is an **additive** analysis package for `squareshorts/colonization-history` (baseline source commit `23cf4552cf2b1f2666e362b03fe20dd3fe015823`). It is not a replacement for the original repository, the v5 nonlinear simulator patch, or the original image data. The scripts and outputs support a narrower methods paper on the limits of causal carrier attribution.

Files provided:

- `empirical/field_level_v7.py` — reads the original derived `tracks_full.csv` files; matched field-level KM and attachment-definition sensitivity.
- `analysis/import_archival_field_v7.py` — independently verified aggregate archival observations with Git blob hashes; regenerates the included audit table/figure without shipping raw derived cell records.
- `analysis/calibration_bounds_v7.py` — 1,000-replicate calibration-aware partial-identification benchmark with miscalibration and nonlinearity negative controls.
- `analysis/selectivity_sweep_v7.py` — 26-cell cross-reactivity/replication sweep (500 synthetic experiments per cell).
- `tests/test_v7.py` — nine regression and numerical tests.
- `results/` — executed results and provenance; `manuscript/V7_METHODS_RESULTS.tex` is a freestanding, explicitly scoped manuscript section draft.
- `REPORT_v7.md` — mathematical assumptions, exact results, interpretation and remaining publication blockers.

Run from the package directory:

```bash
python -m pip install numpy scipy pandas matplotlib pytest
python -m pytest tests -q
python analysis/calibration_bounds_v7.py --replicates 1000 --n 16 64 160
python analysis/import_archival_field_v7.py
python analysis/selectivity_sweep_v7.py --replicates 500 --n 32 128
```

To rerun the complete field-level pipeline from individual tracks, copy the additive scripts into the **original** repository root and execute:

```bash
python empirical/field_level_v7.py --repo-root . --outdir results/field_level_v7
```

No manuscript, GitHub remote, Zenodo record, or original model files were overwritten. All synthetic parameters and scenario settings are declared in code. Experimental data were obtained from the published derived observables of Lee et al. (2018) in the source repository. The images and original full microscopy dataset were not reprocessed in this environment.
