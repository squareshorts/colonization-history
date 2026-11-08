# Bacterial colonization history -- calibrated-discrimination analysis

This reproducibility package accompanies the theoretical bacterial-colonization-history manuscript. It provides the independently seeded 32,000-study check of a repeated 34% result and a no-effect stress test of invalid calibrated missing-outcome probabilities. The earlier full loss-grid, known-parameter repeated-contact, and shared-area-calibration outputs are retained as comparative inputs. The source of all microscopy-derived counts remains the public `empirical/observables` directory of https://github.com/squareshorts/colonization-history.

Install Python 3.10+ with `numpy`, `scipy`, and `matplotlib`:

```
python -m pip install -r requirements.txt
python -m unittest discover -s tests -p 'test_v*.py' -v
python analysis/calibration_stress.py
```

The stress test is seeded and rewrites `results/null_asymmetric_calibration_stress.csv`. For the independently seeded larger check, the command below regenerates the 32,000 synthetic studies in four independent batches:

```
PYTHONPATH=analysis python analysis/repeat_cloglog_check.py
```

The other model and figure scripts remain in `analysis/`. In a checkout of the main public repository the occupancy calibration may be recomputed using:

```
python analysis/shared_calibration.py --repo /path/to/colonization-history --out results
```

The main manuscript, supplement and bibliography are distributed in a separate editable Overleaf ZIP. They are not included in the public GitHub project. Source-format tests skip automatically if that separate `overleaf/` folder has not been supplied.

The previous Zenodo software record https://doi.org/10.5281/zenodo.23990601 predates the added analyses. The current code has not been independently archived with a new DOI.
