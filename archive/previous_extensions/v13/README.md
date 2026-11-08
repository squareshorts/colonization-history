# Bacterial colonization carrier discrimination: reproducibility extension

The scientific manuscript, supplementary material and bibliography are provided in the separate Overleaf project. This research-code package contains the additional independent analysis of two-arm recovery, common cell-area calibration and power under probabilistic missing-outcome assumptions.

## Reproduce

Python 3.10+ and `numpy`, `scipy`, `matplotlib` are required.

```
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python analysis/fine_recovery.py
```

The final command reads the saved fine-grid CSV if present, then independently generates the probability-calibrated assay outcomes and new figures. To regenerate the entire fine loss grid, remove `results/fine_loss_detection.csv` and run the command again. All experiments use fixed seeds. Source simulations are exploratory and their parameters were not fitted to the original bacterial microscopy.

Reproduce the single-cell reference calibration against a checkout of the existing public repository, which already includes the derived microscopy observables:

```
python analysis/shared_calibration.py --repo /path/to/colonization-history --out results
```

The corresponding `counts.csv` files were extracted from the Figshare binary image deposit by `empirical/extract.py` in the public repository. The previously produced `meta.json` files contain the single-cell-area references. The supplied new script reads these released derived records; it does not resegment the original microscopy images.

## Inference scope

The exact finite-sample missing-outcome bounds require only binary outcomes and original assignment counts. The probability-assisted bound uses an *additional* assumption of conditional Bernoulli independence and externally verified upper probabilities for the missing sham outcomes. The model with valid upper bound `u=0.10` is deliberately constructed from hypothetical sham types with `p_H=0.10`, `p_L=0.03`. A deliberately invalid `u=0.03` is included as a failure control. Treat the results as design sensitivity, not empirical assay performance.

The main repository's old Zenodo record https://doi.org/10.5281/zenodo.23990601 predates these analyses. No new DOI was minted.

Source-format checks in `tests/test_v13.py` are automatically skipped if the separate Overleaf project is absent. Scientific and numerical tests run without the manuscript. To exercise all checks, place the extracted Overleaf project in an `overleaf/` subdirectory of this reproducibility folder.
