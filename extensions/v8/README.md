# Colonization history: v8 adversarial analysis and manuscript revision

Authoring status: 2026-11-07. Researcher's existing original GitHub repository is squareshorts/colonization-history. This additive package does not change the original v1.1.0 simulator or data.

## Key files

- `overleaf/main.tex`, `overleaf/supplement.tex`, `overleaf/references.bib`, and `overleaf/figures/` form a complete Overleaf project. Compile with pdfLaTeX -> Biber -> pdfLaTeX -> pdfLaTeX separately for main and supplement.
- `HOSTILE_REVIEW_V8.md` is the point-by-point reviewer audit, with unresolved empirical limitations and CAPES journal status.
- `analysis/` contains source for lag falsification, bounded perturbations, a composition-only selection counterexample, 160-parameter agent ablation, endpoint sensitivity, power sensitivity and tests.
- `results/v8/` contains generated summary JSON/CSV files. Raw images from Lee et al. are not redistributed; original derived tracks remain in the existing GitHub repository.
- `agents/`, `v7/` preserve the predecessor source and reports.

## Reproduction

Python 3.10+ with numpy, scipy, matplotlib, pandas if regenerating plots. From this package root:

```
python analysis/hostile_checks_v8.py --draws 160 --reps 5 --seed 2026110709
python analysis/lag_robustness_v8.py
python analysis/selection_counterexample_v8.py
python analysis/endpoint_robustness_v8.py --n 16
python analysis/power_sensitivity_v8.py
python analysis/figures_v8.py
python -m unittest analysis/test_v8.py -v
```

The original v7 checks pass with `cd v7 && python -m pytest tests -q`; the prior agent checks with `cd agents && python -m unittest discover -s tests -v`. The two article PDFs are generated from Overleaf source, and included in the package. The 37-hour result was not prospectively held out, since the model was developed knowing the published outcome. Fitted selection model parameters are illustratively tuned to the reported points. Publication venue CAPES eligibility is conditional on a final journal-level check.
