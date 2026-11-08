# Colonization carrier discrimination: corrected computational analyses

## Scope

Two complementary analyses support a conservative assay-design argument.

1. **Missing-outcome attribution under two-arm loss.** `analysis/critical_design_v12.py` simulates 280 designs crossing two attachment/residence laws, seven recovery scenarios, five nominal loss levels, and four induced effect sizes, each with 1,200 independent synthetic studies of 32 clone blocks. Full originally assigned outcomes are hidden from the test but generated to validate the sharp identification bounds. A missing-at-random assumption is not used to derive bounds. The alternative recovered-only analysis is provided to reveal the bias-power trade-off.
2. **Repeated-contact information at simulator-scale probabilities.** The repeated-measurement analysis uses high- and low-phenotype attachment probabilities of 0.10 and 0.03. A known-parameter likelihood-ratio benchmark estimates 80% discrimination sample sizes, with a Monte Carlo check and sensitivity to other probabilities.
3. **Independent archival estimator.** `analysis/archival_within_estimator_v12.py` recomputes ratios from the eight repository-derived occupancy-lag estimates. It also verifies that the technical ranges for matched FC1/FC2 recordings remain separated and documents why two minima round to 4.97 h. Full source segmentation was performed by the public `empirical/extract.py` against the archived binary microscopy stacks, with derived `counts.csv` and `meta.json` in the source repo. This extension reuses those count records and has not independently regenerated the image segmentation.

The published Lee et al. 2018 exponential-fit lags and the derived area-occupancy-crossing lags are different estimators. The former are used only for an algebraic model illustration. The latter support a descriptive independent-flow-cell comparison. The 14- and 37-hour holds each have one available recorded run and no between-replicate confidence interval can be estimated from them.

## Reproduction

Requires Python >=3.10 with `numpy`, `scipy`, `matplotlib`, and `pytest`.

```bash
python analysis/critical_design_v12.py
python analysis/archival_within_estimator_v12.py
python -m pytest tests/test_critical_design_v12.py -q
```

The `results/` directory includes the same numerical summaries saved from the tested random seeds. All figures are programmatically created from the simulation data; none is a photograph or experimental measurement.

## Empirical provenance

- Original microscopy and study: Lee et al. (2018), PNAS, DOI 10.1073/pnas.1720071115; raw binary images DOI 10.6084/m9.figshare.5969335.
- Derived area/count extraction: `https://github.com/squareshorts/colonization-history/tree/main/empirical`, `extract.py` and `occupancy.py`.
- The archived 8-series data table in `results/` is a transcription of the prior repository output `results/empirical/occupancy_validation.csv`, with reported original source commit and Git blob IDs in the separate technical-summary table.

The manuscript source is excluded from GitHub and distributed separately as an Overleaf ZIP. The existing repository DOI covers an earlier release; no new Zenodo release has been minted for this extension.
