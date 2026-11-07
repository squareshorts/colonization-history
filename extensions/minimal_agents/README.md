# Synthetic minimal-agent analysis of colonization-history signatures

A complete, independently reproducible **exploratory** agent-based capability-ablation study, dated 7 November 2026. It is complementary to the user's preserved `colonization_history_v7_analyses.zip` and does not replace or change the original `squareshorts/colonization-history` GitHub repository.

## Start here

`REPORT.md` documents the experimental phenomenology, rules, ablations, bootstrapped comparisons, parameter sensitivity, oscillator test, novelty boundaries and publication assessment. The main figure is `figures/agent_capability_ablation.png`.

## Run on Windows PowerShell / any Python 3.10+ runtime

```powershell
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python src/run_capability_sweep.py --reps 48
python -c "import sys;sys.path.insert(0,'src');from run_robustness import run;run(10)"
python src/oscillatory_signatures.py
python src/run_surface_conditioning.py
python src/make_results.py
```

### Important constraints

- Published lag observations are from a fitted exponential curve; the agent analysis reports **first-passage time to 100 persistent agents** under censoring at 40 h. Numbers are not directly interchangeable.
- Parameters were selected with prior knowledge of published phenotypes. This is *not* a pre-registered prospective test, calibration to raw data, or a discovery of real bacterial mechanisms.
- Mean-restricted-time intervals quantify Monte Carlo field variability conditional on the synthetic parameters, not uncertainty across actual experiments.
- A two-state stochastic oscillator demonstrates how a sign-changing cAMP–TFP-type correlation can arise; it does not reproduce actual measured fluorescence.
- To reproduce the entire published microscopic phenomenon, independent experimental or raw-image validation is still required.

## Package contents

- `src/minimal_agents.py`: stochastic agents, 10 capability variants, hidden-state reset, daughter-transfer counterfactual.
- `src/run_capability_sweep.py`: main 48-field-per-condition comparison.
- `src/run_robustness.py`: 27 parameter combinations, 10 fields per arm.
- `src/oscillatory_signatures.py`: independent one- versus two-state stochastic time series.
- `src/run_surface_conditioning.py`: positive control with naive cells challenged on previously conditioned surfaces.
- `src/make_results.py`: field-bootstrap contrasts and figures.
- `tests/test_agents.py`: 11 deterministic checks.
- `results/`: raw synthetic-field CSVs and summary JSON/CSVs.
- `figures/`: publication-size 210-dpi figures.

Primary comparison: Lee et al. (2018), DOI 10.1073/pnas.1720071115. Related environmental memory: Gajrani et al. (2025), DOI 10.1093/ismejo/wraf173.
