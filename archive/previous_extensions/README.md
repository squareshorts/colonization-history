# Additional independently audited analyses (post-v1.1.0 repository update)

These additions complement the original repository and **do not change or refit** the v1.1.0 stochastic simulator or the Zenodo DOI associated with that release.

- `v7/`: independent matched-field derived-track audit and partial-identification Monte Carlo with calibrated sensitivity intervals, invalid-calibration controls and nonlinear misspecification tests. The included outputs were generated in a separately audited workspace. Execute `cd extensions/v7 && python -m pytest tests -q`. The scripts that regenerate derived-track summaries from the original deposited CSVs require the original repository root via `--repo-root ../..`.
- `minimal_agents/`: exploratory nonspatial, stochastic agent capability ablations and a separate oscillator-signature comparison. Execute `cd extensions/minimal_agents && python -m unittest discover -s tests -v` and use the README to reproduce the supplied CSV/JSON outputs.
- `identifiability/`: independent algebraic unit test for the six-arm ordered contrast identity.

All parameter values and simulated error distributions are illustrative, not estimated from the original microbiology images. The published 2018 paper contains original molecular and multigenerational evidence that is not contradicted by the restricted derived-track audit. New source and result files are not part of the earlier Zenodo v1.1.0 archive unless a new version is explicitly released.

Summary text results are included, while the large per-field observations and full agent JSON archives remain in the separate downloadable research package. Figures can be regenerated with the provided scripts and are included separately in the corresponding Overleaf package. The LaTeX manuscript is intentionally not committed to this public repository.
