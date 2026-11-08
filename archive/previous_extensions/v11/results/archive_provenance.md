The following derived observations originate from the public GitHub repository
https://github.com/squareshorts/colonization-history (main at commit
387e47a88a29e1259f1f2b66351f2204e7d6c903, last known at analysis start).

Input `archival_occupancy_validation_input.csv` is an exact transcription of
`results/empirical/occupancy_validation.csv`, including the three matched FC1-FC2
series and the single Tr14 / Tr37 recordings. The derived-record fit uses the
repository `empirical/occupancy.py` estimator, which is NOT the identical
published lag estimator. The source `counts.csv` blobs and their Git object SHAs
appear in `technical_fit_summary.csv`. Counts and per-dataset `meta.json` were
read through the connected GitHub service for the within-record smoothing and
fit-window sensitivity calculation, without copying or reprocessing raw images.
The bundled recomputation script reproduces the calculation from a repository
checkout and emits the complete fit grid. The summary values were rounded to
two decimals and should not be mistaken for independent biological uncertainty.

Both prolonged-hold groups have one run in the distributed observables. A
between-replicate confidence interval for their washout-lag contrast is not
estimable from these materials. The published fitted curve CIs are different.
