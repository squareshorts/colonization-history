# Corrected archival occupancy audit

The 20-observation moving average uses a denominator equal to the number of available samples at the beginning and end of each time series. This removes the zero-padding edge artifact in earlier `np.convolve(..., mode="same")` estimates.

Record-specific paired mean FC2-FC1 = -10.6908 hours, descriptive t interval [np.float64(-18.263805920827494), np.float64(-3.117714208288196)]. Shared-area paired mean = -11.6648 hours, descriptive t interval [np.float64(-17.370429594290734), np.float64(-5.959202589852172)].

Corrected normalized advantages, record-specific: 14 h 1.1257, 37 h 0.2992. Shared: 14 h 1.0420, 37 h 0.2860. One archived recording exists for each prolonged hold.

FC3 (series 3) lag: 0.9445 hours, shared calibration 0.9359 hours. FC1 series 2 uses 152 half-minute samples in the fit and extrapolates from about 20-33 cell equivalents.

The derived measurements were recalculated from the GitHub-hosted `empirical/observables/<dataset>/counts.csv` and `meta.json` files on 2026-11-08. The script in the reproducibility package recomputes all values from an ordinary local repository checkout. No raw microscopy segmentation was repeated.
