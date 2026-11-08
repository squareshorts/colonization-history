# Directly observed establishment time — archival endpoint revision

The binary-image experiment `WT YFP-CFP (T=26.75h;FC1)` ends at 26.75 h and reaches only approximately 33 occupied cell equivalents. Its fitted exponential lag extrapolates from a narrow part of the recorded expansion. These constraints cannot be changed by additional curve fitting or by simulation. A longer recording of **the same historical experiment** is unavailable in the deposited series. No new biological measurements were generated.

## Main endpoint

Time to first observed sustained fivefold occupancy reference, using the edge-corrected 20-observation (~10-minute) smoother, the median of the first two hours as a within-record reference, a one-cell-equivalent minimum baseline, search starting after the two-hour reference period, and a 15-minute persistence requirement. All series crossed within their available recording.

| Series | FC1 (hours) | FC2 (hours) | Difference FC2−FC1 |
|---|---:|---:|---:|
| 1 | 21.0167 | 11.0917 | −9.9250 |
| 2 | 24.3667 | 16.7917 | −7.5750 |
| 3 | 23.6583 | 10.0667 | −13.5917 |

The mean is −10.3639 hours; descriptive paired 95% Student t interval [−17.8964,−2.8314] h, based on only three matched experimental series. A second consecutive transfer (FC3) in series 3 crossed at 4.05 h. The 14- and 37-hour held-cell recordings crossed at 9.4917 and 15.7583 h, respectively, but each has only one archived run and varying early occupancy.

## Sensitivity and interpretation

Retrospective threefold and fivefold threshold checks with 5-, 10- and 20-minute moving averages, and zero-, 15- and 30-minute persistence windows, preserve the sign of all three paired differences. At the fivefold threshold, shared 262-pixel calibration gives mean −10.3778 h; removing the one-equivalent reference floor gives −10.3194 h.

The threshold is an operational trajectory summary. It is sensitive to initial delivery, attached-cell number, component segmentation and subsequent growth. The endpoint does not estimate a single-cell memory or the exponential growth lag of Lee et al. (2018). The threshold was selected retrospectively after inspecting the available data, and the small-sample interval is descriptive. The original exponential-lag estimate for FC1 series 2 remains uncertain; a reliable value requires a longer biological recording or independent biological replication with a sufficiently extended growth range.

## Exact reproducibility

Source image deposit: Lee et al. 2018, Figshare DOI 10.6084/m9.figshare.5969335. Derived occupancy inputs are `empirical/observables/{record}/counts.csv` and `meta.json` in `squareshorts/colonization-history`. These inputs were previously produced by `empirical/extract.py`; they are not supplied as precomputed counts in the original imaging deposit.

The script `analysis/observed_establishment.py` recomputes the endpoint from the original derived records in a checked-out repository, including the full threshold/smoothing/persistence sweep and the descriptive paired mean. `analysis/make_observed_figure.py` rebuilds the primary plot from the CSV. The included CSV files reproduce independently checked results. The earlier fitted-lag sensitivity remains in the supplementary article and unchanged older archive.

No Zenodo DOI was minted or changed.
