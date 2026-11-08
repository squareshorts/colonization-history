# Surface-occupancy observation functions and holding-time profiles

The archival microscopy records support earlier sustained surface occupancy after transfer in three matched FC1/FC2 experiments. The absolute interpretation of the relative fivefold endpoint depends strongly on the number of initially attached cells. Early reference occupancies were 1.00, 3.07 and 2.29 cell-area equivalents in FC1, compared with 33.45, 12.55 and 9.05 in FC2. In FC1 some early arrivals detached or disappeared from view before the later rise in occupancy, making an initial-baseline fold increase an imperfect representation of continuous population growth.

All records were normalized to one 262-pixel reference (the median of three FC1-specific single-cell-area calibrations). Common sustained 20-, 30-, and 50-cell-area-equivalent thresholds were assessed using a ten-minute boundary-corrected moving average and fifteen minutes above the threshold. At the common 30-equivalent threshold, FC2-minus-FC1 differences were -24.8583, -12.2833 and -19.3250 hours. The direction was preserved across every specified threshold, smoothing, and persistence combination. In the FC2-1 record the 20-equivalent level was already maintained at the first available frame, so the time of first reaching that level is left-censored; the earliest observed attainment is 0 hours. In FC1 series 2, the 50-equivalent threshold remained unattained at the end of 26.7417 hours; FC2 series 2 attained it at 15.1083 hours, implying a transfer-associated advancement greater than 11.6333 hours under this readout.

The descriptive 14- and 37-hour normalized advantages depend on the measurement rule. The published fitted-lag point estimates give 0.939 and 0.322; the image-derived extrapolated occupancy lags give 1.126 and 0.299; the within-record fivefold timings give 1.305 and 0.700; the common 30-equivalent timings give 0.931 and 0.484. These are separate observational summaries and not direct measurements of a biological memory fraction. They come from the same archival experimental conditions, but there is only one recorded 14-hour and one recorded 37-hour holding condition available for an independent analysis. The ratios have no between-preparation confidence intervals.

The corrected analysis retains the original 20.24-hour fitted-lag estimate for FC1 series 2 as an explicitly extrapolation-dependent sensitivity outcome. Its growth-fit range extends only to approximately 33 cell-equivalents. The directly observed absolute threshold of 30 equivalents is attained at 26.425 hours within the recording and does not depend on an extrapolated growth curve. Every 30-equivalent threshold, including FC1 series 2, was attained during recording.

## Reproduction

A repository checkout with `empirical/observables` can regenerate the calculations with

```sh
python paper/analysis/observed_establishment.py --repo . --out paper/results/observed_establishment.csv --sensitivity paper/results/observed_establishment_sensitivity.csv
python paper/analysis/endpoint_observation_comparison.py --repo . --records paper/results/endpoint_observation_comparison.csv --summary paper/results/endpoint_observation_summary.json --sensitivity paper/results/absolute_threshold_sensitivity.csv
python paper/analysis/make_endpoint_comparison_figure.py --csv paper/results/endpoint_observation_comparison.csv --out paper/figures/endpoint_comparison.pdf
python -m unittest discover -s paper/tests -v
```

Data provenance: Lee et al. (2018), *PNAS*, doi:10.1073/pnas.1720071115; archival binary images, Figshare doi:10.6084/m9.figshare.5969335; existing derived `counts.csv` and `meta.json` files in `empirical/observables/` produced by `empirical/extract.py`. The current recalculation operates on previously segmented area data. No longer recordings were available, and no additional biological data were generated.
