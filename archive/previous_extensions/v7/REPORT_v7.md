# Limits of carrier attribution in bacterial colonization history

**Analysis date:** 6 November 2026 (Belém time)  
**Research status:** Executed computational and deposited-track analyses; biological mechanism attribution is not validated.

## Scope and provenance

This package adds an independent-field audit of deposited derived cell tracks and a calibration-aware partial-identification benchmark. It supplements the existing v3–v5 computational analyses. It neither re-extracts microscopy images nor validates an actual cell-state erasure, surface cleaning, or residual-cell viability assay. No new bacterial experiments were performed. The public repository was inspected at `squareshorts/colonization-history`, with baseline source commit `23cf4552cf2b1f2666e362b03fe20dd3fe015823`. The six `tracks_full.csv` Git blob identifiers are recorded in `results/field_level_v7/archival_field_summary_v7.json`.

The full field-level script `empirical/field_level_v7.py` reads existing CSV observables from that repository. In this analysis environment, source CSV records were independently read through the connected GitHub repository and aggregated. Aggregate values, without cell-level CSV redistribution, are preserved in `analysis/import_archival_field_v7.py` and the generated results. The field-level script has unit tests; its full end-to-end execution against the original on-disk repository requires a local checkout containing the original derived CSVs.

## 1. Empirical audit on the correct replication unit

Three paired flow-cell series in the published *Pseudomonas aeruginosa* transfer experiment were treated as three independent matched experimental units. Each cell track from the original derived `tracks_full.csv` was included when its `origin` was `arrival`, its initial observation was in the first three hours, and it was visible for at least two consecutive three-second frames. Kaplan–Meier survival to 300 seconds from first visibility was estimated within each flow cell, using the deposited censoring indicator. We repeated the analysis under 3-, 4-, 6- and 10-frame attachment thresholds and in initial 0–1 h and 1–3 h intervals. Differences use FC2 minus FC1. The 95% interval is the paired Student-t interval with two degrees of freedom. An exact eight-sign-flip enumeration is included as a descriptive sensitivity test; treatment assignment was not randomized within this archival reanalysis.

| Minimum frames | Mean FC2–FC1 five-minute survival | 95% field-level Student-t interval | Exploratory exact p |
|---:|---:|---:|---:|
| 2 | +0.0161 | [−0.3134, +0.3456] | 0.75 |
| 3 | +0.0286 | [−0.2749, +0.3321] | 0.75 |
| 4 | −0.0097 | [−0.482, +0.463] | 1.00 |
| 6 | −0.0265 | [−0.634, +0.581] | 1.00 |
| 10 | −0.0733 | [−0.786, +0.640] | 1.00 |

With the two-frame threshold, the three independent pair differences were +0.1362 (34 FC1 versus 606 FC2 attached tracks), +0.0385 (57 versus 275), and −0.1263 (9 versus 110). The 0–1 h mean was +0.1019 across pairs; the 1–3 h mean was −0.0616. The standard deviation across the three complete 0–3 h matched differences was 0.1326. Removing one matched series at a time gave mean contrasts of −0.0439, +0.0050 and +0.0873.

The analysis shows that the archival conditional residence contrast is highly sensitive to field identity, arrival window and attachment classification. It is not a claim that bacterial surface-sensing memory does not exist. The original 2018 study provides molecular and lineage evidence beyond these derived conditional residence endpoints. The present audit only limits conclusions from the specific reprocessed tracks. The assay includes a substantial imbalance in the number of arrivals that qualify for tracking in each flow cell and does not estimate attachment probability among all initial surface contacts.

## 2. Calibration-aware partial identification

An additive transformed-endpoint model represents the cell-associated history contrast as `C=A+S` and the surface-associated contrast as `F=E+R`, where `A` is acquired cellular state, `S` is selection of pre-existing adhesion phenotypes, `E` is persistent environmental material and `R` is viable residual-cell carryover. Two auxiliary readings have means `L_C=q_A A+q_S S` and `L_F=q_E E+q_R R`. The four factorial-arm means also identify a cell-by-surface interaction on the assumed endpoint scale, but no four-arm comparison alone distinguishes `A` from `S` or `E` from `R`.

For nonnegative effects and externally supplied intervals for each sensitivity coefficient, we form joint 95% intervals for `C,F,L_C,L_F` with Bonferroni-adjusted Student-t limits across independent field blocks. The feasible set for the cellular pair is the intersection of `A>=0`, `S>=0`, `C_low<=A+S<=C_high`, `q_A_low A+q_S_low S<=L_C_high`, and `q_A_high A+q_S_high S>=L_C_low`. Its sharp two-dimensional bounds are obtained by convex-polygon clipping and were checked against linear programming. The surface pair is analogous. Feasible bounds remain broad whenever the assay is insufficiently selective. The guarantee requires the additive model, the nonnegativity restriction and valid sensitivity intervals.

The synthetic observation model used matched blocks with a shared normal block effect (SD 0.28), independently generated factorial-arm errors (SD 0.26), and independent assay errors (SD 0.12). There were 16, 64 or 160 independent fields **per arm and per assay**. Each condition was repeated in **1,000 independent simulated datasets** under seed `20261106` using the documented fixed nominal and bounded coefficients. All numerical values are in synthetic transformed-endpoint units; none is fitted to bacterial microscopy.

**True coefficients in the main scenarios** were `q_A=0.78`, `q_S=0.18`, `q_E=0.14`, `q_R=0.67`, with assumed valid sensitivity bounds of `[0.70,0.85]`, `[0,0.20]`, `[0,0.18]`, and `[0.60,0.90]`, respectively. The deliberately miscalibrated point procedure assumed `q_A=0.78`, `q_S=0.075`, `q_E=0.065`, and `q_R=0.78` with no calibration uncertainty.

### False attribution when acquired and environmental effects are absent

For true `A=0`, `S=0.8`, `E=0`, and `R=0.7`, the fraction of trials declaring a strictly positive two-sided 95% lower bound was:

| Fields per arm / assay | False A, point calibration | False A, bounded calibration | False E, point calibration | False E, bounded calibration |
|---:|---:|---:|---:|---:|
| 16 | 73.1% | 0.0% | 14.7% | 0.0% |
| 64 | 99.9% | 0.0% | 47.7% | 0.0% |
| 160 | 100.0% | 0.0% | 85.8% | 0.0% |

The point method converges toward the incorrect pseudo-true parameters implied by a misspecified readout model. The bounded method preserves zero as feasible in these scenarios because the true coefficients lie inside the stipulated bounds. Those rates are conditional on a deliberately adversarial coefficient configuration and should not be treated as universal performance or biological false-positive rates.

### Power–identification tradeoff with all four carriers present

In the mixture `A=0.24`, `S=0.56`, `E=0.18` and `R=0.42`, the bounded method declared A positive in 99.9% of synthetic experiments with 64 fields, but declared E positive in only 4.0%. At 160 fields, E was resolved in 23.3% of experiments. Assay cross-reactivity bounds constrain how far replication can improve carrier identification.

### Failure under invalid calibration and model misspecification

A negative-control mixture with `A=E=0` placed the true selection sensitivity `q_S=0.28` above its supposed upper bound of 0.20, and placed the true carryover sensitivity `q_R=0.48` below its supposed lower bound of 0.60. Even the bounded method falsely declared A positive in 53.8% of simulations with 64 fields and 98.6% with 160 fields. This shows that the validity of calibration intervals is a prerequisite, not an optional robustness feature.

An additional fully mixed synthetic scenario introduced unmodeled within-side endpoint interactions and nonlinear assay responses, while retaining nominally valid first-order coefficients. At 64 fields the bounded cellular-A interval contained its stipulated first-order amplitude in 1.9% of simulated experiments. The method has no coverage guarantee under that altered response-generating process. The nonlinear benchmark is a deliberate falsification test of the additive approximation, not a fit to microbial occupancy.

### Cross-reactivity parameter sweep

A separate sweep varied the true selection cross-reactivity coefficient `q_S` from 0 to 0.30 in 13 steps with 32 or 128 fields per arm. Each grid cell contained 500 independent synthetic experiments, and true acquired cellular adaptation was absent. The conservative analysis stipulated `q_S` in `[0,0.20]` and the naive point analysis incorrectly fixed `q_S=0.075`.

At 128 fields, `q_S=0.150` produced false acquired-state declarations in **99.8%** of point-calibrated runs and **0.0%** of bounded runs. At `q_S=0.250`, outside the stipulated valid interval, the point method was positive in **100%** and the bounded procedure in **43.2%**. At `q_S=0.300`, the bounded false-declaration frequency rose to **99.4%**. The sweep confirms that the failure is systematic in cross-reactivity and is not restricted to the single main scenario. Monte Carlo resolution is 0.2 percentage points per grid cell; these conditional percentages do not estimate biological error frequencies.

Outputs are in `results/selectivity_sweep_v7/`, and the generating code is `analysis/selectivity_sweep_v7.py`.

## 3. Main scientific interpretation

The field-level archival reanalysis does not establish consistent independent-experiment residence enhancement under the selected conditional endpoint. The new simulation demonstrates a precise failure mechanism in causal localization, where increasing experimental replication increases confidence in a false carrier when assay selectivity is misspecified. Bounded sensitivity analysis avoids this under correct, sufficiently inclusive calibration intervals, at the cost of wider identified sets and sometimes low detection power. It can fail when either the calibration ranges or the assumed transformed-endpoint model are invalid.

The model-conditional result is complementary to the v5 nonlinear counterfactual analysis. The latter shows that incomplete carrier erasure can create false attributions in a stochastic colonization model. This v7 analysis supplies an explicit conservative inference procedure under known additive restrictions and tests its failure under broken calibration and response nonlinearities. The methods should be presented as different levels of abstraction; neither analysis validates the selectivity of a real biological intervention.

## 4. Hostile-review disposition after v7

**Computational reproducibility:** improved. Nine automated checks pass, including a comparison of the sharp polygon bounds with linear programming; the parameter sweep was also executed.

**Statistical reporting:** improved. Matched flow-cell experiments, tracking-definition sensitivity, field-level intervals, and partial-identification assumptions are explicit.

**Biological inference:** limited. No independent assays establish which carrier generated the published observations. The original microscopy has not been re-extracted or manually validated here. The new simulations are not prospectively calibrated against experimental instrument errors, cell-state reporters, or residual-cell survival.

**Novelty and likely editorial response:** the cross-reactivity failure and sharp bounded correction are operationally useful, but structural identifiability, set-valued inference and measurement-error bias have substantial precedents. A general methodological paper still needs a clearer comparative advantage against the closest published methods, preferably with independent biological validation. No honest computational audit can guarantee protection from desk rejection.

## 5. Reproduction and independent files

From a Python 3.10+ environment containing `numpy`, `scipy`, `pandas`, `matplotlib`, and `pytest`:

```bash
python -m pytest tests -q
python analysis/calibration_bounds_v7.py --replicates 1000 --n 16 64 160 --seed 20261106
python analysis/import_archival_field_v7.py
python analysis/selectivity_sweep_v7.py --replicates 500 --n 32 128
```

For full empirical rederivation, merge the package into the root of the original `squareshorts/colonization-history` checkout that contains the archived derived CSV tables, then run:

```bash
python empirical/field_level_v7.py --repo-root . --outdir results/field_level_v7
```

The `analysis/import_archival_field_v7.py` script uses the exact aggregate values independently calculated during this review from the connected GitHub blobs, preserving source hashes. It is not a substitute for verifying the underlying image segmentation or rerunning the full source-level empirical Python script on a local checkout. The repository and Zenodo release have not been modified.

## Primary literature

Lee et al. (2018). Multigenerational memory and adaptive adhesion in early bacterial biofilm communities. *PNAS* 115, 4471–4476. DOI 10.1073/pnas.1720071115.

Gajrani et al. (2025). Environment-mediated interactions cause an externalized and collective memory in bacteria. *The ISME Journal*. DOI 10.1093/ismejo/wraf173.

Liu, Maini and Baker (2026). Optimal experiment design for practical parameter identifiability and model discrimination. *Mathematical Biosciences* 399, 109710. DOI 10.1016/j.mbs.2026.109710.
