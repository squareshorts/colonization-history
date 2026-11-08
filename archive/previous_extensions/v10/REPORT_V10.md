# v10 audit: general scalar-output equivalence, recovery bounds and localization

7 November 2026. Computational and mathematical study; **no new bacterial experiments**. v9 remains an independent archived computational phase. v10 scripts and results are reproducible with NumPy, SciPy and Matplotlib, and use a new random seed series.

## Questions answered

1. Does exact equivalence depend on the particular Hill/logistic construction? **At the level of an unconstrained scalar readout, no.** Any injective trajectory x(t) can represent an arbitrary target Y(t) on the observation interval via g(x)=Y(x^{-1}(x)). This is a familiar input-output reparameterization, not a new theorem and not a claim of equivalence for molecular or lineage measurements. Fixing g by independent biological calibration can break it.
2. Can selective loss create spurious mother/daughter inheritance signals? **Yes, but sharp missing-outcome bounds can prevent unwarranted attribution when each founder and one daughter are individually indexed.** They protect the assigned-cell estimand at the price of broad intervals under severe attrition. Sampling inference still relies on independent clone blocks and an appropriate uncertainty model.
3. Can the five-carrier system (intracellular state A, selection S, cell-associated exterior C, surface conditioning E, viable carryover R) be mathematically separated? **Yes, conditionally.** A fifth localization contrast V=gamma*A+kappa*C closes the rank deficiency only when gamma != kappa and their differential sensitivities are independently calibrated. No such assay was experimentally calibrated here.

## Generality audit

Five independently selected synthetic lag families: Weibull, stretched exponential, biexponential, logistic and delayed onset. In 300 independent reparameterization trials, each had 25 time points to train two flexible latent readouts and 600 unseen holding times to evaluate (180,000 test times total). Mean holdout RMSE was 0.01028 synthetic hours for the decaying-state coordinate and 0.01190 for the selection-frequency coordinate. Maximum disagreement between the two fitted readouts across tests was 0.13717 synthetic hours. These are numerical interpolation results using constructed data. Exact equality is a mathematical identity only with unrestricted response functions.

## Recovery-bias correction and independent microscopic observation laws

For N initially assigned and individually indexed cells, n recovered, S observed binary attachment successes, the **sharp finite-population** interval for the original population mean is [S/N, (S+N-n)/N]. For a pulsed arm 1 and sham arm 0, the originally assigned-population contrast is in [S1/N1 - (S0+N0-n0)/N0, (S1+N1-n1)/N1 - S0/N0]. The departure of recovered-arm mean from original-arm mean is at most eps=1-n/N; the differential recovery component is bounded by eps1+eps0. This is an established Manski/Horowitz-style missing-outcome bound, not a new estimator. The bound does NOT eliminate random group imbalance, unmeasured differential growth, unknown founder identity, or cell-associated extracellular carryover.

Each kinetic law used 1,500 independent simulated studies; each study included 32 independent clonal blocks, with 48 originally assigned founders in each pulse/sham arm. Every recovered founder could contribute only one indexed immediate daughter; this one-to-one condition is required for the daughter bounds. Inference uses one-sided paired t contrasts with Bonferroni alpha=0.025 per mother/daughter endpoint. These intervals are approximate at finite clone count, although the missing-outcome identified sets are algebraically sharp.

| Model | Exponential dwell: daughter positive | Weibull dwell: daughter positive | Exponential dwell: joint recovered-only | Weibull dwell: joint recovered-only | Joint attrition-bounded, both laws |
|---|---:|---:|---:|---:|---:|
| Inherited acquired state, complete recovery | 96.9% | 97.8% | 96.9% | 97.8% | same as recovered-only |
| Acquired state, no daughter inheritance | 3.6% | 2.8% | 3.6% | 2.8% | same |
| Pure selection, complete recovery | 2.5% | 2.8% | 0.07% | 0% | same |
| Phenotype-biased recovery: high-adhesive survival .98, low-adhesive survival .45 | 26.0% | 31.6% | 10.5% | 16.5% | 0% |
| Phenotype-biased recovery: .998 vs .98 | 3.3% | 2.9% | 0.13% | 0% | 0% |
| Inherited exterior coating | 96.9% | 98.5% | 96.9% | 98.5% | same |
| Induced adhesive phenotype switching | 96.9% | 98.6% | 96.9% | 98.6% | same |

A high-attrition preparation recovered about 60% of initially assigned pulsed cells on average, while the lower-loss preparation retained about 98.5%. **No conclusion about cell-intrinsic molecular localization follows from positive mother/daughter effects.** The induced-switch comparator shows that even a phenotype change without the postulated biochemical memory carrier can produce a similar observation.

## Exterior-carrier localization and specificity calibration

The five-carrier matrix is:

    [1   1   1   0   0]    cell-history transfer: A+S+C
    [0   0   0   1   1]    conditioned-surface challenge: E+R
    [1   0   1   0   0]    unsorted clonal pulse: A+C
    [0   0   0   1   0]    viable-cleared surface: E
    [g   0   k   0   0]    molecular localization: g*A+k*C

Full rank requires g != k under these ideal additive rows. The condition number is 4.05 at (g,k)=(1,0), 48.58 at (1,.9), 250.45 at (1,.98), and the matrix has rank 4 at (1,1). With n=32 synthetic fields, per-field SD .12, underlying A=.25, C=.18, and **valid** gamma in [.95,1.0] and kappa in [0,.2], simultaneous Gaussian observation intervals gave conservative identified bounds for the expected observation: A in [.154,.331] and C in [.052,.323]. Across 4,000 Gaussian datasets with valid specificity, intervals covered both carriers in 100%, illustrating conservatism. With invalid actual kappa=.65 the coverage dropped to 12.8%; with invalid actual gamma=.60 it dropped to 14.9%. These reflect simulation-specific ranges, not measurements of an actual exterior-specific assay.

## Hostile-review verdict

The result now has three defensible components: a concrete all-time Hill/logistic equality, a general reparameterization explanation of why a single latent readout is insufficient, and an operationally precise intervention benchmark with a sharp attrition correction and localization-specific rank condition. The general reparameterization and missing-data bounds are prior mathematical results (Meshkat et al., 2011; Horowitz and Manski, 2000; Manski, 2005). Novelty should be claimed for the integration into the bacterial colonization-history design and its explicit counterexamples / failure controls, not for inventing these theorems.

Unresolved criticisms: (a) no experimentally calibrated parameters or assays; (b) published growth lags are not the same estimator as synthetic first-passage times; (c) surface/coating and induced phenotypic switching confound a particular intracellular carrier; (d) valid inference after missing daughters requires one-to-one founder lineage enumeration and no selective differential cell reproduction; (e) v7 archived empirical audit has only three matched series and the original raw images were not reprocessed; (f) cross-model holdout test uses *synthetic* time points and flexible readouts, so it is not biological validation. These should remain explicit in the manuscript.

## Run

    python analysis/generality_attrition_v10.py --trials 300 --studies 1500
    python -m unittest discover -s tests -v

`results/v10_summary.json` contains complete results and metadata; CSV files retain full synthetic replicates and kinetic summary rows. Figures are regenerable. Full earlier v7, v8 and v9 methods remain in their own repository extensions and companion archives.

### Primary references

- Horowitz JL, Manski CF (2000). Nonparametric Analysis of Randomized Experiments with Missing Covariate and Outcome Data. JASA 95:77-84. DOI:10.1080/01621459.2000.10473902.
- Manski CF (2005). Partial identification with missing data: concepts and findings. International Journal of Approximate Reasoning 39:151-165. DOI:10.1016/j.ijar.2004.10.006.
- Meshkat N, Anderson C, DiStefano JJ (2011). Finding identifiable parameter combinations in nonlinear ODE models and the rational reparameterization of their input-output equations. Mathematical Biosciences 233:19-31. DOI:10.1016/j.mbs.2011.06.001.
- Lee CK et al. (2018). Multigenerational memory and adaptive adhesion in early bacterial biofilm communities. PNAS 115:4471-4476. DOI:10.1073/pnas.1720071115.
