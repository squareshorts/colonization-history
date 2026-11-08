# Scientific and numerical audit of carrier-discrimination analyses

## Independent replication of the apparently duplicated 34.3%

Three distinct quantities rounded to the same percentage in two separate analyses. The complementary-log-log/Weibull law with half-rate sham loss and 10% pulse loss produced 34.333% bounded joint detection in 1,200 studies from the coarse loss grid. The fine-grid replication of the same condition gave 30.4375% in 1,600 independently seeded studies. An additional four-batch, 32,000-study replication of exactly this design yielded 32.8125% (normal-approximation Monte Carlo SE ~0.262 percentage points). The calibrated-binomial-tail test in a separate fixed-probability population with equal 10% loss gave 34.2857% in 2,800 studies. Coincidental rounded agreement did not reflect a duplicated result.

## Recovery design rule

The expected lower bound for a randomized pulse-minus-sham contrast under independent uniform attrition is

E[ell]=(1-epsilon1)*p1-p0-epsilon0*(1-p0).

For equal attrition, E[ell]>0 precisely when epsilon<(p1-p0)/(1+p1-p0). The logistic/exponential model with e=0.85 gives approximate loss limits of 6.09% for daughter and 9.15% for mother attachment; the daughter is limiting for joint detection. These are thresholds for an expected contrast, not finite-sample detection thresholds.

## External-probability calibration and risk

The main kinetic model uses heterogeneous block-level predictors and allows some high-type sham probabilities to exceed 0.10. A universal bound u0=0.10 would be invalid for that population. The probability-assisted demonstration was instead conducted in a separate synthetic population with p_high=0.10 and p_low=0.03, under the additional assumptions of independent conditional Bernoulli outcomes and independently known u0>=max_i p0i.

Under invalid u0=0.03 and equal 10% loss, 2,800 null studies produced zero joint false declarations while the empirical missing-sham tail coverage was 91.16% (below the nominal 99.5%). Invalid coverage does not imply that every tested null scenario has a large false-positive rate.

A new independent null stress test used 6,000 studies per condition, 32 independently sampled clone blocks, 48 founders per arm per block, and a high-type fraction sampled between 0.80 and 0.99. Under 40% sham-only loss, u0=0.03 produced 929/6000 = 15.483% false joint attributions and effectively zero mother missing-success upper-tail coverage. Correct u0=0.10 yielded 0/6000 false joint attributions. Under 30% and 20% sham-only loss, the invalid cap produced 3.667% and 0.850% false attributions. The simulation validates an identified failure mode, without claiming that the conservative test has a universal finite-sample guarantee.

## Scientific and document checks

The revised main article includes a 175-word abstract; both loss thresholds, including equal-arm loss; a precise distinction between missing-cell tail coverage and false attribution; a clarification of separate calibration populations; and an Imbens-Manski citation adjacent to partial-identification inference. The supplement states the full 96-design factorization (8 loss levels x 2 kinetic laws x 3 allocations x 2 effects), provides the 32,000-study verification, and tabulates adversarial null false attribution. All figure paths and SI section references are checked. Internal version and development references are absent from the article and supplementary scientific prose.

## Outstanding archiving and authorization

The current code revision needs a separately archived, version-specific DOI for the Royal Society submission route. No new DOI has been minted, respecting the author's standing instruction. Funding, CRediT attribution and AI-use declarations require final author confirmation under the target journal's policy.
