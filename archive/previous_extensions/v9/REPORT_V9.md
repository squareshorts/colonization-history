# Exact observational equivalence and discriminating-intervention test

7 November 2026. This is a **completed theoretical and synthetic validation**. It contains no new bacterial measurements. Parameters and simulated effect sizes are deliberately uncalibrated. Original biological observations are from Lee et al., *PNAS* 115:4471-4476 (2018), DOI 10.1073/pnas.1720071115.

## Central result

A simple population-selection model and a decaying single-cell state with a nonlinear Hill readout can produce **identical lag trajectories at every holding time**, not merely at the three observed holding times. Specifically,

    selection fraction       f(t) = [1 + q exp(c t)]^-1
    decaying state           m(t) = m0 exp(-t/tau)
    Hill readout             h(m) = [1 + (K/m)^nu]^-1
    equivalent parameters   nu/tau = c,  (K/m0)^nu = q
    observable advantage    D(t) = D(0) f(t)/f(0)
                                  = D(0) h(m(t))/h(m0)

In both cases the macroscopic growth lag is `L(t) = L_naive - D(t)`. The two latent mechanisms carry fundamentally different biological interpretations. Yet their lag likelihoods are identical whenever the observation noise laws are identical. No number of holding-time points, nor arbitrary additional replication on the same lag endpoint, resolves this equivalence within this model class.

Using the reported growth-lag point values (naive 17.56 h, transfer 0 h 4.68 h, 14 h 5.46 h, 37 h 13.41 h), one convenient representation has selected-phenotype initial fraction 0.9903590, selection relative growth disadvantage 0.1456747 per hour; the identical intracellular representation has chosen memory time constant 10 h, Hill exponent ~1.45675, and state midpoint obtained analytically. The exact numerical gap across a dense 0-60 h grid is below 4e-15 h. These parameters merely interpolate known observations and have no confirmed biological realism.

An additional parameter-global stress test drew 1,000 independent combinations of phenotype fraction (0.05--0.95), relative growth disadvantage (0.04--0.25 per hour), intracellular lifetime (6--25 h) and initial state (0.4--0.95), without using Lee et al. points. The two readouts remained equal at all 121 tested holding times up to numerical error of 4.44e-16 in frequency; the identity is not a special fit to the three observed lags.

A new lag-only logistic classifier, trained and evaluated with separate simulated datasets and independent Gaussian lag measurement noise of 0.75 h, achieves ROC AUC 0.522, near the analytical chance value of 0.5. Its exact chance property follows from the identical likelihoods, not from classifier performance.

## Discriminating prospective intervention

The proposed challenge randomizes transient contact versus sham **within the same clonal lineage**, verifies complete recovery and unchanged phenotype composition, and then challenges equal numbers of cells on virgin substrates before reacquisition of any surface history. Both treated mothers and immediately isolated daughter cells are tested separately. The pre-specified primary endpoint is the fraction of newly introduced cells remaining attached for five minutes. Each independently prepared clone is one paired biological block. The corresponding daughter assay is a test of parent-to-daughter transmission; a nonsignificant daughter result alone does not establish absence of inheritance.

A separately implemented microscopic simulator computes stochastic contact success and exponential detachment survival, marginalizing over pre-existing adhesive phenotypes within each clone. It uses independent clone-specific variation and broad uncalibrated nuisance-parameter distributions; it is **not** a simulation of the published 2018 microscopy or a fit to the colony lag curves. For each of 2,000 simulated studies in each condition, pulse/sham differences are tested across independent clone blocks with one-sided paired Student-t tests at alpha=0.025 per endpoint (Bonferroni for mothers and daughters).

At 32 clonal blocks, 48 cells per arm per block:

| Model | Positive mother response | Positive daughter response | Both | Mean true mother/daughter risk differences, pp |
|---|---:|---:|---:|---|
| Induced inherited intracellular state | 100.0% | 98.4% | 98.4% | 12.39 / 10.15 |
| Induced state with no daughter inheritance | 99.9% | 2.55% | 2.55% | 12.40 / 0.00 |
| Immutable phenotype selection with complete recovery | 3.0% | 2.55% | 0.05% | 0.00 / 0.00 |
| Inherited extracellular coating attached to the cell | 99.95% | 98.45% | 98.4% | 12.13 / 9.83 |
| Selection with phenotypically differential preparation recovery (odds ratio 8) | 87.8% | 89.1% | 79.85% | 3.95 / 3.95 |

**Interpretation.** The randomized full-recovery protocol can distinguish an induced acquired effect from *pure selection of pre-existing states*. It can also provide evidence consistent with transmission to daughters. The test **cannot distinguish intracellular state from inherited extracellular material on the cell**. The coating alternative is a deliberate observationally equivalent negative control, motivated by documented cell-associated extracellular polymers in *P. aeruginosa* (Ma et al., 2009, DOI 10.1371/journal.ppat.1000354). It is not claimed to explain the Lee et al. PA14 data specifically.

## Robustness and failure controls

Across preparation survival odds ratio values 1, 1.5, 2, 4, 8, and 16 in a population with pre-existing phenotypic heterogeneity, the fraction of false *joint* inherited-acquired declarations rises from approximately 0% to 1.1%, 8.7%, 50.7%, 79.0%, and 90.0%. This stress test establishes that **verification of full, phenotype-neutral cell recovery is part of the identifying intervention**. If preparation enriches adhesive cells, a within-clone contrast can be selection-derived, even when the clone has one genotype.

At 32 clonal blocks, weakening the induced response amplitude to 0.25, 0.50, 0.75, 1.0 and 1.5 times the stipulated baseline gives joint detection fractions 33.7%, 83.2%, 95.2%, 97.8% and 99.7%, respectively. Holding induced cells in suspension for 0, 7, 14, and 37 synthetic hours yields 97.8%, 76.3%, 42.1%, and 1.7% joint detection under this independent kinetic model. These are study-design sensitivity results, **not laboratory power estimates**. The true effect and cluster variance must be measured in biological pilot work.

## Minimal carrier discriminability and assay selectivity

For the additive carrier vector `(A,S,E,R)`, standard four cell/surface factorial arms span only two distinct directions, corresponding to `A+S` and `E+R` (rank 2). Adding a randomized no-selection cell pulse measuring `A` and a completely viable-free, trace-preserving conditioned-surface challenge measuring `E` increases the rank to 4. Each added arm is necessary in this restricted design family. This is a rank result for an **additive, calibrated contrast scale**, not a statement about realistic causal identification without intervention validation.

With imperfect pulse selectivity, its observed contrast is `A + kappa*S`; with imperfect surface conditioning it is `q*E + rho*R`. Identifiability requires `kappa != 1` and `q != rho`; the condition number grows rapidly as they approach equality. At 32 illustrative independent fields and a Gaussian contrast SD of 0.12, the full matrix has condition number 2.62 in the ideal case, 18.14 for `(kappa,q,rho)=(0.8,0.75,0.25)`, and 398.0 for `(0.99,0.6,0.5)`. Near-singularity destroys practical precision even if the matrix remains formally full rank.

Adding a fifth carrier, inherited extracellular material attached to cells (`C`), expands the causal ambiguity. Standard cell transfer measures `A+S+C` and the within-clone mother/daughter pulse measures `A+C`. The proposed four-carrier assay has rank 4 for five unknown components; its null direction exchanges `A` and `C`. An additional assay with **independently validated molecular-location selectivity** is necessary to distinguish those carriers. This is not currently available in the analyzed data.

## Independent verification, novelty and limits

Twelve regression and algebra tests pass, including exact equality across held-out hold times, rank deficiency, out-of-model coating counterexample, differential-recovery false positives, independent event-time simulation of the microscopic attachment probability, washout of the proposed mother/daughter contrast, and seed reproducibility.

The novel synthesis is a constructive observational-equivalence proof **plus a smallest-in-class intervention design and an independently coded finite-sample stress test with explicit failure modes**. The underlying Hill/logistic identity and linear rank arguments are elementary established mathematics, and novelty should not be claimed for those mathematical tools alone. The scientific insight is that even arbitrarily well-measured lag trajectories can be unable to choose between mechanistic carriers, and that the nature and validity of additional experiments matter more than adding temporal resolution to the same endpoint.

The result is *not* biological confirmation. The same-cell signaling and molecular-lineage evidence in Lee et al. (2018) provides information beyond the coarse lags and is not nullified by this theorem. This synthetic experiment is prospective **as an experimental design**, not a completed prospective biological validation; no independent laboratory test has been performed.

## Reproduce

From a Python environment with numpy, scipy, scikit-learn and matplotlib:

    python analysis/discriminating_test_v9.py --datasets 2000 --cells-per-arm-clone 48 --seed 202611081
    python analysis/finalize_v9.py
    python analysis/global_equivalence_stress.py
    python -m unittest discover -s tests -v

The outputs are `results/critical_test_v9.json`, `results/intervention_power.csv`, `results/assay_design_condition_numbers.csv` and four figures. The internal original-model and v7/v8 results remain unchanged.
