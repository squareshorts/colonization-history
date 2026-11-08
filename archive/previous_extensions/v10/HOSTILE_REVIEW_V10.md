# Hostile review after v10

## Verdict

The paper has a defensible, more general **observation-specific identification** result plus a statistically honest intervention audit. It still cannot claim that real bacterial intracellular memory is uniquely identified. I would submit as a *theoretical-biology methods* article, not as a mechanistic biological discovery.

### Reviewer objection: The Hill/logistic equivalence is cherry-picked.

**Corrective analysis completed.** A general injective-latent/readout reparameterization proof shows that any scalar lag trajectory can be represented by either latent process on a specified observation interval when its readout is unrestricted. A numerical interpolation check tested five different target curve families with 25 training and 600 disjoint test times, repeated 300 times. It achieved mean RMSE ~0.01 simulated hours for either latent coordinate. The theorem itself is standard identifiability mathematics, and the data are entirely synthetic. This result **does not** validate equivalence of cell-level signaling, physical attachment models with independently constrained response laws or multioutput models.

### Reviewer objection: Randomized clonal assays cannot exclude phenotypic selection when recovery differs.

**Corrective analysis completed.** For individually indexed founder cells with binary outcomes, we implemented sharp missing-outcome bounds accounting for the fraction lost in each treatment arm. We tested arbitrary missing-outcome patterns and sharpness. Severe phenotype-biased recovery produced 10.5% and 16.5% false joint detections with a recovered-only paired test under two distinct kinetic laws; the attrition-robust joint procedure made 0/1500 false declarations in each simulated batch. The procedure is conservative under severe loss, and cluster t intervals remain approximations. A real implementation requires documented founder and daughter identities, recovery denominators and no untracked differential proliferation. Merely reporting 99% recovery does not eliminate finite-sample baseline imbalance.

### Reviewer objection: An inherited extracellular coating mimics an intracellular state.

**Boundary proven; not biologically solved.** The five-carrier additive design has rank four until an independent measurement has different sensitivity to intracellular state and exterior carried material. An ideal specificity contrast produces full rank and moderate conditioning. Near-identical sensitivities produce extreme error amplification. With valid stipulated sensitivity intervals the synthetic procedure was conservative, but deliberately invalid calibration reduced coverage of true carriers to ~13-15%. No real location-specific assay was validated.

### Reviewer objection: Are these experimental predictions independently validated?

**No.** The two kinetic laws and five time-shape holdouts are independently specified *simulations*, not biological tests. A future laboratory study would need to measure phenotype-neutral recovery, one-to-one mother/daughter identity, true assay specificity, five-minute attachment and the corresponding molecular signatures. There is no additional experimental dataset supporting a direct carrier identification in this manuscript.

### Reviewer objection: Are the mathematical tools new?

**No.** The input-output change-of-variable identity and missing-data identified bounds are established mathematical ideas (Meshkat et al., 2011; Horowitz and Manski, 2000; Manski, 2005). The manuscript's incremental novelty is their integration with bacterial colonization-history simulations, an exact interpretable counterexample, a specified 5-carrier probe design and calibrated failure tests. Novelty is methodological synthesis and application, not a new abstract theorem.

## Submission assessment

The manuscript is improved as a theoretical identifiability and experimental-design study. If submitted to *Theory in Biosciences*, the abstract, title and cover letter should emphasize conditional model discrimination and the precise necessary calibration requirements. No amount of further simulation can replace biological validation of assay selectivity. All claims must remain conditional on the observation model and tracking design. The report does not represent the synthetic detection percentages as laboratory power estimates.
