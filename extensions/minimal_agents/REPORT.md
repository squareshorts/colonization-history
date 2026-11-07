# Minimal agent capabilities for colonization-history signatures

**Exploratory simulation and capability ablation — 7 November 2026**

## Question and scope

Which small sets of agent-level capabilities are sufficient to reproduce the *qualitative* population-history signatures reported for *Pseudomonas aeruginosa* surface colonization? The purpose is to identify transparent generative counterexamples and discriminating measurements. Model sufficiency is not identification of the true biological carrier. This is an exploratory, constructed example, not a fit to raw experimental observations.

The completed independent v7 work is preserved separately as `colonization_history_v7_analyses.zip` in the `colonization-history` Library folder. This agent study leaves that package and the public GitHub repository untouched.

## Reference observations

Lee et al. (2018), *PNAS* 115:4471–4476, DOI 10.1073/pnas.1720071115, reported the following:

1. About 95% of early naive surface contacts have residence times in the 3–30 s range, far below the approximately 1-h cell division time.
2. Following transfer to a virgin flow cell, previously exposed cells and their descendants display substantially earlier establishment and more multigenerational family structures than naive cells.
3. Their published exponential-fit colonization lags were approximately 17.56 h (naive FC1), 4.68 h (immediate transfer), 5.46 h (14-h hold), and 13.41 h (37-h hold). The reported CIs refer to the original fitted trajectories, not independent field-level uncertainty.
4. Previously exposed lineages exhibit damped, out-of-phase cAMP and type IV pilus activity correlations, including a sign change around a 5-h lag.

These are separate *observable targets*. The archived v7 reanalysis independently found heterogeneous early residence effects across only three matched tandem experimental series. The published lag observations and the v7 field-level contrast must not be merged into an artificial claim of uniform replication.

## Synthetic agent grammar

The simulator is a **nonspatial**, well-mixed, capacity-limited, stochastic agent model. Attached resident agents are represented by a continuous memory state `m in [0,1]`, founder identity and lineage generation; fixed adhesive phenotype `z` and a shared surface trace `S` can be switched on as comparators. Agents arrive according to a Poisson process, transiently contact the surface, may establish persistent residence, divide, retain daughters, or detach. Short unsuccessful visits receive exponential residence times with a 10-s illustrative mean.

A contact increments intracellular state by `m <- m + a (1-m)`, subject to variation in each agent's number of preparation contacts. Free-cell state decays as `m(t)=m(0) exp(-t/tau)` (default tau = 10 h). With ongoing sensing, resident state increases during surface residence. Both persistent attachment probability and daughter retention depend on a normalized Hill function of `m`. Newborns optionally inherit their parent's state; newborns with inherited state disabled can still reacquire state on the surface. The environment-only comparator includes surface-state accumulation and decay, but all transfer challenges occur on **fresh surfaces**. Immutable-selection agents have fixed adhesive types that can be enriched by prior exposure; fixed types are passed to descendants and do not relax during suspension.

The synthetic dynamics are deliberately simplified. They do not track fluid flow, Psl deposition, type IV pili, actual cAMP concentrations, growth in holding tubes, or detailed microscopic image formation. The transferred donor population is represented by individual stochastic numbers of prior contacts, not simulated from a validated FC1-to-FC2 flow-cell passage.

### Prespecified synthetic endpoint and limitations

The primary endpoint is the **first passage to 100 persistent attached agents**, with a 40-h horizon. Observations that fail to reach 100 are assigned 40 h for the reported restricted-time mean, and the crossing fraction is reported separately. The original paper estimated lag from fitted exponential growth curves. **These two lag definitions are not interchangeable.** Direct numerical closeness between them has no inferential standing.

Parameters were selected illustratively after consulting the known experimental phenomenology. The displayed 0/14/37-h comparisons are **not an external validation or a blinded held-out test**. No parameter has been fitted by maximum likelihood to raw microscopy. Every CI below is conditional on this synthetic simulator and its hand-set parameters.

## 1. Capability ablation

Each mechanism/arm has 48 independently generated stochastic fields. Common seeds are used across model families for reproducibility and partial variance control, but simulation pathways change the number of random variates consumed, so these are not exact paired potential-outcome trajectories. Within-model challenge arms use independent stochastic streams.

**Restricted time (hours) to the operational threshold**

| Agent capabilities | Fresh | Transfer, 0-h hold | Transfer, 14-h hold | Transfer, 37-h hold |
|---|---:|---:|---:|---:|
| Passive | 25.23 | 27.11 | 26.21 | 28.56 |
| Immutable selection | 22.13 | 5.88 | 5.98 | 5.80 |
| Environmental state on a fresh challenge surface | 25.20 | 27.14 | 26.22 | 28.54 |
| One-bit stochastic memory | 20.45 | 5.81 | 8.19 | 14.76 |
| Decaying scalar memory, no daughter inheritance | 19.29 | 5.52 | 5.62 | 19.27 |
| Decaying scalar memory, daughter inheritance | 18.99 | 5.32 | 5.46 | 18.98 |
| Scalar memory without decay | 18.99 | 5.32 | 5.39 | 5.67 |
| Scalar memory with linear feedback | 20.12 | 5.63 | 8.50 | 16.43 |

The simulation recapitulates three qualitative signatures with a **single decaying intracellular scalar and a nonlinear response**. Transferred cells colonize earlier; a short hold does little; a long hold removes much of the advantage. The naive short-contact proportion in the scalar model is about 94.5%, corresponding to the independently stipulated exponential short-contact rule. No mechanistic conclusion about actual bacteria follows from these choices.

A direct positive control verifies that environmental state can create a strong history effect under the appropriate **same-surface/conditioned-surface** design. In 60 independently seeded environment-only fields, naive cells on virgin surfaces had a mean restricted threshold time of **24.10 h** and a mean of **4.3** cells at 8 h; naive cells on conditioned surfaces had a mean restricted threshold time of **8.42 h** and **99.9** cells at 8 h. The environmental state is not transmitted with cells to a virgin challenge surface. These numbers are synthetic and depend on the illustrative initial surface state (`S=0.70`) and its assumed decay. See `results/environment_positive_control.csv`.

A fixed-selection comparator can produce rapid transfer and large lineage growth but retains its advantage after 37 h; the no-decay cellular comparator fails similarly. These failures reflect *their specified immutable assumptions*. A selection model with phenotypic turnover, population dilution, or changing selection pressures might also yield the washout, and would require its own comparison. An environment-only state cannot influence fresh-surface transfer in the absence of other carriers, by construction.

### Field-level Monte Carlo contrasts

Independent-field bootstrap intervals (3,000 bootstrap draws of the 48 simulated fields per arm) for restricted-time differences were:

| Mechanism | Immediate-transfer advantage (fresh minus 0 h) | Washout delay (37 h minus 0 h) |
|---|---:|---:|
| Passive | -1.88 [-5.85, 2.08] | 1.45 [-2.69, 5.59] |
| Immutable selection | 16.24 [13.21, 19.28] | -0.08 [-0.53, 0.35] |
| One-bit memory | 14.65 [11.68, 17.81] | 8.96 [7.20, 11.04] |
| Scalar, no daughter inheritance | 13.76 [10.90, 16.74] | 13.75 [10.88, 16.70] |
| Decaying scalar memory | 13.68 [10.87, 16.73] | 13.67 [10.79, 16.59] |
| Scalar, no decay | 13.68 [10.73, 16.72] | 0.35 [-0.04, 0.74] |

These are *within-model computational precision intervals* and do not account for parameter uncertainty, alternative microscopic rules, or between-experiment biological variability.

### What is and is not minimal

The baseline passive process can generate colonies and family trees without memory, but under fixed incoming composition and virgin surfaces it cannot systematically encode a different exposure history. At least one component retaining historical information is required by that conditional model specification: a state inside agents, persistent phenotype composition, persistent surroundings, or residual agents.

A single intracellular scalar is **sufficient** for the coarse transfer/short-hold/long-hold signatures in this model class. **Direct daughter inheritance is not necessary for those coarse signatures** because rapid post-birth reacquisition can compensate. A one-bit stochastic state also reproduces immediate transfer and long-hold washout with the chosen rules; its 14-h effect is weaker and delayed in this particular parameterization. There is no universal proof that a continuous scalar is necessary or that a one-bit state is insufficient. The specific nonlinear response and holding dynamics were chosen with known phenomenology in mind.

## 2. Local sensitivity

For a fixed model grammar, the decay constant (7, 10, 14 h), memory-response midpoint (0.075, 0.115, 0.160) and Hill exponent (4, 7, 10) were crossed in 27 configurations. Each configuration used 10 independent fields per arm, with different seeds. A configuration was counted successful when the restricted-time difference satisfied all three **illustrative qualitative** rules:

- Immediate transfer at least 4 h sooner than fresh.
- After 37 h, colonization at least 4 h later than immediate transfer.
- A 14-h hold changes the restricted time by between -1 h and +3 h relative to immediate transfer.

**23/27** passed. The four failures occurred when slow decay severely reduced the long-hold washout, or a high response threshold and fast decay weakened the 14-h effect. These are coarse demonstration thresholds, not laboratory-calibrated acceptance ranges. Ten synthetic fields per configuration are inadequate for exact success probabilities; the grid is a local robustness diagnostic.

## 3. Separating post-birth reacquisition from true inheritance

The parent-state model and the no-inheritance model have similar macroscopic restricted-time curves, despite different daughter-state rules. If a daughter could be isolated immediately at birth, transferred to a virgin surface, and measured *before acquiring a new surface signal*, the model predicts contrasting initial stable-attachment probabilities for a parent state `m=0.8`:

- With direct daughter inheritance, approximately **0.0850**.
- Without direct daughter inheritance, **0.0035**.

The approximately **24-fold difference is a synthetic, conditional prediction**, not a validated effect size. Delaying the daughter assay permits post-birth sensing, reducing its discriminating value. An existing fixed adhesive phenotype could also support enhanced daughter attachment and therefore requires prospective phenotype/lineage tracking to distinguish selection from acquired inheritance.

## 4. Microscopic oscillatory signature needs a richer dynamical class

The one-scalar colonization model does not explicitly model cAMP or type IV pilus time series. To explore the **additional** published cross-correlation signature, we simulated two independent signal classes:

- A scalar continuous-time Ornstein–Uhlenbeck state, with motility readout negatively coupled to it.
- A two-state, damped, stochastic rotational system, with the same style of readout.

In 160 independently simulated lineages, the scalar model had a cAMP-proxy/motility-proxy cross-correlation of **-0.31** near +5 h; the coupled model had **+0.13** near +5 h, while both were negative at zero lag. The intervals shown in the figure are lineage-resampled bootstrap intervals. These are free synthetic signals with an illustrative angular frequency of 0.60 rad/h and **were not fitted to the biological fluorescence traces**.

A continuous-time *scalar linear autonomous Markov OU* process with an instantaneous negative readout has a monotonically decaying negative cross-correlation. A second coupled dynamic degree of freedom permits oscillatory sign reversals. This is not a universal claim about all nonlinear scalar processes, time-delayed models, or non-Markov systems. The original Lee et al. paper already used a mechanistically richer coupled stochastic model; this example is not a newly discovered oscillator mechanism.

## 5. Interpretation and publication value

The strongest theoretical result from this exploratory implementation is an **endpoint-dependent minimum**. Coarse colony growth and transfer washout can be produced by elementary agents carrying only one persistent scalar state. In this model class, the same colony-level observables do not distinguish inherited daughter state from rapid reacquisition. Oscillatory biochemical readouts impose stronger dynamical requirements. The synthetic daughter-transfer assay supplies a potential discriminatory test of these distinct explanations.

A publishable contribution should combine this minimal-agent capability map with the existing v5/v7 **intervention-identifiability** results. The generative part shows which observables can be matched by simpler mechanisms, while the causal analysis determines which experiments separate those possibilities. Simply showing that a simple agent model recapitulates an already published bacterial memory phenotype would face a serious novelty objection, especially given the original 2018 stochastic model.

No new bacterial measurements were performed. The original v7 analysis' empirical weaknesses remain. The full set of generated ablations is not a validated model-selection exercise, and **the findings do not make the manuscript desk-rejection-proof**. The highest-value next test is an *independently implemented*, pre-specified holdout of a full empirical observable that was not used to choose parameters, preferably a lineage-resolved transfer or daughter-isolation assay, with replicate-level uncertainty and intervention-specificity calibration.

## Sources

- Lee CK et al. (2018). Multigenerational memory and adaptive adhesion in early bacterial biofilm communities. *PNAS* 115:4471–4476. DOI: https://doi.org/10.1073/pnas.1720071115.
- Gajrani S, Ye X, Ratzke C (2025). Environment-mediated interactions cause an externalized and collective memory in bacteria. *ISME Journal* 19:wraf173. DOI: https://doi.org/10.1093/ismejo/wraf173.

## Code and reproducibility

Run with Python 3.10+ and NumPy, Pandas and Matplotlib installed, from the package root:

```
python -m pip install -r requirements.txt
python src/run_capability_sweep.py --reps 48
python -c "import sys;sys.path.insert(0,'src');from run_robustness import run;run(10)"
python src/oscillatory_signatures.py
python src/run_surface_conditioning.py
python src/make_results.py
python -m unittest discover -s tests -v
```

Model assumptions and all defaults are defined in `src/minimal_agents.py`. The stochastic fields and bootstrap estimates are provided under `results/`, with figures under `figures/`. The reproducibility tests passed in the producing environment. All generated files belong to this exploratory supplement, separate from the original experimental archive.
