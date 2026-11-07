# Adversarial review and analytic response, revision v8

Date: 7 November 2026. Model and archival data provenance: Lee et al. (2018), PNAS 115:4471-4476, DOI 10.1073/pnas.1720071115; derived three matched FC1/FC2 series in squareshorts/colonization-history; v7 and minimal-agent v1 synthetic analyses. This report is an adversarial scientific audit, not a guarantee of journal acceptance.

## Revised editorial decision

Major revision before submission remains appropriate. The strongest positive contribution is a precise falsification of an observationally linear nonnegative exponential-memory explanation of the three published washout point estimates, with a constructive population-composition competitor and a broader test of synthetic-agent assumptions. Biological mechanism identification is still unachieved. The new quantitative benchmark does not turn the previously published observations into independent validation.

## Criticisms tested and outcomes

| Likely hostile criticism | Executed corrective test | Finding | Status |
|---|---|---|---|
| The model was tuned to known washout data | Replaced a single local 27-grid robustness result by 160 LHS parameter vectors and five independently seeded simulated fields per arm, comparing four agent classes | 39/160 scalar, 29/160 one-bit, 39/160 without inheritance, 0/160 immutable-selection variants met all three selected criteria. Only 96/160 scalar settings matched the early pattern | Partially addressed. The parameter ranges and criteria remain researcher-chosen and the 37-h outcome was known |
| Scalar memory may be neither minimal nor necessary | Tested one-bit state and explicit no-daughter-inheritance agent; analyzed immutable selection with changing composition during suspension | Both one-bit and no-inheritance models can generate the pattern. A two-phenotype population turnover model exactly interpolates three lags without acquired intracellular memory | Resolved as a *negative* claim against uniqueness, not demonstrated biological selection |
| The reported lags do not constrain anything specific | Algebraic exponential advantage, nonnegative multi-exponential log-convexity, and bounded ±0.25,0.5,1,2,3,4-h stress tests | Fitted at 0 and 14 h, a single exponential predicts 6.64 h at 37 h versus reported 13.41 h. Any nonnegative exponential mixture with a linear observation mapping requires normalized residual advantage >=0.8478; observed 0.3222. The incompatibility persists under ±1-h worst-case independent perturbations, not ±2-h | New conditional theoretical result; cannot treat fitted curve CIs as biological replication |
| Washout rules out selection | Constructed two immutable phenotypes with fitness-differential population turnover, logistic frequency dynamics and linear lag response | Exact interpolation at f0=0.99036, relative growth-rate disadvantage 0.14567/h, f14=0.9304, f37=0.3191. These parameters are unmeasured and potentially extreme | Contradicts exclusion of selection; biological calibration still lacking |
| Matching threshold times means matching published lags | Generated complete 40-h occupancy trajectories and compared restricted individual threshold time to ensemble log-linear count onset (16 stochastic fields per arm) | Default scalar immediate-transfer restricted threshold 5.31 synthetic hours versus log-linear intercept -0.94 synthetic hours. The latter is extrapolated and cannot equal the 2018 lag fit. | Addressed by openly disallowing numerical cross-estimator validation |
| Thousands of cells provide high power | Kept three matched *experimental units* in empirical analysis; computed hypothetical noncentral-t planning sensitivity based on n=3 observed pairwise SD | Pairwise SD 0.1326. Hypothetical n~16 independent pairs for mean difference 0.10 at 80% power if true SD equals this estimate; n~34-58 if it is 1.5-2x. No valid prospective power is claimable from n=3 heterogeneous pairs | Cannot be resolved without new independent biological experiments |
| Partial identification is guaranteed by more data | Retained v7 valid/invalid assay-sensitivity bounds and nonlinear misspecification tests | The bounded method protects against false attribution only conditional on valid measurement bounds and additive model; invalid calibration and nonlinear response defeat it | Explicit model-dependent conditional limitation |
| Intervention purity is assumed, not demonstrated | Counterfactual algebra kept separate from biological feasibility | The six-arm contrast identifies an *ordered* allocation only if the selective resets can be executed without off-target changes | Unresolved experimentally |
| Signaling oscillator was a disconnected second model | Removed oscillator graph from main paper and relegated it to supplementary mathematical control | No unified molecular signaling and colony simulation is claimed | Framing fixed; unified biology still unavailable |
| Existing random-forest classifier has privileged inputs | Kept original classifier results out of the main inference; repository README already warns about exact residual-cell features | No practical mechanism-recovery accuracy claim from oracle features | Scope fixed; legacy code remains as archived provenance |
| Novelty compared with literature | Explicit comparison to Lee 2018, Zhao 2013, Gajrani 2025 and standard identifiability theory | Novelty is an observation-model constraint and demonstration of the experimental consequences of generative equivalence; generic causal nonidentifiability is established theory | Reduced but meaningful desk-rejection risk |

## What still cannot be fixed analytically

A high-ranking biological journal can reasonably demand biological calibration of arrival and attachment rates, representative independent flow-cell replication, phenotype-frequency measurements during hold, cross-reactivity testing of any proposed cell-state / surface-trace assay, validated selective erasure, and a held-out *new* lineage-resolved prediction. None was collected here, and no number of synthetic agents can replace these observations.

## Biologically discriminating next measurements

A useful follow-up combines (i) high-adhesion phenotype frequency and planktonic growth trajectories at 0, 14, 24 and 37 hours; (ii) daughters sampled immediately and before new surface sensing; and (iii) a fresh-surface challenge with a trace-specific sentinel and viable-residual-cell check. Concentration, inoculum composition and true independent biological culture unit should be held fixed or separately randomized. The 24-hour lag is a forward prediction of the currently fitted *phenomenological* models; it is not a preregistered study outcome. Model-specific synthetic 24-hour predictions are in `results/v8/selection_composition_counterexample.json` and `published_lag_model_audit.json`.

## CAPES journal suggestion (verification status)

**Theory in Biosciences**, Springer Nature, ISSN 1431-7613 / 1611-7530, is a hybrid journal whose 2026 scope includes artificial life, bioinspired modeling, mathematical biology, robustness and population biology. The journal rejects standard analyses without new conceptual insight. Its reported open-access list APC is GBP 2490 / USD 3390 / EUR 2790, but the publisher's CAPES transformative agreement can cover eligible Springer hybrid titles for participating Brazilian corresponding authors accepted 2026-01-01 through 2028-12-31. Universidade Federal do Pará appears by name in the publisher's eligible Springer institution list. CAPES requires the corresponding author's ORCID to be linked to Meus Dados CAPES/Sucupira and CC BY. The source lists which Springer journals are included in an XLSX attachment; the exact inclusion of this title in the 2026 title list was not independently verified because the live XLSX download was inaccessible in the analysis environment. Thus **APC COVERAGE IS CONDITIONAL, NOT CONFIRMED for the specific title**. Verify the title in the live `Springer portfolio` spreadsheet and seek publisher eligibility confirmation before submission or accepting an OA fee.

Sources:
- https://link.springer.com/journal/12064/aims-and-scope
- https://link.springer.com/journal/12064/how-to-publish-with-us
- https://www.springernature.com/br/open-science/oa-agreements/brazil/capes
- https://www.gov.br/capes/pt-br/assuntos/noticias/saiba-como-publicar-artigos-em-acesso-aberto

Alternative candidate: **BioSystems** (Elsevier) has a stronger artificial-life remit; Elsevier's CAPES contract covers selected hybrid titles only, and this specific title's inclusion has not been verified. An unqualified APC promise would be inappropriate.

## Source and manuscript integrity

The reconstructed v8 Overleaf package is an updated version of the previous reconstructed v7 manuscript. It is not a line-by-line change to any unavailable earlier private manuscript. New stochastic parameter priors and all criteria are disclosed. No new biological assay was performed. The v7 archive DOI 10.5281/zenodo.23990601 is **not** asserted to archive new extensions. Executable tests and readme accompany the reproducibility archive.
