# Carrier discrimination critical test v9

This independently reproducible addition to the `colonization-history` project addresses a specific question. Can a **predeclared intervention** distinguish two mechanism classes whose colonization-lag histories are exactly observationally equivalent?

The analytic counterexample is exact on the full lag curve. The prospective within-clone pulse challenge produces strong model-conditional separation under complete phenotype-neutral recovery in an independent single-cell attachment/detachment simulation; biased recovery or inherited cell-associated extracellular polymers defeat the intracellular attribution.

The documented numerical results and all interpretation boundaries are in `REPORT_V9.md`. The script is `analysis/discriminating_test_v9.py`. Run from this folder:

    python analysis/discriminating_test_v9.py --datasets 2000
    python analysis/finalize_v9.py
    python analysis/global_equivalence_stress.py
    python -m unittest discover -s tests -v

The study contains no new biological measurements, no measured bacterial effect sizes, and no claim of uniquely identifying bacterial memory from published lag curves. This is a versioned extension, not a replacement for the original v1.1.0 Zenodo release.
