# Scientific and technical audit — corrective revision

## Scientific corrections

- The reanalysed image-derived lags, all using the same operational estimator, give naive mean 15.503 h, immediate-transfer mean 5.177 h, single 14-h hold 3.71 h and single 37-h hold 12.59 h. Relative advantages 1.142 and 0.282. Because 14-h advantage is larger than the immediate mean, these data do not show monotone washout. Neither held condition has independent replication. The source study's published curve-fit point estimates are separate and used only in a model-conditional algebraic comparison.
- Exact independent image-derived source: the original Lee et al. binary stacks were segmented with the public `empirical/extract.py` (connected components >=40 pixels, 30-s sampling). `empirical/occupancy.py` performs smoothing and occupancy-crossing fits. The present revision recalculates from the stored derived occupancy/lag results; it does not repeat raw-image segmentation.
- The identical 4.97-h minima are a rounding coincidence: FC1 series 1 yields 4.97405 h at smoothing 40, fit window 40–400; FC2 series 2 yields 4.97246 h at smoothing 20, fit window 10–250. Maximal technical FC2-minus-FC1 differences for the matched pairs remain negative: -1.39, -4.79 and -4.96 h. These are sensitivity checks, not statistical confidence intervals.
- Two-arm attrition is implemented, including identical loss in both arms, loss in sham at half the pulse rate, greater loss in sham, and phenotype-specific attrition. The sharp full-assignment contrast remains valid and can lose virtually all detection power under 10% uniform loss in both arms. The recovered-only contrast retains more power under uniform missingness but is susceptible to differential selective recovery. Monotone-selection trimming is excluded from the new scenario suite because cross-world recovery monotonicity was not established.
- Nominal 40% dropout among low-adhesion cells caused approximately 29% population attrition. Recovered-only false joint declarations were 4.42% and 6.08% under the two kinetic laws, respectively, in 1,200-study null batches; bounded joint false declarations were zero there.
- Repeated-contact discrimination at p_high=.10, p_low=.03 requires 22,633,507, 459,080 and 14,901 cells at the 0-, 14- and 37-hour holding intervals, respectively, for 80% known-parameter power. Monte Carlo detection at those sample sizes was approximately 79–81%. These remain optimistic because the phenotype probabilities are assumed known and the acquisition process is assumed homogeneous.

## Mechanical fixes

- Rewritten and tightened Introduction, Methods, Results and Discussion according to the supplied scientific-prose card.
- Abstract under 200 words; numbered bibliography; correct TeX accent commands for Delbrück, Åström and Klingmüller; section-specific citations.
- Source microscopy derivation and independent estimation procedure clearly documented.
- All figures in main text and supplement are referenced; `localization_selectivity.pdf` appears as a cited supplementary figure. Figure legends were revised for consistent solid bound versus dashed recovered-only semantics. The previous excessive agent-model supplement section was removed.
- The article has CRediT-formatted author roles, a plain funding statement, a competing-interest statement, and a detailed declaration naming OpenAI ChatGPT (GPT-6). The author must check these personal declarations and any external AI-review tool whose identity was not supplied.
- No internal version numbers, audit process, Overleaf workflow, review language, or similar development bookkeeping appear in the article or supplement scientific prose.

## Remaining author actions / limits

1. The author must review and independently validate the research in light of the Royal Society's AI policy, which restricts AI from replacing researchers in scientific inference and data analysis. Disclosure alone may not make a submission eligible. The published policy is https://royalsociety.org/journals/ethics-policies/artificial-intelligence/.
2. Confirm CRediT roles, funding, any CAPES/CNPq support or fellowships, and competing interests before the final submission. No grant or fellowship information specific to this project was supplied.
3. Confirm whether any additional identifiable generative AI tool (used to produce the external audit) requires disclosure. Do not invent a tool name.
4. No new Zenodo DOI has been minted at the user's explicit instruction. The new code is in GitHub and the underlying image deposit has a DOI, but a version-specific code archive DOI is outstanding if the publisher requires one.
5. The source is a revision of the reconstructed earlier LaTeX package. A line-by-line comparison against the author's unavailable original private Overleaf project could not be performed. Original raw binary microscopy images have not been independently reprocessed in this revision. Single-run holding conditions prevent biological confidence intervals and direct validation of the washout kinetic law.

The revised manuscript is a coherent theoretical modeling study with experimentally unverified carrier-specific tests. Bibliographic and technical integrity do not resolve those biological and editorial-policy limitations.
