# Pre-specified analysis plan, empirical extension (v2)

Written 2026-11-05 before any cAMP-linked or full-record residence result was computed.
Results that follow this plan are reported whatever their direction; any deviation is
listed in the manuscript as a deviation.

## Data
- Wild-type YFP-CFP datasets in the transfer groups FC1, FC2, FC3, Tr14, Tr37 (10 fields).
- Bright-field binary frames every 3 s; fluorescence (YFP cAMP reporter, constitutive
  second channel stored as `dapi`, here called CFP) every 15 min. Fluorescence frame k is
  registered to binary frame 300(k+1)-1 by FFT cross-correlation, per dataset; a dataset
  whose registration correlation is below 0.3 is excluded from cAMP analyses.
- Tracking at full frame rate over the whole record (extension of v1, which tracked only
  1-h windows).

## Observables
- Attached arrival: visible in >= 2 consecutive frames (unchanged from v1).
- cAMP of a tracked cell: background-subtracted mean YFP over background-subtracted mean
  CFP inside the cell's mask at the first fluorescence frame at which the cell is present.
  Background = median of the frame outside a 5-px dilation of all masks.
- Exposure analysis window: arrivals during the first 3 h after inoculation (as v1), and,
  as the pre-specified extension, all arrivals before the population exceeds 100 cell
  equivalents (crowding changes detachment and tracking accuracy beyond that).

## Primary tests
P1 (replication with more data). Offset family versus head-start family (v1 definitions,
    lognormal heterogeneity, acquisition) on the extended arrival set. Report logL, AIC,
    LRTs against null and nested families, 200 bootstrap resamples.
P2 (state predicts retention). Among attached arrivals with a cAMP reading, add
    beta * z(log cAMP) to the individual log-hazard (left truncation at the age of the
    reading). Test beta = 0 by LRT. Fitted jointly across groups with group offsets.
P3 (mediation). Compare the FC2 offset mu_FC2 with and without the cAMP term on the same
    cAMP-read subset. Report the proportion of the offset explained, with bootstrap CI.

## Secondary
S1. Distribution of cAMP at first reading by group (Kruskal-Wallis; FC2 vs FC1 Mann-Whitney).
S2. Daughter retention: for splits on the surface, correlation between mother residence
    before division and daughter residence (Spearman), by group.
S3. Decay of the offset with holding time (Tr14, Tr37) on the extended data.

## Interpretation rules fixed in advance
- P2 significant and P3 attenuation >= 50%: exposure effect carried by a measured internal
  state (cAMP); supports cell history set before arrival over hidden adhesive variants.
- P2 significant, P3 attenuation < 50%: cAMP predicts retention but does not explain the
  exposure effect; a second carrier is required.
- P2 not significant: cAMP at the reported resolution does not predict residence; the
  selection-versus-pre-set-state question stays open, and the paper says so.

## Results and deviations log (appended 2026-11-05 after running the plan)

Implementation fixes made before any outcome was inspected: registration frames chosen
by cell content (mid-record FC1 frames are empty); file indexing by frame ranges (one
dataset stores 10 frames per file). All 10 datasets registered (correlation 0.40-0.71).

P1, first 3 h (continuous tracking, 2642 attached arrivals): offset vs null chi2_4 = 10.5,
p = 0.034; head start vs null chi2_4 = 4.2, p = 0.38; offset - head start = +3.1 logL.
Replicates v1.

P1, extended window (22888 attached arrivals): offset vs null chi2_4 = 142.6; head start
vs null chi2_4 = 21.8; offset - head start = +60.4 logL. The FC2 offset changes sign
(-0.28). Exploratory follow-up: residence drifts with time since inoculation in every
flow cell (FC2: S(5 min) 0.33 at 0-1 h, 0.12 at 5-8 h; FC1: 0.14 at 0-3 h, 0.05-0.07 at
3-20 h, 0.18 at 20-30 h). 93% of extended FC1 arrivals fall at 20-30 h on a populated
surface, so the pre-specified extension compares groups at different times and surface
states; it is reported but not interpreted as an exposure effect.

P2/P3: beta = -0.04 (LRT p = 0.36), mediated fraction -0.02. Validity check failed: the
constitutive channel is 2-10 counts above background (noise floor), and YFP reads near
zero for a large fraction of tracked cells. The per-cell ratio is not a valid cAMP
measurement in this deposit; P2/P3 are uninformative, not negative.

S1: not interpretable (channel calibration differs between datasets: 'cfp' vs 'dapi').
S2: mother pre-split residence vs daughter residence rho = -0.05 to -0.18; splits include
segmentation separations, so not interpreted.

Analyses added after the plan (exploratory, motivated by the drift above):
- Shared piecewise arrival-time effect on log hazard (0-1, 1-3, 3-8, 8-20, >20 h), all
  arrivals: offset - head start = +42.0 logL (bootstrap 24.5-63.5; offset better in 200/200);
  FC2 offset 0.58 (0.34-0.88). Held-cell offsets 1.02 and 1.29, unstable across analyses.
- Log-occupancy covariate: coefficient 0.02; offset - head start = +60.5 (42.0-81.4).
- Prospective power from the fitted residence law (design_power.py).
