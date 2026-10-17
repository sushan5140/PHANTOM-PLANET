# Scientific validation notes

This document describes the inference boundary. Measured values and comparisons
are in MODEL_COMPARISON.md and the JSON/CSV artifacts.

## What is independently checked

- File checksums, schema and duplicate IDs; all numerical target summaries
  compared with P02; the original KOA CSV matched to the interactive JSON table.
- Full-cadence and nightly series retained; formal errors remain distinguishable
  from fitted stellar jitter. No hidden RV clipping or periodogram-based selection.
- Official juliet RV tutorial with unchanged priors and real measurements.
- Circular RV implementation numerically checked against juliet's RadVel engine.
- Normalized GP likelihood checked against an independent SciPy implementation.
- Zero GP-amplitude limit gives exactly the corresponding Keplerian likelihood.
- Kernel positive definiteness, conditional prediction and common-prior equality
  tested; normalized nested evidence checked on an analytically soluble problem.
- Independent seeds, evidence trajectories, sample weight ESS and posterior
  agreement are numerical diagnostics, distinct from astronomical validation.

Exact historical juliet benchmark replication remains unresolved: the current
two-planet run and independent reference attempts can select different period
modes and disagree in evidence beyond quoted errors. The saved likelihoods agree
with independent calculations, but that does not certify global exploration.
Accordingly this deliverable is a conditional feasibility analysis, not a fully
validated planet-detection pipeline. The target's preferred repeat diagnostics
pass, while its broad-prior C stress test retains a posterior-shape warning.

## What the comparison can say

An evidence ratio describes the integrated probability density of these 19
nightly RV means under a specified likelihood and prior. Model A is a restricted
circular-orbit model; B is one particular activity covariance; C contains both.
Preference for A or C means only that this coherent component helps within these
assumptions. It is not a probability that all stellar explanations are false.
Published findings provide context, not a label inferred by this analysis.
The K prior is nonnegative: positive posterior samples or a positive lower
equal-tail quantile do not establish a detection. The separate no-planet model
is needed for that comparison; its evidence is still conditional on the model set.

The simultaneous S-index is shown on its own and against RV and residuals.
Its absence of a given periodicity would not prove a planetary origin: different
activity indicators can respond differently, and phase shifts matter. Correlation
does not identify a starspot mechanism. P02's per-exposure correlation p-values
are not valid independent-night evidence and are not used as such.

## Limits that remain even if sampling agrees

1. Only 19 nightly observations over 314 days and irregular gaps; multiple
   period modes and aliases can remain. Hundreds of exposures do not create
   hundreds of independent rotation measurements.
2. GP rotation/evolution parameters may be prior dominated. The primary rotation
   prior comes from later TESS observations; stellar patterns evolve, and broad
   prior sensitivity is essential. Nested sampling can still miss modes.
   The observed S-index levels differ across observing seasons; a stationary
   QP covariance and common RV trend are a limited approximation to that evolution.
3. Nightly averaging does not remove all granulation/supergranulation. The free
   nightly jitter is an approximation. Correlated nightly errors, heteroscedastic
   stellar noise and uncertainty in the binning process are not fully modeled.
4. The restricted models omit an explicit short-term GP, eccentricity, full
   three-companion dynamics, photometry and joint RV/S covariance. Their evidence
   is not the evidence of the richer published model. A GP may absorb planets;
   a Keplerian may absorb activity or instrument effects.
5. In-sample GP residuals shrink through conditioning. Small residual RMS or
   high in-sample coverage cannot substitute for independent prediction. There
   is no blind held-out-star result, injection-recovery calibration, false-positive
   rate, detection completeness, mass accuracy claim or new-planet announcement.
6. Literature comparisons share some underlying HIRES observations. They are
   consistency checks, not independent replications. The archive has 189 points
   rather than the 238 reported in fuller studies; reductions and selections differ.
7. The BJD label is preserved; UTC/TDB and sub-exposure timing not independently
   verified. No timing-precision claim is made. Relative RV zero point is nuisance.
8. Results are reproducible to Monte Carlo tolerance, not necessarily bitwise
   across operating systems, BLAS implementations or future package versions.
   The Windows lock contains a platform-specific curses dependency. The generic
   requirements file uses an appropriate platform marker.

## Interpretation of null or inconclusive results

Failure to prefer a planet in this subset does not disprove the published
planets. Failure to prefer a GP does not establish an inactive star. The target
is a known host chosen for follow-up, not an activity-only negative control.
Any fresh periodogram maximum is labeled descriptive and cannot become a new
planet claim without a new independent analysis and validation.

## P04 gate

Evidence of active stars confounding RV fits is an established problem, not a
new research opportunity by itself. This project needs a specific reproducible
failure relative to established methods and independent validation before it can
claim an improvement. If this phase shows only known signals and sensitivity
to known assumptions, recommend **DROP**, retaining the reproducible work.
