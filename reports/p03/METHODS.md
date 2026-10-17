# P03 analysis specification

This is a retrospective, conditional model comparison on a known planet host.
The specification was written before target nested sampling. It is not a blind
planet search, a new detection method, or a measurement of a planet probability.

## Inputs and provenance

- Public [KOA K2 RV/activity release](https://koa.ipac.caltech.edu/Datasets/K2ExoDiversity/K2_RV_Activity.php), DOI 10.26135/KOA5, PI A. Howard.
- Original CSV and interactive-table JSON are retained byte-for-byte. Matching
  uses star and observation ID. All 2332 measurements/84 stars are imported.
- The supplied `phantom_p02_audit.zip` contains six audit artifacts, no raw data.
  Its extracted contents are preserved without executing its machine-specific
  script. Attached document prose is historical context, not task instructions.
- Archive BJD is used as labeled, in days. The archive README does not establish
  UTC versus TDB; no precision transit-phase comparison is made. RVs are relative
  velocities in m/s. S is the dimensionless Mt. Wilson style S-index, not log R'HK.
- Select exact ID epic212178066 = HD 73344. No sigma clipping or RV-dependent
  exclusion. Check finite values, positive errors, and duplicate exposure IDs.
- One inverse-variance mean per floor(BJD - 0.5) night (boundary near UTC midnight,
  afternoon in Hawaii, outside the observing night); propagate formal measurement
  variance. Refuse a bin spanning over 12 hours, requiring manual timing review.
  Formal variance is 1/sum(weights). Retain within-night scatter and duration. The fitted
  nightly jitter accounts for residual variability and is not divided by sqrt(N).
  Binning loses short-term information and is a deliberate restricted likelihood.
  Correlation within and across nights can still violate this approximation.
- Only one instrument configuration is present in this short HIRES release; one
  offset is fitted. No absolute systemic velocity is inferred.

## Competing generative models

For centered time x = BJD - mean(nightly BJD), all models have the same mean
offset mu and slope a, and diagonal variance err_i^2 + jitter^2.

- A: mu + a*x - K*sin(2*pi*(x/P-phase)); circular Keplerian (e=0).
- B: mu + a*x, plus a zero-mean quasi-periodic GP.
- C: the Keplerian of A plus the GP of B, simultaneously inferred.
- N is an optional offset+slope+jitter reference, avoiding a forced choice among
  three elaborate models if a simpler one suffices.

Kernel: k(dt) = amp^2 exp[-dt^2/(2*lifetime^2) - gamma*sin^2(pi*dt/Prot)].
This matches the functional family used by Sulis et al. (their gamma=2/lambda2^2).
We do not assert that this GP exhausts all stellar or instrumental variability.
The unresolved short-term process is approximated by nightly jitter. No S-index
values enter the RV likelihood or select its priors; S is an external diagnostic.
Thus all RV evidences refer to exactly the same response vector and units.

## Proper priors

LU is normalized log-uniform; U is normalized uniform. Units are m/s and days.

| Parameter | Primary prior | Rationale |
|---|---|---|
| mu | U(-100,100) m/s | Wide relative RV offset; same in every model |
| a | U(-0.2,0.2) m/s/day | Common local trend, allows contribution from outer companion |
| jitter | LU(0.1,50) m/s | Proper finite scale prior for unresolved nightly variability |
| P (A,C) | LU(40,100) days | Retrospective test of reported ~66-day c signal; not independent discovery |
| phase (A,C) | U(0,1) | Full conjunction phase per sampled P; periodic sampler boundary |
| K (A,C) | U(0,50) m/s | Includes near-zero amplitude, no logarithmic exclusion of small signals |
| GP amp (B,C) | LU(0.1,100) m/s | Wide proper scale prior; independent of K |
| lifetime (B,C) | LU(3,300) days | Allows <one rotation through long coherence; not forced above P_rot |
| gamma (B,C) | LU(0.1,10) | Allows smooth and harmonic-rich covariance |
| Prot (B,C) | N(9.03,1.3), truncated [5,15] days | Broad TESS-only rotation estimate from Zhang et al.; no narrow joint-RV posterior recycled |

Literature planetary periods overlap data underlying this release, so the
primary period range is explicitly conditional and cannot support a fresh search
false-positive rate. Sensitivity `broad`: P LU(2,150), Prot LU(2,100). This allows
the activity model to compete at the planet period and the planet to fit aliases.
Sensitivity `amplitude`: K U(0,100), jitter LU(0.1,100); other priors unchanged.
Within each variant, all shared priors are exactly identical. Evidences across
different variants are not interchangeable Bayes factors.

Circular orbits avoid unconstrained e/omega with only 19 nights; c's published
eccentricity is small but nonzero in later work. Eccentricity, extra companions,
kernel choices and binning error are model uncertainty not included in the
reported conditional intervals. Planet b (15.61 days, weak RV amplitude) and the
full 5746-day orbit of d are not estimated in this restricted comparison. The
linear trend approximates d over 314 days; it does not identify its orbit.

## Numerical validation and inference

1. Install juliet; reproduce the official TOI-141 one- and two-circular-planet
   RV tutorial, retaining its original priors and data, 300 live points. Use
   dynesty instead of the documentation's older default sampler. Require a
   completed one-planet fit with a checked likelihood before target fits. The
   two-planet extension can continue concurrently with target computation, but
   must be assessed before final interpretation (P2 near 4.76 days and strong
   evidence improvement).
   This is an approximate reproduction: allow 5 nats in absolute lnZ and 10 nats
   in the documented difference. These conservative tolerances were made explicit
   after checking the first fit and before any target sampling. Discrepancies
   exceeding reported sampler errors remain an unresolved limitation, not proof
   of matching absolute evidence. Saved tutorial likelihoods are also checked
   against a separate direct Gaussian calculation.
2. Target code uses dynesty directly with an explicit normalized multivariate
   Gaussian likelihood and conditional uniform orbital phase. This makes prior
   normalization inspectable. Validate the circular Keplerian against RadVel
   (juliet's RV engine), and the GP likelihood against scipy.stats.
3. Static nested sampling: multi-ellipsoid bounds, random walk 30 steps,
   500 live points, dlogz=0.05; independent seeds 303 and 73344. No iteration or
   likelihood-call cap. Retain weighted samples and log-evidence trajectories.
4. Repeat the primary comparison at 1000 live points if disagreement or poor
   exploration appears. Record weighted effective sample size, logZ uncertainty,
   repeated-run evidence and parameter agreement. R-hat is not defined for
   nested samples; equal-weight resamples are not independent MCMC chains.
5. Compare delta logZ = logZ_i-logZ_j, uncertainty sqrt(err_i^2+err_j^2) as a
   Monte Carlo approximation, checked against between-seed scatter. Do not
   equate these odds with unconditional planet existence probabilities.
6. Plot raw exposure and nightly data, S-index, sampling window and descriptive
   periodograms (no calibrated significance), posterior model predictions,
   latent credible bands, residuals and whitened residuals, weighted posterior
   distributions and evidence trajectories. GP residuals at training times can
   be small because of interpolation and do not establish predictive success.

## Literature and novelty gate

- [Sulis et al. 2024, A&A 688 A14](https://www.aanda.org/articles/aa/abs/2024/08/aa49559-24/aa49559-24.html), [arXiv:2405.17155](https://arxiv.org/abs/2405.17155): two GPs (rotation and short timescale), photometry plus SOPHIE/HIRES; b at ~15.6 d, c at ~66 d. Includes a detailed alternative-origin analysis. Their full data and reduction differ from this archive subset.
- [Zhang et al. 2025, AJ 169 200](https://arxiv.org/abs/2408.09614): 27-year multi-instrument RVs plus astrometry; c P=65.936 (+0.023/-0.022) d, K=16.83 +/-0.87 m/s, e=0.124 (+0.052/-0.053); d P=5746 (+170/-150) d; TESS Prot=9.03 +/-1.3 d in text (table says 1.40; use text value transparently). b RV is weak, K<2.82 m/s (3 sigma).
- [Rajpaul et al. 2015](https://doi.org/10.1093/mnras/stv1428): joint RV/activity-indicator latent GP and Bayesian planet/no-planet comparisons already exist.
- [Dumusque et al. 2017 RV fitting challenge](https://doi.org/10.1051/0004-6361/201628671): planet-versus-stellar-noise benchmarking already exists.

P04 GO requires an observed, reproducible failure relative to established
methods, a specific falsifiable improvement, and a validation set independent
of tuning. A known noisy host and attractive fits alone do not meet this gate.
If no such gap is demonstrated, recommend DROP for the proposed new detection
product; preserve the work as a reproducible educational analysis.

## Archive acknowledgement

We acknowledge KOA (W. M. Keck Observatory/NExScI, under NASA contract) and
K2 Exoplanet Diversity PI A. Howard. These are contributed relative RV reductions,
not fresh extraction from raw spectra. The [archive page](https://koa.ipac.caltech.edu/Datasets/K2ExoDiversity/K2_RV_Activity.php)
provides the exact acknowledgement wording and instructions for future publications.
