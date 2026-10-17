# P03 assumptions, attempts and decisions — 2026-10-17

The user's P03 request controls scope. P02 prose is treated as historical evidence,
not as instructions. No LLM training, Laya integration, website, merge or deployment.

1. Started in an empty projectless workspace. Cavemem compact search returned no
   relevant observations. Primary graph inventory had no repo checkout; cloned
   requested repository, indexed the P01 branch, checked coverage, read exact
   P01 source, tests, requirements and README. Graph freshness metadata for P01
   source was uncertain; direct source was used as authoritative.
2. Initial sandbox GitHub CLI network calls failed (`error connecting to
   api.github.com`). Repeated the requested public repo clone and PR reads with
   network permission; successful. PR #1 is OPEN, no reviews or comments, head
   e4d0eaca52c7b57d92c40c09f942d67f38311fe5, base main
   0547f97. New branch research/p03-hd73344-model-comparison starts at that P01
   head. P01 source is preserved; no merge of PR #1 is performed.
3. Reviewed PR #1: four toy cases and unweighted strongest-period baseline are
   explicitly exploratory. Optional `fit_one` is a circular single-planet fit
   with narrow synthetic K/period bounds, no activity competitor, and no prior
   validated execution. It must not be reused as a real-data detection claim.
   No review was posted externally. P03 is additive and uses a separate pipeline.
4. Expanded six P02 audit files without executing `analyze.py`; it hardcodes
   `/mnt/data` and uses unweighted periodograms and naive correlation p-values.
   Downloaded actual public original KOA CSV and JSON table plus README.
5. Agent Reach GitHub routing used gh CLI. Doctor ran; social accounts were
   irrelevant. Exa search ran for primary literature; official web search and
   direct official downloads supplied cross-checks. Jina returned a misplaced
   SVG page for the Sulis HTML URL; not used as scientific evidence. Direct A&A
   access through the web tool returned 403, so arXiv primary manuscript was read.
   No source contents were treated as executable task instructions.
6. `juliet==2.2.10`, RadVel and UltraNest compiled successfully on Windows.
   Initial resolver picked dynesty 3.1.0; pinned 2.1.5 before fits to keep the
   sampler version explicit and compatible. The initial tutorial import failed
   with `ModuleNotFoundError: No module named '_curses'`; installed
   `windows-curses==2.4.2`. Exact traceback retained in tutorial_attempt1_curses.log.
7. First tests: 9 passed, one RadVel import failure from `_curses`, one pytest
   temporary-directory access error. Fixed dependency; moved pytest temporary
   directory into workspace. Next run: 11 passed, including the RadVel numerical
   circular-orbit and independent scipy GP likelihood cross-checks.
8. Initial global ingestion rejected a nonpositive S-index error. This exposed
   an unrelated archive star's bad activity uncertainty. Adjusted import to
   retain the full archive, record invalid counts, and enforce positive errors
   on the selected target. One nonpositive S error globally, none in target.
   No exposure discarded or uncertainty imputed for HD 73344.
9. Environment freeze in the sandbox initially failed due to uv interpreter
   cache access. Repeated with a workspace-local cache and recorded exact versions.
10. Input audit: 2332 rows, 84 stars, 189 target exposures, 21 floor(BJD) nights,
    314.1899984199554-day exposure span. Five target summaries match P02 to
    numerical tolerance; every RV/activity number in CSV matches JSON by ID.
    The public release has fewer exposures than the papers' 238 HIRES points;
    no claim is made that it is the identical reduced full literature dataset.
11. Scientific specification is METHODS.md. Circular c comparison and nightly
    QP GP are restricted baseline models. Broad priors and repeated runs probe
    prior/sampling dependence. Full joint indicators, short-term GP, eccentric
    c and all three companions are not claimed to have been fitted.
12. Juliet prints optional-package warnings for absent zeus/jax/kelp, which are
    unused in these RV fits, and a Windows `cp` shell warning when copying the
    input RV file into its output folder. The original file is retained in
    data/raw; actual fitting data are read directly from that file.
13. Optimized the rotation-prior transform with scipy.special ndtr/ndtri, keeping
    the same normalized truncated normal. Added an independent scipy.stats
    truncnorm comparison at five quantiles; final unit suite has 12 tests.
    The known-answer sampler result is lnZ=-3.013119592233819 +/-0.06704565384945606,
    compared with the analytic -2.995732273553991 (passed).
14. Read official tutorial figures in the browser because their transparent
    backgrounds hid black labels in the local image viewer. The documented
    one/two-planet lnZ values are -805.302 +/-0.027 and -690.892 +/-0.023.
    Our first fit gives -803.8608078471582 +/-0.35809480620758816: a +1.4422-nat
    discrepancy, larger than the quoted sampler errors. Cause is unresolved;
    sampler/version differences are possibilities, not an established explanation.
    The tutorial gate is an approximate reproduction (5-nat absolute and 10-nat
    delta tolerance, plus the second period). It is NOT a validation of absolute
    juliet evidence to the documentation's small error bars. The independent
    known-answer sampler test separately checks our evidence normalization.
15. Checked the TOI-141 input against the documentation build commit 763f49c4:
    byte-identical SHA-256 a0e0e9d0d218ba971fed61bcabe01c3c368180c2528ed8c2824b6f72b5291506.
    Changed the reproducible fetch URL to that commit. Data drift does not explain
    the observed evidence discrepancy.
16. Compared 50 saved juliet likelihood evaluations with direct circular sine
    plus normalized Gaussian calculations. An initial absolute-only tolerance
    1e-5 failed on extreme low-likelihood prior samples (max absolute difference
    6.148044485598803e-05). Added relative tolerance 1e-10 alongside atol=1e-5;
    maximum relative difference 1.1133794395921466e-10, all points pass the combined
    tolerance. Initial failure retained in tutorial_likelihood_attempt1.log.
17. After the completed first documented RV example and independent likelihood
    check, allowed target computation while the second example continues. The
    user requested one documented example before target analysis; requiring both
    fits before any computation was an extra sequencing restriction. Both
    examples will be assessed before final interpretation. No target sampling
    preceded the first example's completion and independent check.
18. Caught a scientific binning error before any target fit completed: floor(BJD)
    turns over at noon UTC (~02:00 Hawaii) and split two observing nights. The
    21-bin counts and preliminary description above are superseded. Correct rule
    floor(BJD - 0.5) turns over near UTC midnight (~14:00 Hawaii) and yields 19
    nights, consistent with the literature's night count (not proof of identical
    reductions). Stopped the four target processes by exact inspected process IDs,
    preserving the invalid input and logs in invalid_binning_runs/. No target
    evidence was completed or used from those attempts. The tutorial continued
    unaffected. Added a regression test across the Julian-noon boundary. Final
    target fits use the corrected 19-night data and fresh hashes exclusively.
    The old bins also merged observations from adjacent physical nights, with
    a maximum duration of 0.846 days; corrected maximum is 0.22387352026998997 days
    (5.37 hours). Added a guard rejecting bins over 12 hours for manual review.
    The first two primary seeds now run in separate processes; all workers read
    an audited checksum-verified immutable CSV and write distinct run files.
19. Added portable handling for juliet's Windows `cp` warning (copy the identical
    input file with shutil before loading), workspace-local Matplotlib cache,
    and tutorial cache settings checks. These affect file handling/provenance,
    not priors or likelihoods. Git attributes preserve third-party data bytes.
    Routine ad hoc shell formatting/glob errors were corrected without altering
    measurements or fitted results. Scientific failures and cancelled runs are
    explicitly retained above and in their logs.
20. Placed P03's pytest tests in tests/p03/ (no package __init__) so P01's
    original unittest discovery still runs only its four original tests and
    does not require P03 dependencies. Full pytest discovery runs all 13 tests.
21. Primary A/B seed comparisons pass the specified repeat diagnostics. Model C's
    log-evidence agrees across two 500-live runs and the 1000-live check, but
    nuisance posterior shape differences exceed the 0.3 pooled-SD Wasserstein
    diagnostic (up to 0.463 for offset). Retained these flags and launched a
    second independent 1000-live C run, seed 2303. Evidence agreement alone is
    not declared full posterior convergence; no threshold is relaxed after seeing
    the result. Broader-period sensitivity and tutorial extension continue.
22. Juliet's two-planet tutorial completed: lnZ=-691.1300383358323 +/-0.4515772381285053,
    P2=4.785173968726813 days (16/84%: 4.784423751069262, 4.785945611151349).
    Delta lnZ=112.73076951132588, versus the documented ~114.4. Browser inspection
    of the original corner plot confirms a historical posterior near 4.760 days,
    so this is an actual numeric mismatch, not just a prose typo. Our result is
    close to the published 4.78503 days mentioned by that tutorial, but this does
    not resolve the historical mismatch. Explicitly do not claim exact historical
    fit replication. Launched a direct normalized Gaussian/sine implementation
    with the same original priors, 500 live points and independent seeds 9142/9143,
    to cross-check the current juliet result. The additional check is distinct
    from the physical target conclusions, and all discrepancies remain recorded.
23. A two-planet direct likelihood check initially failed on one sampled point
    because juliet rejects P2<=P1 with -1e101. Verified the exact source of the
    increasing-period constraint and included it in the independent check; both
    50-point checks then pass. The first direct sampler attempt omitted that
    small prior region restriction; retained it under tutorial_direct_attempt1_unordered/
    and reran the two-planet reference with the exact ordering constraint.
    The initial direct run found the historical ~4.760-day mode while current
    juliet found ~4.785 days; full historical posterior replication is not claimed.
    The benchmark metadata distinguishes an approximate sanity-check gate from
    the explicit false exact_historical_fit_reproduced flag.
24. The ordered direct two-planet reference also found the historical ~4.760-day
    mode, but lnZ=-695.1872769911744 +/-0.33837129300439917 differs greatly from
    current juliet. This demonstrates unresolved mode/evidence sensitivity in the
    tutorial extension, rather than a trustworthy exact replication. Retained
    both results. As a final global-exploration check, launched a direct 2000-live
    uniform-within-multi-ellipsoid run (seed 19143) with identical priors and
    ordering constraint. This does not change the target results or their priors.

Subsequent completed run outputs, diagnostics and final interpretation are
recorded in MODEL_COMPARISON.md and the machine-readable JSON/CSV artifacts.
