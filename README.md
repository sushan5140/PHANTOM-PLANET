# PHANTOM-PLANET — P03 real-data feasibility analysis

Research only. No new planet or new detection method is claimed.

P03 imports the public Keck/HIRES K2 RV/activity release and the supplied P02
audit, then compares a circular Keplerian, a rotation Gaussian process, and their
combination for **HD 73344 / EPIC 212178066**. The correct observing-night
boundary gives **189 exposures on 19 nights**, across 314.19 days. A fourth
offset/trend/jitter reference is included. Literature-conditioned and broader
priors are explicitly distinguished.

- [Model comparison and P04 recommendation](reports/p03/MODEL_COMPARISON.md)
- [Methods, priors and assumptions](reports/p03/METHODS.md)
- [Scientific validation and limitations](reports/p03/SCIENTIFIC_VALIDATION.md)
- [PR #1 source review](reports/p03/PR1_REVIEW.md)
- [Attempts and corrections](reports/p03/WORK_LOG.md)
- [Diagnostic figures](reports/p03/figures/)

## Reproduce P03

Tested on Windows with CPython 3.10.17. Install into a local environment:

```powershell
uv venv --python 3.10 .venv
uv pip install --python .venv/Scripts/python.exe -r environment-p03.lock.txt
.venv/Scripts/python.exe -m pytest -q --basetemp .p03-cache/pytest
.venv/Scripts/ruff.exe check --isolated --select E,F --ignore E501 src/phantom/p03_data.py src/phantom/p03_models.py scripts tests/p03
```

On other platforms use `requirements-p03.txt` (platform marker for curses),
`python` from the activated environment, and the platform's environment path.
Compiled dependencies may need a C/C++ compiler when no wheel is available.
The lock records the complete environment used here; it is Windows-specific.

From a clean checkout, run the complete sequence with
`./scripts/reproduce_p03.ps1`, or issue its Python commands individually.
It downloads missing public numeric inputs only after a checksum check; the
checked-in data allow offline numerical reproduction after dependency installation.
The user-supplied P02 contents are retained, so the original zip is not required.

The sequence runs an analytic evidence check, both official TOI-141 `juliet` RV
examples, independent likelihood checks, two seeds for A/B/C and the simpler N
reference, broad-period and amplitude-prior sensitivities, and a 1000-live-point
combined-model check. It then generates figures, CSVs and posterior diagnostics.
Individual fits take minutes; complete serial reproduction can take several hours.
Independent target runs can use `--reuse-audit` concurrently after the audit.
Results are reproducible to Monte Carlo tolerance;
fixed seeds do not promise bitwise equality across platforms or library versions.

Target sampling requires a completed official example and checked likelihood.
The documented tutorial's absolute evidence is only approximately reproduced;
the discrepancy is retained, not silently declared numerical agreement.
Full details are in the tutorial JSONs and work log.

Existing target result files are not overwritten: use a fresh `results/` directory
or a new seed/live-point count. `--reuse-audit` reads a checksum-verified processed
CSV without rewriting it, for independent concurrent runs. Original measurements,
configuration, exact versions, posterior samples and diagnostics are retained.
`results/` is an ignored working directory; portable deliverable samples and
summaries are in `reports/p03/`.

## P01 historical snapshot — minimal feasibility experiment

The section below describes the earlier P01 state; its statements about data
not yet being downloaded or methods not implemented apply to that snapshot.
P01 Python sources and tests are preserved unchanged from PR #1.

Research only. **Not a planet detector.** No model trained; no astronomical discovery claimed.

## Start

```bash
python -m pip install -r requirements.txt
python run_p01.py --seeds 10
python -m unittest discover -s tests -v
```

Produces 40 synthetic radial-velocity series (four conditions, 10 seeds) and strongest *candidate* periods. The starspot signal is an intentionally simplified quasi-periodic toy model, not SOAP or verified HARPS data. A periodogram peak is **not** a detection.

## juliet: optional sanity check

```bash
python -m pip install -r requirements-juliet.txt
python -c "import sys; sys.path.insert(0,'src'); from phantom.juliet_fit import fit_one; print(fit_one('results/synthetic/planet_000.csv','results/juliet_planet'))"
```

This fits **one circular Keplerian model**, not a validated activity-vs-planet comparison. Do not interpret the evidence as a planet probability. Follow juliet's official RV and GP tutorials before building model comparisons with carefully matched priors.

## Next gate (NOT IMPLEMENTED)

1. Run official juliet tutorial and verify reproducibility.
2. Add planet-only, quasi-periodic GP activity-only, and planet+activity models with documented priors and sampler convergence checks.
3. Import a small set of real RV observations and activity indicators, with licensing/target audit.
4. Evaluate blinded hold-out stars and more realistic simulations.
5. Identify a systematic failure relative to existing astrophysics methods or **DROP**.

## References

- https://juliet.readthedocs.io/en/latest/tutorials/rvfits.html
- https://juliet.readthedocs.io/en/latest/user/priorsnparameters.html
- https://www.eso.org/qi/catalog/show/416
- https://koa.ipac.caltech.edu/Datasets/K2ExoDiversity/K2_RV_Activity.php

No real observational data has been downloaded.
