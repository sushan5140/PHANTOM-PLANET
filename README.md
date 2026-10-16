# PHANTOM-PLANET P01 — minimal feasibility experiment

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
