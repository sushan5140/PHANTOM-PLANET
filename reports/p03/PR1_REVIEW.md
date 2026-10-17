# Read-only review of PR #1

Reviewed 2026-10-17: [P01: synthetic RV baseline and optional juliet fitting](https://github.com/sushan5140/PHANTOM-PLANET/pull/1).
Head `e4d0eaca52c7b57d92c40c09f942d67f38311fe5`, base main `0547f97`.
It was open with no reviews/comments when inspected. No review comment was
posted, no PR state was changed, and P03 does not merge it.

Source review covered all P01 Python files, requirements, tests and README.

| Finding | P03 treatment |
|---|---|
| `sim.py` creates deterministic toy activity and circular sinusoidal planets | Preserve as a teaching baseline; not physical stellar validation |
| `baseline.py` ignores heteroscedastic errors and chooses the largest scanned power | Preserve; no FAP, detection or planet classification inferred |
| `juliet_fit.py::fit_one` has synthetic period/K bounds, one circular planet, no activity model, no convergence verification | Preserve unchanged; separate real-data likelihood and proper shared priors |
| `run_p01.py` runs simulations and writes descriptive candidates | Preserve unchanged; real data have a separate entrypoint |
| Four P01 tests verify reproducibility/shapes/basic range, not astrophysical identification | Retain and execute alongside independent P03 likelihood and ingestion tests |
| README explicitly says real data and official tutorial have not been run | Preserve historical README; append a clearly dated P03 section |

No defect in P01's stated exploratory scope required rewriting its source.
Its outputs cannot be promoted into scientific validation. P03 addresses that
boundary with documented data provenance, a real-data comparison, reference
reproduction and uncertainty checks. `git diff` against the P01 head confirms
that its Python sources and existing tests are unchanged.
