"""Audit inputs and run nested sampling. Target inference requires tutorial gate."""
import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import dynesty
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from phantom.p03_data import prepare  # noqa: E402
from phantom.p03_models import RVModel  # noqa: E402


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit():
    import pandas as pd

    data, target, nightly = prepare(ROOT)
    previous = pd.read_csv(ROOT / "data/p02_audit/phantom_p02/all_stars_quality.csv")
    old = previous.loc[previous.star == "epic212178066"].iloc[0]
    current = {"n": len(target), "span_days": float(target.bjd.max()-target.bjd.min()),
               "rv_std_ms": float(target.mnvel.std(ddof=1)),
               "median_rv_error_ms": float(target.errvel.median()),
               "rv_s_pearson_r": float(target.mnvel.corr(target.svalue))}
    checks = {key: bool(np.isclose(value, old[key], rtol=1e-9, atol=1e-9))
              for key, value in current.items()}
    if not all(checks.values()):
        raise ValueError(f"P02 audit disagreement: {checks}")
    # Match archive table and original CSV by observation ID, not row order.
    table = pd.read_json(ROOT / "data/raw/koa_k2_rv_activity.json")
    joined = data.merge(table, on=["star", "observation_id"], validate="one_to_one")
    agrees = {key: bool(np.allclose(joined[key+"_x"], joined[key+"_y"], atol=1e-8, rtol=0))
              for key in ["bjd", "mnvel", "errvel", "svalue", "svalue_err"]}
    maximum_difference = {key: float(np.max(np.abs(joined[key+"_x"]-joined[key+"_y"])))
                          for key in agrees}
    if len(joined) != len(data) or not all(agrees.values()):
        raise ValueError("KOA table and CSV disagree")
    summary = {"archive_rows": len(data), "archive_stars": data.star.nunique(),
               "archive_nonpositive_rv_error": int((data.errvel <= 0).sum()),
               "archive_nonpositive_s_error": int((data.svalue_err <= 0).sum()),
               "target": "HD 73344 / EPIC 212178066", "target_summary": current,
               "p02_checks": checks, "archive_table_checks": agrees,
               "archive_table_max_abs_difference": maximum_difference,
               "nightly_rows": len(nightly), "bjd_min": float(target.bjd.min()),
               "night_bin_rule": "floor(BJD - 0.5), near UTC midnight / Hawaii afternoon",
               "bjd_max": float(target.bjd.max()), "excluded_exposures": 0,
               "max_nightly_duration_days": float(nightly.duration_days.max()),
               "nightly_rv_std_ms": float(nightly.mnvel.std(ddof=1)),
               "nightly_formal_error_median_ms": float(nightly.errvel.median()),
               "hashes": {str(p.relative_to(ROOT)): digest(p)
                          for p in sorted((ROOT / "data").rglob("*")) if p.is_file()}}
    (ROOT / "reports/p03/data_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    return nightly


def fit(model, variant, seed, nlive, dlogz, nightly):
    out = ROOT / "results/p03" / f"{variant}_{model}_s{seed}_n{nlive}"
    if out.with_suffix(".json").exists():
        raise FileExistsError(f"Run exists; use a different seed/live-point count: {out}")
    center = float(nightly.bjd.mean())
    problem = RVModel(nightly.bjd.to_numpy()-center, nightly.mnvel.to_numpy(),
                      nightly.errvel.to_numpy(), model, variant)
    started = time.time()
    sampler = dynesty.NestedSampler(
        problem.loglike, problem.prior, len(problem.names), nlive=nlive,
        bound="multi", sample="rwalk", walks=30,
        rstate=np.random.default_rng(seed), periodic=([4] if model in ("A", "C") else None),
    )
    sampler.run_nested(dlogz=dlogz, print_progress=False)
    result = sampler.results
    weights = np.exp(result.logwt-result.logz[-1])
    weights /= weights.sum()
    equal = dynesty.utils.resample_equal(result.samples, weights,
                                         rstate=np.random.default_rng(seed+100000))
    summary = {"model": model, "variant": variant, "seed": seed, "nlive": nlive,
               "dlogz_requested": dlogz, "sampler": "dynesty static / multi / rwalk / walks=30",
               "dynesty_version": dynesty.__version__, "python": platform.python_version(),
               "center_bjd": center, "n_data": len(nightly), "n_dim": len(problem.names),
               "data_sha256": digest(ROOT / "data/processed/hd73344_nightly.csv"),
               "lnZ": float(result.logz[-1]), "lnZerr": float(result.logzerr[-1]),
               "weighted_ess": float(1/np.sum(weights**2)), "niter": int(result.niter),
               "ncall": int(result.ncall.sum()), "efficiency_percent": float(result.eff),
               "seconds": time.time()-started,
               "parameters": {name: np.quantile(equal[:, i], [.025, .16, .5, .84, .975]).tolist()
                              for i, name in enumerate(problem.names)}}
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out.with_suffix(".npz"), samples=result.samples, weights=weights,
                        equal_samples=equal, names=np.array(problem.names),
                        logl=result.logl, logwt=result.logwt, logz=result.logz,
                        logzerr=result.logzerr, ncall=result.ncall)
    out.with_suffix(".json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--reuse-audit", action="store_true", help="Read-only verified processed data for concurrent independent runs")
    parser.add_argument("--models", nargs="+", choices=["A", "B", "C", "N"], default=["A", "B", "C"])
    parser.add_argument("--variant", choices=["primary", "broad", "amplitude"], default="primary")
    parser.add_argument("--seeds", nargs="+", type=int, default=[303, 73344])
    parser.add_argument("--nlive", type=int, default=500)
    parser.add_argument("--dlogz", type=float, default=.05)
    args = parser.parse_args()
    if args.nlive < 100 or args.dlogz <= 0:
        parser.error("Use at least 100 live points and positive dlogz")
    if args.audit_only:
        audit()
        return
    gate_path = ROOT / "reports/p03/tutorial_validation.json"
    first_gate = ROOT / "reports/p03/tutorial_likelihood_check.json"
    full_pass = gate_path.exists() and json.loads(gate_path.read_text()).get("approximate_sanity_check_passed", False)
    first_pass = first_gate.exists() and json.loads(first_gate.read_text()).get("first_example_complete", False)
    if not (full_pass or first_pass):
        raise RuntimeError("Reproduce a documented juliet RV example and independently check its likelihood first")
    if args.reuse_audit:
        import pandas as pd
        path = ROOT / "data/processed/hd73344_nightly.csv"
        audit_record = json.loads((ROOT / "reports/p03/data_audit.json").read_text())
        expected = audit_record["hashes"][str(path.relative_to(ROOT))]
        if digest(path) != expected:
            raise RuntimeError("Processed observations changed since audit")
        nightly = pd.read_csv(path)
    else:
        nightly = audit()
    for seed in args.seeds:
        for model in args.models:
            fit(model, args.variant, seed, args.nlive, args.dlogz, nightly)


if __name__ == "__main__":
    main()
