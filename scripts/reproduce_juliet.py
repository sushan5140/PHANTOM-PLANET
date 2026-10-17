"""Official TOI-141 RV example, original priors; dynesty replaces MultiNest.

Run both circular 1- and 2-planet examples before fitting HD 73344.
Source: https://juliet.readthedocs.io/en/latest/tutorials/rvfits.html
"""
import argparse
import hashlib
import json
import os
import pickle
import shutil
from pathlib import Path

import numpy as np


def tutorial_priors(nplanet):
    priors = {}

    def add(name, dist, hyper):
        priors[name] = {"distribution": dist, "hyperparameters": hyper}

    add("P_p1", "normal", [1.007917, 0.000073])
    add("t0_p1", "normal", [2458325.5386, 0.0011])
    for instrument in ["CORALIE14", "CORALIE07", "HARPS", "FEROS"]:
        add("mu_" + instrument, "uniform", [-100., 100.])
        add("sigma_w_" + instrument, "loguniform", [1e-3, 100.])
    add("K_p1", "uniform", [0., 100.])
    add("ecc_p1", "fixed", 0.)
    add("omega_p1", "fixed", 90.)
    if nplanet == 2:
        for name, dist, hyper in [
            ("P_p2", "uniform", [1., 10.]),
            ("t0_p2", "uniform", [2458325., 2458330.]),
            ("K_p2", "uniform", [0., 100.]),
            ("ecc_p2", "fixed", 0.), ("omega_p2", "fixed", 90.),
        ]:
            add(name, dist, hyper)
    return priors


def main():
    root = Path(__file__).resolve().parents[1]
    cache = root / ".p03-cache/matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import juliet

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nlive", type=int, default=300)
    p.add_argument("--seed", type=int, default=141)
    args = p.parse_args()
    summaries = []
    for planets in [1, 2]:
        out = root / "results" / "tutorial" / f"{planets}planet"
        out.mkdir(parents=True, exist_ok=True)
        config = {"nlive": args.nlive, "seed": args.seed+planets,
                  "data_sha256": hashlib.sha256((root / "data/raw/rvs_toi141.dat").read_bytes()).hexdigest(),
                  "priors": tutorial_priors(planets)}
        config_path = out / "run_config.json"
        if config_path.exists() and json.loads(config_path.read_text()) != config:
            raise RuntimeError(f"Tutorial cache has different inputs/settings: {out}")
        if (out / "_dynesty_NS_posteriors.pkl").exists() and not config_path.exists():
            with (out / "result_for_plot.pkl").open("rb") as stream:
                legacy = pickle.load(stream)["summary"]
            if legacy["seed"] != config["seed"] or legacy["nlive"] != config["nlive"]:
                raise RuntimeError(f"Legacy tutorial cache has different settings: {out}")
        config_path.write_text(json.dumps(config, indent=2)+"\n")
        # juliet 2.2.10 otherwise uses a Unix `cp` shell command on Windows.
        shutil.copyfile(root / "data/raw/rvs_toi141.dat", out / "rvs.dat")
        np.random.seed(args.seed + planets)
        dataset = juliet.load(
            priors=tutorial_priors(planets),
            rvfilename=str(root / "data/raw/rvs_toi141.dat"),
            out_folder=str(out),
        )
        results = dataset.fit(
            sampler="dynesty", n_live_points=args.nlive,
            dlogz=0.1, print_progress=False,
            rstate=np.random.default_rng(args.seed + planets),
        )
        samples = results.posteriors["posterior_samples"]
        summary = {"planets": planets, "nlive": args.nlive,
                   "seed": args.seed + planets,
                   "lnZ": float(results.posteriors["lnZ"]),
                   "lnZerr": float(results.posteriors["lnZerr"]),
                   "posterior": {k: np.quantile(samples[k], [.16, .5, .84]).tolist()
                                 for k in ["K_p1"] + (["K_p2", "P_p2"] if planets == 2 else [])}}
        with (out / "result_for_plot.pkl").open("wb") as f:
            pickle.dump({"summary": summary, "samples": samples}, f)
        summaries.append(summary)
        print(json.dumps(summary), flush=True)
    delta = summaries[1]["lnZ"] - summaries[0]["lnZ"]
    period = summaries[1]["posterior"]["P_p2"][1]
    documented_z, documented_error = [-805.302, -690.892], [.027, .023]
    evidence_matches = [bool(abs(row["lnZ"]-value) <= 3*np.hypot(row["lnZerr"], error))
                        for row, value, error in zip(summaries, documented_z, documented_error)]
    period_matches = bool(abs(period-4.76) < .005)
    validation = {"models": summaries, "delta_lnZ_2_minus_1": delta,
                  "documented_delta_lnZ": 114.4,
                  "documented_figure_lnZ": [-805.302, -690.892],
                  "documented_figure_lnZerr": [.027, .023],
                  "lnZ_minus_documented": [summaries[0]["lnZ"]+805.302,
                                            summaries[1]["lnZ"]+690.892],
                  "documented_rv_only_P2_approx_days": 4.76,
                  "P2_minus_documented_approx_days": float(period-4.76),
                  "historical_evidence_within_3sigma": evidence_matches,
                  "historical_period_rounding_match": period_matches,
                  "historical_check_definition": "Numerical consistency within 3 combined evidence errors and 0.005 d of documented P2, not bitwise equality",
                  "exact_historical_fit_reproduced": bool(all(evidence_matches) and period_matches),
                  "approximate_sanity_check_passed": bool(abs(delta-114.4) < 10 and abs(period - 4.76) < .1
                                 and abs(summaries[0]["lnZ"]+805.302) < 5
                                 and abs(summaries[1]["lnZ"]+690.892) < 5),
                  "scope": "Approximate reproduction of first two circular RV examples. Conservative 5-nat absolute/10-nat delta tolerances; not validation to the tiny documentation figure errors. Sampler/version differ and any discrepancy remains unresolved."}
    report = root / "reports/p03/tutorial_validation.json"
    report.write_text(json.dumps(validation, indent=2) + "\n")
    if not validation["approximate_sanity_check_passed"]:
        raise RuntimeError("Tutorial reproduction gate failed; do not fit target")


if __name__ == "__main__":
    main()
