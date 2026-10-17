"""Independent fast circular-Gaussian fit to the same TOI-141 tutorial priors.

Checks current juliet evidence/periods; does not explain historical doc differences.
"""
import argparse
import json
from pathlib import Path

import dynesty
import numpy as np
import pandas as pd
from scipy.special import ndtri

from reproduce_juliet import tutorial_priors

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", nargs="+", type=int, choices=[1, 2], default=[1, 2])
    parser.add_argument("--nlive", type=int, default=500)
    parser.add_argument("--sample", choices=["rwalk", "unif"], default="rwalk")
    parser.add_argument("--seed-base", type=int, default=9141)
    parser.add_argument("--tag", default="")
    args = parser.parse_args()
    if args.nlive < 100 or any(not (ch.isalnum() or ch in "_-") for ch in args.tag):
        parser.error("Use nlive>=100 and an alphanumeric output tag")
    suffix = "_"+args.tag if args.tag else ""
    (ROOT / "reports/p03/posteriors").mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(ROOT / "data/raw/rvs_toi141.dat", sep=r"\s+", header=None,
                       names=["t", "y", "err", "instrument"])
    t, y, err = [data[name].to_numpy() for name in ["t", "y", "err"]]
    runs = []
    for planets in args.models:
        priors = {name: p for name, p in tutorial_priors(planets).items() if p["distribution"] != "fixed"}
        names = list(priors)
        indices = {name: i for i, name in enumerate(names)}
        mu_indices = np.array([indices["mu_"+name] for name in data.instrument])
        jitter_indices = np.array([indices["sigma_w_"+name] for name in data.instrument])

        def transform(u):
            theta = np.empty(len(names))
            for i, prior in enumerate(priors.values()):
                lo, hi = prior["hyperparameters"]
                dist = prior["distribution"]
                if dist == "normal":
                    theta[i] = lo+hi*ndtri(u[i])
                elif dist == "uniform":
                    theta[i] = lo+(hi-lo)*u[i]
                elif dist == "loguniform":
                    theta[i] = lo*(hi/lo)**u[i]
                else:
                    raise ValueError(dist)
            return theta

        def likelihood(theta):
            if planets == 2 and theta[indices["P_p1"]] >= theta[indices["P_p2"]]:
                return -1e101
            prediction = theta[mu_indices].copy()
            variance = err**2+theta[jitter_indices]**2
            for planet in range(1, planets+1):
                period, t0, amplitude = [theta[indices[f"{key}_p{planet}"]] for key in ["P", "t0", "K"]]
                prediction -= amplitude*np.sin(2*np.pi*(t-t0)/period)
            return -.5*np.sum((y-prediction)**2/variance+np.log(2*np.pi*variance))

        walk_settings = {"walks": 40} if args.sample == "rwalk" else {}
        sampler = dynesty.NestedSampler(likelihood, transform, len(names), nlive=args.nlive,
                                        sample=args.sample, bound="multi", **walk_settings,
                                        rstate=np.random.default_rng(args.seed_base+planets))
        sampler.run_nested(dlogz=.05, print_progress=False)
        result = sampler.results
        weights = np.exp(result.logwt-result.logz[-1])
        equal = dynesty.utils.resample_equal(result.samples, weights,
                                             rstate=np.random.default_rng(args.seed_base+10000+planets))
        summary = {"planets": planets, "nlive": args.nlive, "seed": args.seed_base+planets,
                   "sample": args.sample, "walks": 40 if args.sample == "rwalk" else None,
                   "dlogz": .05, "lnZ": float(result.logz[-1]),
                   "lnZerr": float(result.logzerr[-1]),
                   "P2": (np.quantile(equal[:, indices["P_p2"]], [.16, .5, .84]).tolist() if planets == 2 else None)}
        runs.append(summary)
        np.savez_compressed(ROOT / f"reports/p03/posteriors/tutorial_direct{suffix}_{planets}planet.npz",
                            samples=result.samples, weights=weights, names=np.array(names),
                            logl=result.logl, logz=result.logz, logwt=result.logwt,
                            logzerr=result.logzerr)
        (ROOT / f"reports/p03/tutorial_direct{suffix}_{planets}planet.json").write_text(json.dumps(summary, indent=2)+"\n")
        print(json.dumps(summary), flush=True)
    available = [json.loads((ROOT / f"reports/p03/tutorial_direct{suffix}_{p}planet.json").read_text())
                 for p in [1, 2] if (ROOT / f"reports/p03/tutorial_direct{suffix}_{p}planet.json").exists()]
    (ROOT / f"reports/p03/tutorial_direct{suffix}_validation.json").write_text(json.dumps({"models": available}, indent=2)+"\n")


if __name__ == "__main__":
    main()
