"""Independently check juliet's saved RV likelihood; export portable samples."""
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

from reproduce_juliet import tutorial_priors

ROOT = Path(__file__).resolve().parents[1]


def main():
    data = pd.read_csv(ROOT / "data/raw/rvs_toi141.dat", sep=r"\s+", header=None,
                       names=["t", "y", "err", "instrument"])
    checks = []
    destination = ROOT / "reports/p03/posteriors"
    destination.mkdir(parents=True, exist_ok=True)
    for nplanet in [1, 2]:
        path = ROOT / f"results/tutorial/{nplanet}planet/_dynesty_NS_posteriors.pkl"
        if not path.exists():
            continue
        # Only our own locally generated pickle, never an untrusted downloaded pickle.
        with path.open("rb") as stream:
            saved = pickle.load(stream)
        r = saved["dynesty_output"]
        names = [name for name, prior in tutorial_priors(nplanet).items()
                 if prior["distribution"] != "fixed"]
        differences = []
        relative_differences = []
        matches = []
        ordering_rejections = 0
        for index in np.linspace(0, len(r.samples)-1, 50, dtype=int):
            pars = dict(zip(names, r.samples[index]))
            prediction = np.zeros(len(data))
            for p in range(1, nplanet+1):
                prediction -= pars[f"K_p{p}"]*np.sin(2*np.pi*(data.t-pars[f"t0_p{p}"])/pars[f"P_p{p}"])
            mu = np.array([pars["mu_"+str(inst)] for inst in data.instrument])
            jitter = np.array([pars["sigma_w_"+str(inst)] for inst in data.instrument])
            variance = data.err.to_numpy()**2+jitter**2
            loglike = -.5*np.sum((data.y.to_numpy()-mu-prediction)**2/variance+np.log(2*np.pi*variance))
            if nplanet == 2 and pars["P_p1"] >= pars["P_p2"]:
                # juliet's documented increasing-period constraint returns this sentinel.
                loglike = -1e101
                ordering_rejections += 1
            differences.append(float(loglike-r.logl[index]))
            relative_differences.append(float(abs(loglike-r.logl[index])/(1+abs(r.logl[index]))))
            matches.append(bool(np.isclose(loglike, r.logl[index], atol=1e-5, rtol=1e-10)))
        maximum = float(np.max(np.abs(differences)))
        checks.append({"planets": nplanet, "points_checked": len(differences),
                       "max_abs_lnlike_difference": maximum,
                       "max_relative_lnlike_difference": max(relative_differences),
                       "ordering_rejections_checked": ordering_rejections,
                       "tolerance": "atol=1e-5, rtol=1e-10, including extreme low-likelihood prior samples",
                       "passed": all(matches)})
        np.savez_compressed(destination / f"tutorial_{nplanet}planet.npz", samples=r.samples,
                            weights=np.exp(r.logwt-r.logz[-1]), names=np.array(names),
                            logl=r.logl, logwt=r.logwt, logz=r.logz, logzerr=r.logzerr)
    first_path = ROOT / "results/tutorial/1planet/_dynesty_NS_posteriors.pkl"
    first_lnz = None
    if first_path.exists():
        with first_path.open("rb") as stream:
            first_lnz = float(pickle.load(stream)["lnZ"])
    first_check = next((c for c in checks if c["planets"] == 1), None)
    result = {"checks": checks,
              "first_example_complete": bool(first_check and first_check["passed"] and
                                               first_lnz is not None and abs(first_lnz+805.302) < 5),
              "first_example_lnZ": first_lnz,
              "first_example_documented_lnZ": -805.302,
              "scope": "First documented example completed with approximate evidence agreement and independently checked likelihood. Two-planet extension may continue concurrently with target sampling; no exact evidence equivalence claim."}
    (ROOT / "reports/p03/tutorial_likelihood_check.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result), flush=True)
    if not all(c["passed"] for c in checks):
        raise RuntimeError("Juliet likelihood differs from independent calculation")


if __name__ == "__main__":
    main()
