"""Summarize evidence and seed agreement without selecting a discovery claim."""
import itertools
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance

ROOT = Path(__file__).resolve().parents[1]


def main():
    directory = ROOT / "results/p03"
    runs = [json.loads(p.read_text()) for p in sorted(directory.glob("*.json"))]
    destination = ROOT / "reports/p03/posteriors"
    destination.mkdir(parents=True, exist_ok=True)
    for path in directory.glob("*.npz"):
        shutil.copy2(path, destination / path.name)
    comparisons, consistency = [], []
    for variant in sorted({r["variant"] for r in runs}):
        subset = [r for r in runs if r["variant"] == variant and r["nlive"] == 500]
        for seed in sorted({r["seed"] for r in subset}):
            group = [r for r in subset if r["seed"] == seed]
            for left, right in itertools.combinations(group, 2):
                comparisons.append({"variant": variant, "seed": seed,
                                    "left": left["model"], "right": right["model"],
                                    "delta_lnZ_left_minus_right": left["lnZ"]-right["lnZ"],
                                    "monte_carlo_sigma": np.hypot(left["lnZerr"], right["lnZerr"])})
        for model in sorted({r["model"] for r in subset}):
            group = [r for r in runs if r["variant"] == variant and r["model"] == model]
            for left, right in itertools.combinations(group, 2):
                delta = abs(left["lnZ"]-right["lnZ"])
                sigma = np.hypot(left["lnZerr"], right["lnZerr"])
                a = np.load(directory / f"{variant}_{model}_s{left['seed']}_n{left['nlive']}.npz")
                b = np.load(directory / f"{variant}_{model}_s{right['seed']}_n{right['nlive']}.npz")
                distances = {}
                for i, name in enumerate(a["names"]):
                    if name == "phase":  # Circular phase is not an ordinary real parameter.
                        continue
                    scale = np.sqrt((np.var(a["equal_samples"][:, i])+np.var(b["equal_samples"][:, i]))/2)
                    distances[str(name)] = float(wasserstein_distance(a["equal_samples"][:, i], b["equal_samples"][:, i])/scale)
                consistency.append({"variant": variant, "model": model,
                                    "runs": [[left["seed"], left["nlive"]], [right["seed"], right["nlive"]]],
                                    "abs_delta_lnZ": delta, "combined_sigma": sigma,
                                    "difference_in_sigma": delta/sigma,
                                    "wasserstein_over_pooled_sd": distances,
                                    "diagnostic_pass": bool(delta/sigma < 3 and max(distances.values()) < .3),
                                    "caution": "Agreement is a diagnostic, not proof that all modes were found."})
    preferred = []
    for model in sorted({r["model"] for r in runs if r["variant"] == "primary"}):
        resolution = max(r["nlive"] for r in runs if r["variant"] == "primary" and r["model"] == model)
        selected = [r for r in runs if r["variant"] == "primary" and r["model"] == model and r["nlive"] == resolution]
        checks = [r for r in consistency if r["variant"] == "primary" and r["model"] == model and
                  all(setting[1] == resolution for setting in r["runs"])]
        preferred.append({"model": model, "nlive": resolution, "seeds": [r["seed"] for r in selected],
                          "at_least_two_seeds": len(selected) >= 2,
                          "repeat_diagnostic_pass": len(selected) >= 2 and all(r["diagnostic_pass"] for r in checks)})
    output = {"evidence_differences": comparisons, "repeat_consistency": consistency,
              "preferred_primary_resolution": preferred,
              "all_preferred_primary_checks_pass": all(r["repeat_diagnostic_pass"] for r in preferred),
              "all_primary_repeat_checks_pass": all(r["diagnostic_pass"] for r in consistency if r["variant"] == "primary"),
              "all_runs_ess_at_least_1000": all(r["weighted_ess"] >= 1000 for r in runs),
              "warning": "Nested-sample ESS is weight concentration; not independent MCMC draws. No R-hat claim."}
    (ROOT / "reports/p03/convergence_and_comparisons.json").write_text(json.dumps(output, indent=2)+"\n")
    pd.DataFrame(comparisons).to_csv(ROOT / "reports/p03/evidence_differences.csv", index=False)
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
