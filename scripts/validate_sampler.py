"""Independent known-answer normalized-evidence check for dynesty."""
import json
from pathlib import Path

import dynesty
import numpy as np
from scipy.stats import norm


def main():
    # p(y=0 | mu)=Normal(0 | mu, 1), mu~Uniform(-10,10).
    exact = float(np.log(norm.cdf(10)-norm.cdf(-10))-np.log(20))
    sampler = dynesty.NestedSampler(lambda x: norm.logpdf(x[0]),
                                    lambda u: 20*u-10, 1, nlive=400,
                                    rstate=np.random.default_rng(2105))
    sampler.run_nested(dlogz=.02, print_progress=False)
    r = sampler.results
    result = {"exact_lnZ": exact, "sampled_lnZ": float(r.logz[-1]),
              "lnZerr": float(r.logzerr[-1]),
              "passed": bool(abs(r.logz[-1]-exact) < 4*r.logzerr[-1]),
              "note": "Known-answer check of normalized evidence, not target scientific validation."}
    root = Path(__file__).resolve().parents[1]
    (root / "reports/p03/sampler_validation.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result), flush=True)
    if not result["passed"]:
        raise RuntimeError("Known-answer evidence check failed")


if __name__ == "__main__":
    main()
