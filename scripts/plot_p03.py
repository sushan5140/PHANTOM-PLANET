"""Generate transparent diagnostics from completed P03 runs (no inference)."""
# ruff: noqa: E402 -- configure a writable cache and headless backend before pyplot.
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
cache = ROOT / ".p03-cache/matplotlib"
cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from scipy.linalg import cholesky, solve_triangular
from scipy.stats import norm

sys.path.insert(0, str(ROOT / "src"))
from phantom.p03_models import RVModel  # noqa: E402

FIG = ROOT / "reports/p03/figures"
COLORS = {"A": "#ad543a", "B": "#327c83", "C": "#584b91", "N": "#747474"}


def save(fig, name):
    fig.savefig(FIG / (name + ".png"), dpi=170, bbox_inches="tight")
    fig.savefig(FIG / (name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def data_figures(exposures, nightly):
    first = exposures.bjd.min()
    fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    ax[0].errorbar(exposures.bjd-first, exposures.mnvel, yerr=exposures.errvel,
                   fmt=".", alpha=.35, color="#666666", label="189 exposures; formal errors")
    ax[0].errorbar(nightly.bjd-first, nightly.mnvel, yerr=nightly.errvel,
                   fmt="o", color="#163f59", label=f"{len(nightly)} nightly means; formal errors")
    ax[0].set(ylabel="Relative RV [m/s]", title="HD 73344: sparse nights, dense exposures")
    ax[0].legend(fontsize=9)
    ax[1].errorbar(exposures.bjd-first, exposures.svalue, yerr=exposures.svalue_err,
                   fmt=".", alpha=.35, color="#327c83")
    ax[1].errorbar(nightly.bjd-first, nightly.svalue, yerr=nightly.svalue_err,
                   fmt="o", color="#163f59")
    ax[1].set(ylabel="S-index", title="Simultaneous activity proxy; not a starspot-only label")
    ax[2].bar(nightly.bjd-first, nightly.n_exposures, width=2, color="#ad543a")
    ax[2].set(xlabel=f"Days since BJD {first:.5f}", ylabel="Exposures/night")
    fig.tight_layout()
    save(fig, "01_measurements_activity")

    frequencies = np.linspace(1/150., 1/2., 12000)
    periods = 1/frequencies
    fig, ax = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    for panel, column, error, label in [(0, "mnvel", "errvel", "Nightly RV"),
                                        (1, "svalue", "svalue_err", "Nightly S-index")]:
        power = LombScargle(nightly.bjd-first, nightly[column], nightly[error]).power(frequencies)
        ax[panel].plot(periods, power, color="#163f59")
        ax[panel].set(ylabel="GLS power", title=label + ": descriptive only; no calibrated FAP")
    window = np.abs(np.exp(2j*np.pi*frequencies[:, None]*(nightly.bjd.to_numpy()-first)).mean(axis=1))**2
    ax[2].plot(periods, window, color="#ad543a")
    ax[2].set(ylabel="Window power", xlabel="Period [days]", title="Nightly sampling window (2–150 day range)")
    for axis in ax:
        for period, label, color in [(9.03, "rotation", "#327c83"),
                                      (15.611, "known b", "#888888"), (65.936, "known c", "#584b91")]:
            axis.axvline(period, color=color, ls="--", alpha=.7, label=label)
        axis.set_xscale("log")
        axis.set_xlim(2, 150)
    ax[2].legend(ncol=3, fontsize=8)
    fig.tight_layout()
    save(fig, "02_descriptive_periodograms_window")

    fig, ax = plt.subplots(figsize=(7, 5))
    scatter = ax.scatter(nightly.svalue, nightly.mnvel,
                         c=nightly.bjd-first, cmap="viridis", s=45)
    ax.errorbar(nightly.svalue, nightly.mnvel, yerr=nightly.errvel,
                xerr=nightly.svalue_err, fmt="none", color="grey", alpha=.4, zorder=0)
    fig.colorbar(scatter, ax=ax, label="Days since first exposure")
    ax.set(xlabel="Nightly S-index", ylabel="Nightly relative RV [m/s]",
           title="Activity/RV association is not causal attribution")
    fig.tight_layout()
    save(fig, "03_activity_association")


def prediction_diagonal(model, theta, t):
    """Marginal conditional latent variance without allocating a grid covariance."""
    from scipy.linalg import cho_factor, cho_solve
    from phantom.p03_models import qp_kernel

    pars, base, planet = model.components(theta, t)
    if "gp_amp" not in pars:
        return base+planet, np.zeros(len(t))
    cross = qp_kernel(t[:, None]-model.t[None, :],
                      *[pars[k] for k in ["gp_amp", "lifetime", "gamma", "rotation"]])
    factor = cho_factor(model.covariance(theta), lower=True)
    _, train_base, train_planet = model.components(theta)
    mean = base+planet+cross@cho_solve(factor, model.y-train_base-train_planet)
    variance = pars["gp_amp"]**2 - np.einsum("ij,ji->i", cross, cho_solve(factor, cross.T))
    return mean, np.maximum(variance, 0)


def fitted_figures(nightly, runs):
    rng = np.random.default_rng(20261017)
    fig, axes = plt.subplots(3, 2, figsize=(14, 11), gridspec_kw={"width_ratios": [2, 1]})
    diagnostic = []
    for row, letter in enumerate(["A", "B", "C"]):
        summary, raw = runs[letter]
        center = summary["center_bjd"]
        x = nightly.bjd.to_numpy()-center
        model = RVModel(x, nightly.mnvel.to_numpy(), nightly.errvel.to_numpy(), letter)
        grid = np.linspace(x.min()-2, x.max()+2, 600)
        equal = raw["equal_samples"]
        samples = equal[rng.choice(len(equal), 300, replace=False)]
        draws, at_data, pp_draws = [], [], []
        for sample in samples:
            mean, variance = prediction_diagonal(model, sample, grid)
            draws.append(mean+rng.normal(size=len(grid))*np.sqrt(variance))
            observed_mean, observed_var = prediction_diagonal(model, sample, x)
            at_data.append(observed_mean)
            jitter = sample[2]
            pp_draws.append(observed_mean+rng.normal(size=len(x))*np.sqrt(observed_var+model.err**2+jitter**2))
        q = np.quantile(draws, [.025, .16, .5, .84, .975], axis=0)
        axis = axes[row, 0]
        axis.fill_between(grid, q[0], q[4], color=COLORS[letter], alpha=.12, label="95% pointwise latent band")
        axis.fill_between(grid, q[1], q[3], color=COLORS[letter], alpha=.25, label="68% pointwise latent band")
        axis.plot(grid, q[2], color=COLORS[letter], lw=1.3)
        axis.errorbar(x, model.y, yerr=model.err, fmt="o", color="#222222", ms=4, label="Nightly data; formal errors")
        axis.set(ylabel="Relative RV [m/s]", title=f"{letter}: " + {"A": "circular Keplerian", "B": "rotation GP", "C": "Keplerian + rotation GP"}[letter])
        axis.set_ylim(-85, 85)
        # Report residuals at training data explicitly as an interpolation diagnostic.
        residual = model.y-np.median(at_data, axis=0)
        axes[row, 1].errorbar(x, residual, yerr=np.sqrt(model.err**2+np.median(equal[:, 2])**2),
                             fmt="o", ms=4, color=COLORS[letter])
        axes[row, 1].axhline(0, color="grey", lw=1)
        axes[row, 1].set(ylabel="Residual [m/s]", title="In-sample residual; errors include median jitter")
        best = raw["samples"][np.argmax(raw["logl"])]
        _, base, planet = model.components(best)
        factor = cholesky(model.covariance(best), lower=True)
        whitened = solve_triangular(factor, model.y-base-planet, lower=True)
        coverage = np.mean((model.y >= np.quantile(pp_draws, .025, axis=0)) &
                           (model.y <= np.quantile(pp_draws, .975, axis=0)))
        diagnostic.append({"model": letter, "residual_rms_ms": float(np.sqrt(np.mean(residual**2))),
                           "in_sample_95_predictive_fraction": float(coverage),
                           "max_likelihood_whitened_rms": float(np.sqrt(np.mean(whitened**2))),
                           "whitened": whitened.tolist(),
                           "warning": "In-sample conditional residuals/coverage; not a held-out validation."})
    for axis in axes[-1]:
        axis.set_xlabel(f"BJD - {center:.5f}")
    axes[0, 0].legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    save(fig, "04_models_and_residuals")
    (ROOT / "reports/p03/predictive_diagnostics.json").write_text(json.dumps(diagnostic, indent=2)+"\n")

    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    for axis, item in zip(ax, diagnostic):
        w = np.sort(item["whitened"])
        expected = norm.ppf((np.arange(len(w))+.5)/len(w))
        axis.scatter(expected, w, color=COLORS[item["model"]])
        axis.plot([-3, 3], [-3, 3], color="grey", ls="--")
        axis.set(title=f"{item['model']}: Cholesky-whitened residuals", xlabel="Normal quantile", ylabel="Whitened residual")
    fig.tight_layout()
    save(fig, "05_whitened_residuals")

    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    for name, axis in zip(["period", "K", "jitter", "rotation", "gp_amp", "lifetime"], axes.flat):
        logarithmic = name in ["jitter", "gp_amp", "lifetime"]
        for letter in ["A", "B", "C"]:
            _, raw = runs[letter]
            names = raw["names"].tolist()
            if name in names:
                values = raw["samples"][:, names.index(name)]
                if logarithmic:
                    values = np.log10(values)
                axis.hist(values, weights=raw["weights"],
                           bins=45, density=True, histtype="step", color=COLORS[letter], label=letter)
        bounds = {"period": (40., 100.), "K": (0., 50.), "jitter": (.1, 50.),
                  "rotation": (5., 15.), "gp_amp": (.1, 100.), "lifetime": (3., 300.)}
        low, high = bounds[name]
        prior_x = np.linspace(low, high, 600)
        if name == "K":
            prior_y = np.ones_like(prior_x)/(high-low)
        elif name == "rotation":
            prior_y = norm.pdf(prior_x, 9.03, 1.3)/(norm.cdf(high, 9.03, 1.3)-norm.cdf(low, 9.03, 1.3))
        elif logarithmic:
            prior_x = np.linspace(np.log10(low), np.log10(high), 600)
            prior_y = np.ones_like(prior_x)/np.log10(high/low)
        else:
            prior_y = 1/(prior_x*np.log(high/low))
        axis.plot(prior_x, prior_y, color="grey", ls=":", lw=1, label="Prior")
        labels = {"period": "Orbital period [days]", "K": "RV semi-amplitude [m/s]",
                  "jitter": "Nightly RV jitter [m/s]", "rotation": "Rotation period [days]",
                  "gp_amp": "GP amplitude [m/s]", "lifetime": "GP coherence time [days]"}
        axis.set(xlabel=labels[name], ylabel="Density per dex" if logarithmic else "Posterior density")
        if logarithmic:
            ticks = ({"jitter": [.1, 1, 10, 50], "gp_amp": [.1, 1, 10, 100],
                      "lifetime": [3, 10, 30, 100, 300]})[name]
            axis.set_xticks(np.log10(ticks), labels=[str(value) for value in ticks])
            axis.set_xlim(np.log10(low), np.log10(high))
        axis.legend(fontsize=8)
    axes[0, 0].axvline(65.936, color="grey", ls="--")
    axes[0, 1].axvline(16.83, color="grey", ls="--")
    axes[1, 0].axvline(9.03, color="grey", ls="--")
    settings = ", ".join(f"{letter}: n={runs[letter][0]['nlive']}, s={runs[letter][0]['seed']}" for letter in ["A", "B", "C"])
    fig.suptitle("Posteriors and dotted priors; dashed vertical lines: published values\n"+settings, fontsize=11)
    fig.tight_layout()
    save(fig, "06_parameter_posteriors")


def convergence_figures(all_summaries):
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for letter in ["A", "B", "C", "N"]:
        matches = [r for r in all_summaries if r["model"] == letter and r["variant"] == "primary"]
        for i, summary in enumerate(matches):
            label = f"{letter}, s{summary['seed']}, n{summary['nlive']}"
            stem = f"primary_{letter}_s{summary['seed']}_n{summary['nlive']}"
            raw = np.load(ROOT / "results/p03" / (stem + ".npz"))
            tail = np.arange(len(raw["logz"])) / len(raw["logz"])
            ax[0].plot(tail, raw["logz"], color=COLORS[letter], alpha=.7,
                        ls=["-", "--", ":"][i % 3], label=label)
            ax[1].errorbar(summary["lnZ"], len(ax[1].lines), xerr=summary["lnZerr"],
                           fmt="o", color=COLORS[letter], label=label)
    ax[0].set(xlim=(.5, 1), xlabel="Fraction of nested samples", ylabel="Accumulated ln Z", title="Evidence accumulation (includes final live points)")
    ax[1].set(xlabel="Final ln Z with sampler error", title="Independent seeds and live-point checks", yticks=[])
    ax[0].legend(fontsize=7)
    ax[1].legend(fontsize=7, loc="best")
    fig.tight_layout()
    save(fig, "07_nested_sampling_diagnostics")


def tutorial_figure():
    import pickle
    import radvel

    data = np.genfromtxt(ROOT / "data/raw/rvs_toi141.dat", dtype=None, encoding="utf-8")
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for nplanet, axis in zip([1, 2], axes):
        with (ROOT / f"results/tutorial/{nplanet}planet/result_for_plot.pkl").open("rb") as stream:
            result = pickle.load(stream)
        samples, summary = result["samples"], result["summary"]
        times = np.linspace(2458365, 2458435, 1000)
        predicted = np.zeros(len(times))
        for p in range(1, nplanet+1):
            period, t0, amplitude = [np.median(samples[f"{name}_p{p}"]) for name in ["P", "t0", "K"]]
            predicted += radvel.kepler.rv_drive(times, [period, t0, 0., np.pi/2, amplitude])
        for inst, color in [("HARPS", "#327c83"), ("FEROS", "#ad543a")]:
            rows = [r for r in data if r[3] == inst]
            axis.errorbar([r[0]-2457000 for r in rows],
                           [r[1]-np.median(samples["mu_"+inst]) for r in rows],
                           yerr=[r[2] for r in rows], fmt="o", ms=4, color=color, label=inst)
        axis.plot(times-2457000, predicted, color="#222222")
        axis.set(ylabel="Offset-corrected RV [m/s]", ylim=(-35, 30),
                 title=f"TOI-141: {nplanet} circular planet(s); lnZ={summary['lnZ']:.2f} +/- {summary['lnZerr']:.2f}")
        axis.legend(fontsize=8)
    axes[-1].set_xlabel("BJD - 2457000")
    fig.tight_layout()
    save(fig, "00_juliet_tutorial")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    exposures = pd.read_csv(ROOT / "data/processed/hd73344_exposures.csv")
    nightly = pd.read_csv(ROOT / "data/processed/hd73344_nightly.csv")
    data_figures(exposures, nightly)
    all_summaries = [json.loads(path.read_text()) for path in sorted((ROOT / "results/p03").glob("*.json"))]
    runs = {}
    for letter in ["A", "B", "C"]:
        candidates = [r for r in all_summaries if r["model"] == letter and r["variant"] == "primary"]
        resolution = max(r["nlive"] for r in candidates)
        summary = min((r for r in candidates if r["nlive"] == resolution), key=lambda r: r["seed"])
        stem = f"primary_{letter}_s{summary['seed']}_n{summary['nlive']}"
        runs[letter] = (summary, np.load(ROOT / "results/p03" / (stem+".npz")))
    fitted_figures(nightly, runs)
    convergence_figures(all_summaries)
    tutorial_figure()
    pd.DataFrame([{k: r[k] for k in ["variant", "model", "seed", "nlive", "lnZ", "lnZerr", "weighted_ess", "ncall", "seconds"]}
                  for r in all_summaries]).to_csv(ROOT / "reports/p03/model_comparison.csv", index=False)
    (ROOT / "reports/p03/all_run_summaries.json").write_text(json.dumps(all_summaries, indent=2)+"\n")
    print(f"Wrote figures to {FIG}")


if __name__ == "__main__":
    main()
