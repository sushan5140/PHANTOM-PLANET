from pathlib import Path
import json
import numpy as np

def fit_one(csv_path: str, output_dir: str, period_hint=17.3, n_live_points=150):
    try:
        import juliet
    except ImportError as exc:
        raise RuntimeError('juliet not installed; install from requirements-juliet.txt') from exc
    arr = np.genfromtxt(csv_path, delimiter=',', names=True)
    t=np.asarray(arr['time_days']); y=np.asarray(arr['rv_ms']); dy=np.asarray(arr['rv_error_ms'])
    if len(t)<20 or not(np.isfinite(t).all() and np.isfinite(y).all() and np.isfinite(dy).all() and np.all(dy>0)):
        raise ValueError('Invalid input RV data')
    priors={
        'P_p1':{'distribution':'uniform','hyperparameters':[max(1.,period_hint*.7),period_hint*1.3]},
        't0_p1':{'distribution':'uniform','hyperparameters':[float(np.min(t)),float(np.min(t)+period_hint)]},
        'K_p1':{'distribution':'uniform','hyperparameters':[0.,12.]},
        'ecc_p1':{'distribution':'fixed','hyperparameters':0.},
        'omega_p1':{'distribution':'fixed','hyperparameters':90.},
        'mu_SIM':{'distribution':'uniform','hyperparameters':[-15.,15.]},
        'sigma_w_SIM':{'distribution':'loguniform','hyperparameters':[0.01,10.]},
    }
    ds=juliet.load(priors=priors,t_rv={'SIM':t},y_rv={'SIM':y},yerr_rv={'SIM':dy},out_folder=output_dir)
    results=ds.fit(sampler='dynesty',n_live_points=n_live_points)
    summary={'log_evidence':float(results.posteriors['lnZ']), 'warning':'Single-planet model only. Evidence is not a planet probability or a planet/activity comparison.'}
    Path(output_dir).mkdir(parents=True,exist_ok=True)
    (Path(output_dir)/'phantom_summary.json').write_text(json.dumps(summary,indent=2))
    return summary
