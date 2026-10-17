"""Strict ingestion and fixed, reproducible nightly aggregation for KOA RVs."""
from pathlib import Path

import numpy as np
import pandas as pd

TARGET = "epic212178066"
FIELDS = ["bjd", "mnvel", "errvel", "svalue", "svalue_err"]


def read_koa(path, require_positive_errors=True):
    data = pd.read_csv(path)
    needed = {"star", "observation_id", *FIELDS}
    if not needed.issubset(data.columns):
        raise ValueError(f"Missing columns: {needed - set(data.columns)}")
    if data.duplicated(["star", "observation_id"]).any():
        raise ValueError("Duplicate star/observation_id; resolve before fitting")
    for field in FIELDS:
        data[field] = pd.to_numeric(data[field], errors="raise")
    if not np.isfinite(data[FIELDS].to_numpy()).all():
        raise ValueError("Non-finite measurement; explicit exclusion audit required")
    if require_positive_errors and (data[["errvel", "svalue_err"]] <= 0).any().any():
        raise ValueError("Non-positive measurement uncertainty")
    return data


def nightly_bin(target):
    """BJD bins floor(BJD - 0.5), boundary near UTC midnight / Hawaii afternoon.

    Formal errors alone are propagated here. A free nightly jitter in every
    model accounts for unaveraged short-term noise; it is NOT divided by sqrt(n).
    The short-term covariance within a night is not inferable from this summary.
    """
    rows = []
    for night, group in target.groupby(np.floor(target.bjd - 0.5).astype(int)):
        if group.bjd.max()-group.bjd.min() > .5:
            raise ValueError("A night bin spans over 12 hours; inspect dates and grouping")
        w = 1 / group.errvel.to_numpy() ** 2
        sw = 1 / group.svalue_err.to_numpy() ** 2
        rows.append({"night": night, "n_exposures": len(group),
                     "bjd": np.average(group.bjd, weights=w),
                     "mnvel": np.average(group.mnvel, weights=w),
                     "errvel": np.sqrt(1 / w.sum()),
                     "svalue": np.average(group.svalue, weights=sw),
                     "svalue_err": np.sqrt(1 / sw.sum()),
                     "rv_scatter": group.mnvel.std(ddof=1),
                     "duration_days": group.bjd.max() - group.bjd.min()})
    return pd.DataFrame(rows).sort_values("bjd").reset_index(drop=True)


def prepare(root):
    root = Path(root)
    data = read_koa(root / "data/raw/koa_k2_rv_activity.csv", require_positive_errors=False)
    target = data.loc[data.star == TARGET].sort_values("bjd").copy()
    if target.empty:
        raise ValueError("HD 73344 target not found")
    if (target[["errvel", "svalue_err"]] <= 0).any().any():
        raise ValueError("Non-positive target uncertainty; explicit audit required")
    out = root / "data/processed"
    out.mkdir(parents=True, exist_ok=True)
    target.to_csv(out / "hd73344_exposures.csv", index=False)
    nightly = nightly_bin(target)
    nightly.to_csv(out / "hd73344_nightly.csv", index=False)
    return data, target, nightly
