# PHANTOM-PLANET P02 — exploratory quality + catalog cross-check

Input: KOA K2 Exoplanet Diversity RV Activity CSV, 2332 observations and 84 unique stars.

Outputs: all_stars_quality.csv; three diagnostic figures.

`analyze.py` generates quality summaries, Pearson/Spearman RV vs S-value correlations and *unweighted* Lomb–Scargle exploratory peaks over 2–150 days. These are **NOT** exoplanet detections, calibrated significances, GP fits, or juliet results. Naive correlation p-values ignore autocorrelation, planets, and other confounders. No RV zero-point offsets or activity GP were fitted.

Catalog cross-check:
- EPIC 220383386 = HD 3167, known multi-planet system (Christiansen et al. 2017, AJ 154 122).
- EPIC 228724232 = WASP-107 (SIMBAD identifier cross match).
- EPIC 212178066 = HD 73344, K2 planet-candidate literature (Yu et al. 2018; MNRAS 508, 195).

Important: This KOA dataset is primarily follow-up of known planet candidates/hosts, **not** a labeled starspot-only negative dataset. Clean negative controls and robust multi-model fitting require independent work. juliet not installed here; no fits are claimed.

Reproduction: `python analyze.py` with original CSV in /mnt/data or update input path.

KOA original: https://koa.ipac.caltech.edu/Datasets/K2ExoDiversity/K2_RV_Activity.php ; DOI 10.26135/KOA5
