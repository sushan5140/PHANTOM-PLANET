# Imported measurements and audit

`raw/` contains unmodified public contributed KOA numeric tables and metadata,
plus the official TOI-141 `juliet` tutorial measurements. See
`scripts/fetch_p03_data.py` for exact URLs and expected SHA-256 digests. Fetching
refuses changed content rather than silently altering the scientific inputs.
No raw HIRES FITS spectra were downloaded or re-reduced.

`p02_audit/phantom_p02/` contains the six original artifacts extracted from the
user-supplied `phantom_p02_audit.zip` without executing their script. The original
zip is not needed for numerical reproduction because these contents are retained.

`processed/hd73344_exposures.csv` retains every selected target exposure.
`processed/hd73344_nightly.csv` contains fixed floor(BJD - 0.5) nightly averages,
formal propagated errors, exposure counts, duration and within-night RV scatter.
Formal nightly errors understate remaining stellar variability; every fitted
model includes an additional nightly jitter. See reports/p03/METHODS.md.

All hashes and counts are in reports/p03/data_audit.json. A nonpositive S error
on one unrelated star is retained and counted globally. Target errors are all
positive. RV units are m/s; BJD in days; S-index dimensionless. Archive BJD time
standard beyond its label is not independently established.

The original CSV RA is labeled hours in its README; the JSON exposes degrees.
The numeric RV/activity join uses star and observation ID, not coordinate units.

KOA acknowledgement and project PI attribution are in reports/p03/METHODS.md.
The P02 audit has no separate license notice; its provenance is the user-provided
archive. No blanket new license is asserted over imported third-party data.
