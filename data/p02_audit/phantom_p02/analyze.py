import pandas as pd,numpy as np
from pathlib import Path
from scipy.stats import pearsonr,spearmanr
from scipy.signal import lombscargle
import matplotlib.pyplot as plt
src=Path('/mnt/data/KOA_K2ExoDiversity_RV_Actvitiy_export (1).csv')
d=pd.read_csv(src); out=Path('/mnt/data/phantom_p02')
rows=[]
for star,g in d.groupby('Star'):
 g=g.sort_values('BJD (days)'); rv=g['RV (m/s)'].to_numpy(); s=g['S-value'].to_numpy(); t=g['BJD (days)'].to_numpy(); e=g['RV Error (m/s)'].to_numpy()
 r,p=(pearsonr(rv,s) if len(g)>=3 and np.std(rv)>0 and np.std(s)>0 else (np.nan,np.nan))
 rho,ps=(spearmanr(rv,s) if len(g)>=3 else (np.nan,np.nan))
 rows.append(dict(star=star,n=len(g),span_days=t.max()-t.min(),rv_std_ms=np.std(rv,ddof=1),median_rv_error_ms=np.median(e),s_std=np.std(s,ddof=1),rv_s_pearson_r=r,pearson_naive_p=p,rv_s_spearman_rho=rho, spearman_naive_p=ps, s_error_nonpositive=int((g['S-value Error']<=0).sum())))
res=pd.DataFrame(rows).sort_values('n',ascending=False);res.to_csv(out/'all_stars_quality.csv',index=False)
print(res.head(10).round(4).to_string(index=False))
for star in ['epic220383386','epic212178066','epic228724232']:
 g=d[d.Star==star].sort_values('BJD (days)');t=g['BJD (days)'].to_numpy(); y=g['RV (m/s)'].to_numpy();s=g['S-value'].to_numpy(); err=g['RV Error (m/s)'].to_numpy()
 f=np.linspace(1/150,1/2,12000); w=2*np.pi*f
 pow=lombscargle(t-t.min(),y-y.mean(),w,normalize=True)
 idx=np.argmax(pow); best=1/f[idx]
 fig,ax=plt.subplots(1,3,figsize=(13,3.4));ax[0].errorbar(t-t.min(),y,yerr=err,fmt='.',alpha=.7);ax[0].set(xlabel='Days since first observation',ylabel='RV (m/s)',title=f'{star}: radial velocity');ax[1].scatter(s,y,s=13,alpha=.7);ax[1].set(xlabel='S-value',ylabel='RV (m/s)',title='Activity vs RV');ax[2].plot(1/f,pow);ax[2].axvline(best,ls='--',alpha=.5);ax[2].set(xlabel='Period (days)',ylabel='Normalized power',title=f'Unweighted peak {best:.2f} d');fig.tight_layout();fig.savefig(out/f'{star}_diagnostics.png',dpi=150);plt.close(fig)
 print(star,'unweighted strongest period',round(best,3),'power',round(float(pow[idx]),3))
(out/'README.md').write_text('''# PHANTOM-PLANET P02 — exploratory quality + catalog cross-check\n\nInput: KOA K2 Exoplanet Diversity RV Activity CSV, 2332 observations and 84 unique stars.\n\nOutputs: all_stars_quality.csv; three diagnostic figures.\n\n`analyze.py` generates quality summaries, Pearson/Spearman RV vs S-value correlations and *unweighted* Lomb–Scargle exploratory peaks over 2–150 days. These are **NOT** exoplanet detections, calibrated significances, GP fits, or juliet results. Naive correlation p-values ignore autocorrelation, planets, and other confounders. No RV zero-point offsets or activity GP were fitted.\n\nCatalog cross-check:\n- EPIC 220383386 = HD 3167, known multi-planet system (Christiansen et al. 2017, AJ 154 122).\n- EPIC 228724232 = WASP-107 (SIMBAD identifier cross match).\n- EPIC 212178066 = HD 73344, K2 planet-candidate literature (Yu et al. 2018; MNRAS 508, 195).\n\nImportant: This KOA dataset is primarily follow-up of known planet candidates/hosts, **not** a labeled starspot-only negative dataset. Clean negative controls and robust multi-model fitting require independent work. juliet not installed here; no fits are claimed.\n\nReproduction: `python analyze.py` with original CSV in /mnt/data or update input path.\n\nKOA original: https://koa.ipac.caltech.edu/Datasets/K2ExoDiversity/K2_RV_Activity.php ; DOI 10.26135/KOA5\n''')
