import numpy as np
from scipy.signal import lombscargle

def periodic_candidate(t,rv,period_min=3.0,period_max=45.0):
    """Deliberately weak baseline: strongest Lomb-Scargle period; no FAP/probability claim."""
    t,rv = np.asarray(t), np.asarray(rv)
    frequencies=np.linspace(1/period_max,1/period_min,2500)
    powers=lombscargle(t,rv-np.mean(rv),2*np.pi*frequencies,normalize=True)
    idx=int(np.argmax(powers))
    return {'period_days':float(1/frequencies[idx]),'normalized_power':float(powers[idx])}
