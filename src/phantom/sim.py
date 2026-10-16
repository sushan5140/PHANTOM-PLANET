import numpy as np

CASES = ('planet', 'activity', 'both', 'noise')

def simulate(case: str, seed: int = 42, n: int = 90):
    """Toy radial-velocity series; not a physically complete starspot simulator."""
    if case not in CASES:
        raise ValueError(f'case must be one of {CASES}')
    if n < 20:
        raise ValueError('n must be >= 20')
    rng = np.random.default_rng(seed)
    times = np.sort(rng.uniform(0, 140, n))
    errors = rng.uniform(0.7, 1.3, n)
    planet = 2.4*np.sin(2*np.pi*times/17.3 + 0.6) if case in ('planet','both') else np.zeros(n)
    activity = (3.1*np.sin(2*np.pi*times/18.1+0.3) + 0.9*np.sin(4*np.pi*times/18.1+0.1)) * (0.65+0.35*np.cos(2*np.pi*times/98)) if case in ('activity','both') else np.zeros(n)
    rv = planet + activity + rng.normal(0, errors)
    return times, rv, errors, planet, activity

def write_csv(path, case: str, seed: int, n: int):
    t, rv, err, planet, activity = simulate(case, seed, n)
    data = np.column_stack([t,rv,err,planet,activity])
    np.savetxt(path,data,delimiter=',',header='time_days,rv_ms,rv_error_ms,hidden_planet_ms,hidden_activity_ms',comments='')
