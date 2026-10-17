"""Normalized Gaussian likelihood and proper priors for the P03 comparison.

    A: circular Keplerian + offset + slope + nightly jitter
    B: rotation quasi-periodic GP + offset + slope + nightly jitter
    C: both of the above, with identical shared priors

The circular orbit is a Keplerian special case. Its phase has a uniform prior
conditional on each period, avoiding a fixed t0 window that biases phase.
No activity proxy is used to select the RV prior or subtracted from the RV.
"""
from dataclasses import dataclass

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.special import ndtr, ndtri

ROTATION_CDF_LOW = ndtr((5.-9.03)/1.3)
ROTATION_CDF_HIGH = ndtr((15.-9.03)/1.3)


def loguniform(u, low, high):
    return low * (high / low) ** u


def qp_kernel(dt, amplitude, lifetime, gamma, rotation):
    return amplitude ** 2 * np.exp(
        -0.5 * (dt / lifetime) ** 2 - gamma * np.sin(np.pi * dt / rotation) ** 2
    )


def circular_rv(t, period, phase, amplitude):
    """phase is conjunction phase in [0,1); stellar RV, omega=90 degrees."""
    return -amplitude * np.sin(2 * np.pi * (np.asarray(t) / period - phase))


@dataclass
class RVModel:
    t: np.ndarray
    y: np.ndarray
    err: np.ndarray
    model: str
    variant: str = "primary"

    def __post_init__(self):
        if self.model not in ("A", "B", "C", "N"):
            raise ValueError("Unknown model")
        if self.variant not in ("primary", "broad", "amplitude"):
            raise ValueError("Unknown prior variant")
        self.t = np.asarray(self.t)
        self.y = np.asarray(self.y)
        self.err = np.asarray(self.err)
        if not (self.t.shape == self.y.shape == self.err.shape) or np.any(self.err <= 0):
            raise ValueError("Inconsistent observations or nonpositive errors")
        self.dt = self.t[:, None] - self.t[None, :]
        self.names = ["offset", "slope", "jitter"]
        if self.model in ("A", "C"):
            self.names += ["period", "phase", "K"]
        if self.model in ("B", "C"):
            self.names += ["gp_amp", "lifetime", "gamma", "rotation"]

    def prior(self, u):
        high = 100. if self.variant == "amplitude" else 50.
        values = [-100 + 200 * u[0], -.2 + .4 * u[1], loguniform(u[2], .1, high)]
        index = 3
        if self.model in ("A", "C"):
            lo, hi = (2., 150.) if self.variant == "broad" else (40., 100.)
            values += [loguniform(u[index], lo, hi), u[index + 1], high * u[index + 2]]
            index += 3
        if self.model in ("B", "C"):
            rotation = (loguniform(u[index + 3], 2., 100.) if self.variant == "broad"
                        else 9.03 + 1.3*ndtri(ROTATION_CDF_LOW + u[index + 3] *
                                               (ROTATION_CDF_HIGH-ROTATION_CDF_LOW)))
            values += [loguniform(u[index], .1, 100.),
                       loguniform(u[index + 1], 3., 300.),
                       loguniform(u[index + 2], .1, 10.), rotation]
        return np.array(values)

    def components(self, theta, t=None):
        pars = dict(zip(self.names, theta))
        times = self.t if t is None else np.asarray(t)
        base = pars["offset"] + pars["slope"] * times
        planet = (circular_rv(times, pars["period"], pars["phase"], pars["K"])
                  if "K" in pars else np.zeros_like(times))
        return pars, base, planet

    def covariance(self, theta):
        pars, _, _ = self.components(theta)
        cov = np.diag(self.err ** 2 + pars["jitter"] ** 2)
        if "gp_amp" in pars:
            cov += qp_kernel(self.dt, pars["gp_amp"], pars["lifetime"],
                             pars["gamma"], pars["rotation"])
        return cov

    def loglike(self, theta):
        _, base, planet = self.components(theta)
        residual = self.y - base - planet
        cov = self.covariance(theta)
        try:
            factor = cho_factor(cov, lower=True, check_finite=False)
            alpha = cho_solve(factor, residual, check_finite=False)
            logdet = 2 * np.log(np.diag(factor[0])).sum()
            return -.5 * (residual @ alpha + logdet + len(residual) * np.log(2 * np.pi))
        except np.linalg.LinAlgError:
            return -np.inf

    def predict(self, theta, t):
        """Conditional latent mean/covariance; excludes future measurement noise."""
        pars, base, planet = self.components(theta, t)
        if "gp_amp" not in pars:
            return base + planet, np.zeros((len(t), len(t))), np.zeros(len(t))
        kernel_args = [pars[k] for k in ["gp_amp", "lifetime", "gamma", "rotation"]]
        cross = qp_kernel(np.asarray(t)[:, None] - self.t[None, :], *kernel_args)
        prior_cov = qp_kernel(np.asarray(t)[:, None] - np.asarray(t)[None, :], *kernel_args)
        factor = cho_factor(self.covariance(theta), lower=True)
        _, train_base, train_planet = self.components(theta)
        gp_mean = cross @ cho_solve(factor, self.y - train_base - train_planet)
        cov = prior_cov - cross @ cho_solve(factor, cross.T)
        return base + planet + gp_mean, cov, gp_mean
