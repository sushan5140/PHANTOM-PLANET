"""Scientific failure-mode tests, independent likelihood and data checks."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import multivariate_normal

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from phantom.p03_data import nightly_bin, read_koa  # noqa: E402
from phantom.p03_models import RVModel, circular_rv, qp_kernel  # noqa: E402


def test_gaussian_normalization_and_gp_likelihood():
    t = np.array([0., .2, 2.1, 9.])
    y = np.array([1., 2., -3., 1.])
    model = RVModel(t, y, np.ones(4)*.7, "C")
    theta = model.prior(np.full(len(model.names), .5))
    _, base, planet = model.components(theta)
    expected = multivariate_normal.logpdf(y, mean=base+planet, cov=model.covariance(theta))
    assert model.loglike(theta) == pytest.approx(expected, abs=1e-10)


def test_qp_is_positive_definite_and_zero_amplitude_limit():
    t = np.random.default_rng(4).uniform(0, 300, 50)
    matrix = qp_kernel(t[:, None]-t[None, :], 12., 20., 2., 9.03)
    assert np.linalg.eigvalsh(matrix).min() > -1e-10
    np.testing.assert_allclose(np.diag(matrix), 144.)
    a = RVModel(t, np.sin(t), np.ones(50), "A")
    c = RVModel(t, np.sin(t), np.ones(50), "C")
    pa = a.prior(np.full(6, .5))
    pc = np.r_[pa, [0., 20., 2., 9.03]]
    assert a.loglike(pa) == pytest.approx(c.loglike(pc))


def test_circular_keplerian_matches_radvel():
    import radvel

    t = np.linspace(0, 300, 100)
    period, phase, amplitude = 66., .35, 17.
    # For omega=pi/2, circular transit and periastron reference coincide.
    expected = radvel.kepler.rv_drive(t, [period, period*phase, 0., np.pi/2, amplitude])
    np.testing.assert_allclose(circular_rv(t, period, phase, amplitude), expected, atol=1e-10)


def test_shared_priors_and_uniform_phase():
    args = (np.arange(3.), np.ones(3), np.ones(3))
    for variant in ["primary", "broad", "amplitude"]:
        a, b, c = [RVModel(*args, name, variant) for name in ["A", "B", "C"]]
        u = np.linspace(.15, .85, 10)
        aa = a.prior(u[:6])
        bb = b.prior(np.r_[u[:3], u[6:]])
        cc = c.prior(u)
        np.testing.assert_array_equal(aa, cc[:6])
        np.testing.assert_array_equal(bb[:3], cc[:3])
        np.testing.assert_array_equal(bb[3:], cc[6:])
        assert cc[4] == u[4]


def test_rotation_prior_matches_normalized_truncated_normal():
    from scipy.stats import truncnorm

    model = RVModel(np.arange(3.), np.ones(3), np.ones(3), "B")
    for probability in [.001, .1, .5, .9, .999]:
        u = np.full(7, .5)
        u[-1] = probability
        expected = truncnorm.ppf(probability, (5-9.03)/1.3, (15-9.03)/1.3,
                                 loc=9.03, scale=1.3)
        assert model.prior(u)[-1] == pytest.approx(expected, abs=1e-12)


def test_nightly_errors_do_not_include_fictitious_extra_points():
    data = pd.DataFrame({"bjd": [2450000.6, 2450000.7, 2450001.6],
                         "mnvel": [2., 6., 1.], "errvel": [1., 2., 1.],
                         "svalue": [.2, .21, .19], "svalue_err": [.001]*3})
    nightly = nightly_bin(data)
    assert len(nightly) == 2
    assert nightly.iloc[0].mnvel == pytest.approx(2.8)
    assert nightly.iloc[0].errvel == pytest.approx(np.sqrt(.8))
    assert nightly.n_exposures.sum() == 3


def test_hawaii_night_is_not_split_at_julian_noon_boundary():
    data = pd.DataFrame({"bjd": [2450000.95, 2450001.05, 2450001.95],
                         "mnvel": [1., 3., 2.], "errvel": [1., 1., 1.],
                         "svalue": [.2, .2, .2], "svalue_err": [.001]*3})
    nightly = nightly_bin(data)
    assert len(nightly) == 2
    assert nightly.iloc[0].n_exposures == 2
    assert nightly.iloc[0].mnvel == pytest.approx(2.)


def test_bad_and_duplicate_archive_rows_are_rejected(tmp_path):
    data = pd.DataFrame({"star": ["x", "x"], "observation_id": ["a", "a"],
                         "bjd": [1., 2.], "mnvel": [1., 2.], "errvel": [1., 1.],
                         "svalue": [.2, .2], "svalue_err": [.01, .01]})
    path = tmp_path / "data.csv"
    data.to_csv(path, index=False)
    with pytest.raises(ValueError, match="Duplicate"):
        read_koa(path)
    data.loc[1, "observation_id"] = "b"
    data.loc[1, "errvel"] = 0
    data.to_csv(path, index=False)
    with pytest.raises(ValueError, match="Non-positive"):
        read_koa(path)


def test_gp_prediction_agrees_with_conditional_gaussian():
    t = np.arange(5.)
    model = RVModel(t, np.sin(t), np.ones(5), "B")
    theta = model.prior(np.full(7, .5))
    mean, cov, gpmean = model.predict(theta, t)
    pars, base, planet = model.components(theta)
    k = model.covariance(theta)-np.diag(model.err**2+pars["jitter"]**2)
    expected = k @ np.linalg.solve(model.covariance(theta), model.y-base-planet)
    np.testing.assert_allclose(gpmean, expected)
    np.testing.assert_allclose(mean, base+planet+expected)
    assert np.linalg.eigvalsh(cov).min() > -1e-9
