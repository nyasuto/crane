import numpy as np

from crane import references as ref
from crane.models.simplest import SimplestParams, make_simplest
from crane.search import find_limit_cycle

MODEL = make_simplest(SimplestParams(gamma=ref.GAMMA_REF))


def test_period1_is_also_fixed_point_of_s2():
    """period-1 不動点は S² の不動点でもあり、multiplier は2乗になる。"""
    fp1 = find_limit_cycle(MODEL, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT]))
    fp2 = find_limit_cycle(MODEL, fp1.y, n_strides=2)
    assert fp2.converged
    assert np.allclose(fp2.y, fp1.y, atol=1e-9)
    assert np.isclose(
        np.max(np.abs(fp2.eigenvalues)), np.max(np.abs(fp1.eigenvalues)) ** 2, rtol=1e-3
    )
