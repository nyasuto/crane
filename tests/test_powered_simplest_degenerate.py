# tests/test_powered_simplest_degenerate.py
import numpy as np

from crane import references as ref
from crane.models.powered_simplest import PoweredSimplestParams, make_powered_simplest
from crane.models.simplest import SimplestParams, make_simplest
from crane.search import find_limit_cycle
from crane.stride import stride


def test_single_stride_matches_passive_at_pushoff_zero():
    gamma = ref.GAMMA_REF
    passive = make_simplest(SimplestParams(gamma=gamma))
    powered = make_powered_simplest(PoweredSimplestParams(gamma=gamma, push_off=0.0))
    x0 = passive.lift(np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT]))
    r_passive = stride(passive, x0)
    r_powered = stride(powered, x0)
    np.testing.assert_allclose(r_powered.x_end, r_passive.x_end, atol=1e-12)


def test_limit_cycle_matches_phase1_at_pushoff_zero():
    gamma = ref.GAMMA_REF
    powered = make_powered_simplest(PoweredSimplestParams(gamma=gamma, push_off=0.0))
    fp = find_limit_cycle(powered, np.array([ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT]))
    assert fp.converged
    np.testing.assert_allclose(fp.y, [0.2003109, -0.1998325], atol=1e-5)
    assert np.max(np.abs(fp.eigenvalues)) < 1.0
