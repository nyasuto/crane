# tests/test_powered_rocker_degenerate.py
import numpy as np

from crane.models.rocker_compass import RockerCompassParams, make_rocker_compass
from crane.models.powered_rocker_compass import (
    PoweredRockerCompassParams,
    make_powered_rocker_compass,
)
from crane.search import find_limit_cycle
from crane.stride import stride

NOMINAL = dict(m=1.0, m_h=0.0, c=0.37, rho=0.32, L=1.0, R=0.3, gamma=0.030, g=9.81)
GUESS = np.array([0.30, -1.0, -0.40])


def test_single_stride_matches_passive_at_pushoff_zero():
    passive = make_rocker_compass(RockerCompassParams(**NOMINAL))
    powered = make_powered_rocker_compass(PoweredRockerCompassParams(**NOMINAL, push_off=0.0))
    x0 = passive.lift(np.array([0.30844, -1.26256, -0.87914]))
    np.testing.assert_allclose(stride(powered, x0).x_end, stride(passive, x0).x_end, atol=1e-10)


def test_limit_cycle_matches_phase35_at_pushoff_zero():
    powered = make_powered_rocker_compass(PoweredRockerCompassParams(**NOMINAL, push_off=0.0))
    fp = find_limit_cycle(powered, GUESS)
    assert fp.converged
    np.testing.assert_allclose(fp.y, [0.30844, -1.26256, -0.87914], atol=1e-4)
    assert np.max(np.abs(fp.eigenvalues)) < 1.0
