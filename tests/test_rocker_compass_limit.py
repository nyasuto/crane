# tests/test_rocker_compass_limit.py
import numpy as np

from crane import references_goswami as gref
from crane.models.compass import CompassParams, make_compass
from crane.models.rocker_compass import RockerCompassParams, make_rocker_compass
from crane.search import find_limit_cycle

P_DEG = RockerCompassParams(
    m=gref.M_LEG,
    m_h=gref.M_HIP,
    c=gref.B,
    rho=1e-9,
    L=gref.A + gref.B,
    R=1e-9,
    gamma=gref.GAMMA_GAIT,
    g=gref.G,
)
P_COMPASS = CompassParams(
    m=gref.M_LEG,
    m_h=gref.M_HIP,
    a=gref.A,
    b=gref.B,
    gamma=gref.GAMMA_GAIT,
    g=gref.G,
)


def test_full_cycle_reduces_to_compass():
    """退化 rocker の全 stride 不動点・固有値が Phase 2 compass に一致。"""
    fp_c = find_limit_cycle(make_compass(P_COMPASS), np.array(gref.SECTION_GUESS))
    assert fp_c.converged
    fp_r = find_limit_cycle(make_rocker_compass(P_DEG), fp_c.y.copy())
    assert fp_r.converged
    assert np.allclose(fp_r.y, fp_c.y, atol=1e-5)
    mags_r = np.sort(np.abs(fp_r.eigenvalues))[::-1]
    mags_c = np.sort(np.abs(fp_c.eigenvalues))[::-1]
    assert np.allclose(mags_r, mags_c, atol=2e-3)
