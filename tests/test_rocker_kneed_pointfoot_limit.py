# tests/test_rocker_kneed_pointfoot_limit.py
import numpy as np

from crane import references_kneed as kref
from crane.models.kneed import KneedParams, make_kneed
from crane.models.rocker_kneed import RockerKneedParams, make_rocker_kneed
from crane.search import find_limit_cycle

KW = dict(
    m_h=kref.M_HIP,
    m_t=kref.M_THIGH,
    m_s=kref.M_SHANK,
    l_t=kref.L_THIGH,
    l_s=kref.L_SHANK,
    b_t=kref.B_THIGH,
    b_s=kref.B_SHANK,
    gamma=kref.GAMMA_GAIT,
    g=kref.G,
)


def test_full_cycle_reduces_to_pointfoot_kneed():
    """R→0 の全 stride 不動点・固有値が Phase 3 検証済み点足 kneed に一致。"""
    fp_k = find_limit_cycle(make_kneed(KneedParams(**KW)), np.array(kref.SECTION_GUESS))
    assert fp_k.converged
    fp_r = find_limit_cycle(make_rocker_kneed(RockerKneedParams(R=1e-9, **KW)), fp_k.y.copy())
    assert fp_r.converged
    assert np.allclose(fp_r.y, fp_k.y, atol=1e-5)
    mags_r = np.sort(np.abs(fp_r.eigenvalues))[::-1]
    mags_k = np.sort(np.abs(fp_k.eigenvalues))[::-1]
    assert np.allclose(mags_r, mags_k, atol=2e-3)
