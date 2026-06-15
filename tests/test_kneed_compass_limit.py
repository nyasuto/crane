import numpy as np

from crane import references_goswami as gref
from crane.models.compass import CompassParams, make_compass
from crane.models.kneed import KneedParams, dynamics_locked, heelstrike_map, make_kneed
from crane.search import find_limit_cycle

# compass (m=5, m_h=10, a=b=0.5) と等価な退化 kneed:
#   大腿質量 5 kg を hip から b_t=0.5 に、脛はほぼ無質量。
#   l_t + l_s = 1 を維持（l_t=0.5, l_s=0.5）。
P_DEG = KneedParams(
    m_h=gref.M_HIP,
    m_t=gref.M_LEG,
    m_s=1e-9,
    l_t=0.5,
    l_s=0.5,
    b_t=gref.B,
    b_s=0.25,
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
C_MODEL = make_compass(P_COMPASS)


def test_locked_dynamics_reduces_to_compass():
    """m_s→0 の locked 相力学が compass 力学と一致（ランダム10状態）。"""
    rng = np.random.default_rng(0)
    for _ in range(10):
        th_st, th_sw = rng.uniform(-0.4, 0.4, 2)
        w_st, w_sw = rng.uniform(-2.0, 2.0, 2)
        x6 = np.array([th_st, th_sw, th_sw, w_st, w_sw, w_sw])
        x4 = np.array([th_st, th_sw, w_st, w_sw])
        d6 = dynamics_locked(0.0, x6, P_DEG)
        d4 = C_MODEL.phases[0].dynamics(0.0, x4)
        assert np.allclose([d6[3], d6[4]], [d4[2], d4[3]], rtol=0, atol=1e-5)


def test_heelstrike_reduces_to_compass():
    """m_s→0 の heel-strike が compass 衝突に一致。"""
    rng = np.random.default_rng(1)
    for _ in range(10):
        th = rng.uniform(-0.3, -0.1)
        w_st, w_sw = rng.uniform(-2.0, -0.2, 2)
        x6 = np.array([th, -th, -th, w_st, w_sw, w_sw])
        x4 = np.array([th, -th, w_st, w_sw])
        got = heelstrike_map(x6, P_DEG)
        want = C_MODEL.phases[0].impact(x4)
        assert np.allclose([got[0], got[1]], [want[0], want[1]], atol=1e-12)
        assert np.allclose([got[3], got[4]], [want[2], want[3]], rtol=0, atol=1e-5)


def test_full_cycle_reduces_to_compass():
    """退化 kneed の全 stride 不動点が compass 不動点 (Phase 2 実測) に一致。

    断面はどちらも 3D (θ_st, θ̇_st, θ̇_sw) で同型。"""
    model = make_kneed(P_DEG)
    fp_c = find_limit_cycle(C_MODEL, np.array(gref.SECTION_GUESS))
    assert fp_c.converged
    fp_k = find_limit_cycle(model, fp_c.y.copy())
    assert fp_k.converged
    assert np.isclose(fp_k.y[0], fp_c.y[0], atol=1e-5)
    assert np.isclose(fp_k.y[1], fp_c.y[1], atol=1e-4)
    assert np.isclose(fp_k.y[2], fp_c.y[2], atol=1e-4)
    # 上位 multiplier も compass と一致（slave 脛方向が 1 本増える）
    mags_k = np.sort(np.abs(fp_k.eigenvalues))[::-1]
    mags_c = np.sort(np.abs(fp_c.eigenvalues))[::-1]
    assert np.allclose(mags_k[:2], mags_c[:2], atol=2e-3)
