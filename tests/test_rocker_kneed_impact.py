import numpy as np

from crane import references_kneed as kref
from crane.models.kneed import KneedParams
from crane.models.kneed import dynamics_unlocked as kn_unlocked
from crane.models.kneed import dynamics_locked as kn_locked
from crane.models.kneed import kneestrike_map as kn_ks
from crane.models.kneed import heelstrike_map as kn_hs
from crane.models.rocker_kneed import (
    RockerKneedParams,
    dynamics_unlocked,
    dynamics_locked,
    heelstrike_map,
    kneestrike_map,
    kinetic_energy,
)

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
P = RockerKneedParams(R=0.2, **KW)


def test_kneestrike_locks_and_dissipates():
    x_pre = np.array([0.1, -0.2, -0.35, -1.0, -0.3, 0.5])  # θ_sh<θ_th（膝曲げ）
    x_post = kneestrike_map(x_pre, P)
    assert np.allclose(x_post[:3], x_pre[:3])
    assert x_post[4] == x_post[5]  # 膝ロック
    assert kinetic_energy(x_post, P) < kinetic_energy(x_pre, P)


def test_heelstrike_swaps_and_dissipates():
    th = -0.18
    x_pre = np.array([th, -th, -th, -1.3, -0.5, -0.5])  # locked, strike 面
    x_post = heelstrike_map(x_pre, P)
    assert x_post[0] == -th and x_post[1] == th and x_post[2] == th
    assert x_post[4] == x_post[5]
    assert kinetic_energy(x_post, P) < kinetic_energy(x_pre, P)


def test_rocker_kneed_reduces_to_pointfoot_kneed():
    """R→0 で力学（両相）・衝突（knee/heel）が点足 kneed に一致。"""
    Pd = RockerKneedParams(R=1e-9, **KW)
    Pc = KneedParams(**KW)
    rng = np.random.default_rng(0)
    for _ in range(6):
        x = np.array([*rng.uniform(-0.4, 0.4, 3), *rng.uniform(-1.5, 1.5, 3)])
        du, ku = dynamics_unlocked(0.0, x, Pd), kn_unlocked(0.0, x, Pc)
        assert np.allclose(du[3:], ku[3:], rtol=0, atol=1e-4)
        xl = np.array([x[0], x[1], x[1], x[3], x[4], x[4]])  # locked: θ_th=θ_sh
        dl, kl = dynamics_locked(0.0, xl, Pd), kn_locked(0.0, xl, Pc)
        assert np.allclose(dl[3:], kl[3:], rtol=0, atol=1e-4)
        assert np.allclose(kneestrike_map(x, Pd), kn_ks(x, Pc), rtol=0, atol=1e-4)
    for _ in range(6):
        th = rng.uniform(-0.3, -0.1)
        x = np.array([th, -th, -th, *rng.uniform(-1.5, -0.2, 3)])
        assert np.allclose(heelstrike_map(x, Pd), kn_hs(x, Pc), rtol=0, atol=1e-4)
