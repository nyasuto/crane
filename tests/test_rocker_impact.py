import numpy as np

from crane.models.compass import CompassParams, heelstrike_map as compass_strike
from crane.models.compass import dynamics as compass_dyn
from crane.models.rocker_compass import (
    RockerCompassParams, dynamics, heelstrike_map, kinetic_energy,
)

P = RockerCompassParams(m=5.0, m_h=10.0, c=0.5, rho=0.15, L=1.0, R=0.3,
                        gamma=0.05, g=9.81)


def test_heelstrike_swaps_and_dissipates():
    th = -0.18
    x_pre = np.array([th, -th, -1.4, -0.5])   # strike 面 θ_st+θ_sw=0 上
    x_post = heelstrike_map(x_pre, P)
    assert x_post[0] == -th          # 脚交換
    assert x_post[1] == th
    assert kinetic_energy(x_post, P) < kinetic_energy(x_pre, P)


def test_heelstrike_zero_velocity_is_fixed():
    th = -0.18
    x_pre = np.array([th, -th, 0.0, 0.0])
    x_post = heelstrike_map(x_pre, P)
    assert np.allclose(x_post[2:], np.zeros(2), atol=1e-12)


def test_rocker_reduces_to_compass_dynamics_and_impact():
    """R→0, ρ→0, c=b で rocker が点足 compass の力学・衝突に一致。"""
    Pd = RockerCompassParams(m=5.0, m_h=10.0, c=0.5, rho=1e-9, L=1.0, R=1e-9,
                             gamma=0.05, g=9.81)
    Pc = CompassParams(m=5.0, m_h=10.0, a=0.5, b=0.5, gamma=0.05, g=9.81)
    rng = np.random.default_rng(0)
    for _ in range(8):
        th_st, th_sw = rng.uniform(-0.4, 0.4, 2)
        w_st, w_sw = rng.uniform(-2.0, 2.0, 2)
        x = np.array([th_st, th_sw, w_st, w_sw])
        dr = dynamics(0.0, x, Pd)
        dc = compass_dyn(0.0, x, Pc)
        assert np.allclose([dr[2], dr[3]], [dc[2], dc[3]], rtol=0, atol=1e-4)
    for _ in range(8):
        th = rng.uniform(-0.3, -0.1)
        w_st, w_sw = rng.uniform(-2.0, -0.2, 2)
        x = np.array([th, -th, w_st, w_sw])
        gr = heelstrike_map(x, Pd)
        gc = compass_strike(x, Pc)
        assert np.allclose(gr, gc, rtol=0, atol=1e-4)
