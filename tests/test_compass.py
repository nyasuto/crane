# tests/test_compass.py
import numpy as np
from scipy.integrate import solve_ivp

from crane.models.compass import CompassParams, energy, kinetic_energy, make_compass

P = CompassParams(m=5.0, m_h=10.0, a=0.5, b=0.5, gamma=0.05)
MODEL = make_compass(P)


def test_equilibrium_aligned_with_gravity():
    """両脚が重力方向 (絶対角 γ) に揃った静止状態では加速度ゼロ。"""
    x = np.array([P.gamma, P.gamma, 0.0, 0.0])
    xdot = MODEL.phases[0].dynamics(0.0, x)
    assert np.allclose(xdot, [0.0, 0.0, 0.0, 0.0], atol=1e-12)


def test_swing_phase_conserves_energy():
    """連続相は保存系: E = T + V が一定（受動歩行の本質）。"""
    x0 = np.array([0.2, -0.3, -1.0, 0.5])
    sol = solve_ivp(
        MODEL.phases[0].dynamics, (0.0, 0.5), x0, rtol=1e-11, atol=1e-13, dense_output=True
    )
    e0 = energy(sol.y[:, 0], P)
    drift = max(abs(energy(sol.y[:, k], P) - e0) for k in range(sol.y.shape[1]))
    assert drift < 1e-7 * abs(e0)


def test_impact_dissipates_kinetic_energy():
    """plastic 衝突は運動エネルギーを増やさない（一般に厳密減少）。"""
    x_pre = np.array([-0.15, 0.15, -1.2, -0.8])  # strike 面上 (θ_st+θ_sw=0)、下降接地
    assert abs(MODEL.phases[0].event_value(x_pre)) < 1e-15
    assert MODEL.phases[0].event_accept(x_pre)
    x_post = MODEL.phases[0].impact(x_pre)
    assert kinetic_energy(x_post, P) < kinetic_energy(x_pre, P)


def test_impact_swaps_leg_labels():
    """衝突で配置は脚ラベル交換: θ_st⁺=θ_sw⁻, θ_sw⁺=θ_st⁻。"""
    x_pre = np.array([-0.15, 0.15, -1.2, -0.8])
    x_post = MODEL.phases[0].impact(x_pre)
    assert x_post[0] == x_pre[1]
    assert x_post[1] == x_pre[0]


def test_impact_zero_velocity_is_fixed():
    """静止状態の衝突は静止のまま（線形写像の自明解）。"""
    x_pre = np.array([-0.15, 0.15, 0.0, 0.0])
    x_post = MODEL.phases[0].impact(x_pre)
    assert np.allclose(x_post[2:], [0.0, 0.0], atol=1e-12)


def test_section_lift_project_roundtrip():
    y = np.array([0.2, -1.1, 0.3])
    x = MODEL.lift(y)
    assert np.isclose(x[1], -x[0])
    assert np.allclose(MODEL.project(x), y)
