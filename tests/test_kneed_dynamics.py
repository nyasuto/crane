# tests/test_kneed_dynamics.py
import numpy as np
from scipy.integrate import solve_ivp

from crane import references_kneed as ref
from crane.models.kneed import KneedParams, dynamics_locked, dynamics_unlocked, energy

P = KneedParams(
    m_h=ref.M_HIP,
    m_t=ref.M_THIGH,
    m_s=ref.M_SHANK,
    l_t=ref.L_THIGH,
    l_s=ref.L_SHANK,
    b_t=ref.B_THIGH,
    b_s=ref.B_SHANK,
    gamma=ref.GAMMA_GAIT,
    g=ref.G,
)


def test_equilibrium_aligned_with_gravity():
    """全リンクが重力方向 (絶対角 γ) に揃った静止状態では加速度ゼロ（両相）。"""
    x = np.array([P.gamma, P.gamma, P.gamma, 0.0, 0.0, 0.0])
    assert np.allclose(dynamics_unlocked(0.0, x, P), np.zeros(6), atol=1e-11)
    assert np.allclose(dynamics_locked(0.0, x, P), np.zeros(6), atol=1e-11)


def test_unlocked_phase_conserves_energy():
    """unlocked 相は保存系。"""
    x0 = np.array([0.2, -0.2, -0.45, -1.2, 1.0, 2.5])
    sol = solve_ivp(lambda t, x: dynamics_unlocked(t, x, P), (0.0, 0.4), x0, rtol=1e-11, atol=1e-13)
    e0 = energy(sol.y[:, 0], P)
    drift = max(abs(energy(sol.y[:, k], P) - e0) for k in range(sol.y.shape[1]))
    assert drift < 1e-7 * abs(e0)


def test_locked_phase_conserves_energy_and_keeps_alignment():
    """locked 相は保存系で、θ_th = θ_sh が維持される（6D 埋め込みの整合）。"""
    x0 = np.array([0.25, -0.3, -0.3, -1.5, 0.8, 0.8])
    sol = solve_ivp(lambda t, x: dynamics_locked(t, x, P), (0.0, 0.4), x0, rtol=1e-11, atol=1e-13)
    e0 = energy(sol.y[:, 0], P)
    drift = max(abs(energy(sol.y[:, k], P) - e0) for k in range(sol.y.shape[1]))
    assert drift < 1e-7 * abs(e0)
    assert np.allclose(sol.y[1], sol.y[2], atol=1e-9)  # θ_th ≡ θ_sh
    assert np.allclose(sol.y[4], sol.y[5], atol=1e-9)  # θ̇_th ≡ θ̇_sh


def test_locked_dynamics_is_unlocked_restricted():
    """locked 力学は「θ_th=θ_sh 拘束付き unlocked」と stance 加速度が異なる
    （拘束力が働くため一般に一致しない）が、m_s→0 では一致するはず。
    ここでは整合性のみ: locked の θ̈_th と θ̈_sh が厳密に等しいこと。"""
    x = np.array([0.2, -0.35, -0.35, -1.0, 0.5, 0.5])
    xdot = dynamics_locked(0.0, x, P)
    assert xdot[4] == xdot[5]
