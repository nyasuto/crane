import numpy as np
from scipy.integrate import solve_ivp

from crane import references_kneed as kref
from crane.models.rocker_kneed import RockerKneedParams, dynamics_unlocked, dynamics_locked, energy

P = RockerKneedParams(
    m_h=kref.M_HIP,
    m_t=kref.M_THIGH,
    m_s=kref.M_SHANK,
    l_t=kref.L_THIGH,
    l_s=kref.L_SHANK,
    b_t=kref.B_THIGH,
    b_s=kref.B_SHANK,
    R=0.2,
    gamma=kref.GAMMA_GAIT,
    g=kref.G,
)


def test_unlocked_phase_conserves_energy():
    """unlocked 相（膝曲げ可、3 DOF）は保存系（drift < 1e-7）。"""
    x0 = np.array([0.2, -0.2, -0.35, -1.0, -0.3, 0.4])
    sol = solve_ivp(lambda t, x: dynamics_unlocked(t, x, P), (0.0, 0.3), x0, rtol=1e-11, atol=1e-13)
    e0 = energy(sol.y[:, 0], P)
    drift = max(abs(energy(sol.y[:, k], P) - e0) for k in range(sol.y.shape[1]))
    assert drift < 1e-7 * abs(e0)


def test_locked_phase_conserves_energy():
    """locked 相（膝ロック、θ_th=θ_sh）も保存系。"""
    th = 0.2
    x0 = np.array([th, -0.25, -0.25, -1.0, -0.3, -0.3])
    sol = solve_ivp(lambda t, x: dynamics_locked(t, x, P), (0.0, 0.3), x0, rtol=1e-11, atol=1e-13)
    e0 = energy(sol.y[:, 0], P)
    drift = max(abs(energy(sol.y[:, k], P) - e0) for k in range(sol.y.shape[1]))
    assert drift < 1e-7 * abs(e0)
