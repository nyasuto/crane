import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from crane.models.rocker_compass import RockerCompassParams, dynamics, energy

P = RockerCompassParams(m=5.0, m_h=10.0, c=0.5, rho=0.15, L=1.0, R=0.3, gamma=0.05, g=9.81)


def test_static_equilibrium_recovers_pointfoot_as_R_to_zero():
    """転がり足の静的平衡角は γ ではないが、R→0 で点足の平衡角 γ に縮退する。

    点足 compass では「両脚を重力方向 (絶対角 γ) に揃えた静止」が固定点だが、
    転がり足では接触点 P_st=[−R·θ,0] が θ で移動するため平衡角は γ からずれる。
    平衡角は qdd_st(θ,θ,0,0)=0 の根として定義し、R→0 で γ に収束することを確認。
    """

    def equilibrium_angle(R: float) -> float:
        p = RockerCompassParams(m=5.0, m_h=10.0, c=0.5, rho=0.15, L=1.0, R=R, gamma=0.05, g=9.81)
        return brentq(
            lambda th: dynamics(0.0, np.array([th, th, 0.0, 0.0]), p)[2],
            -0.5,
            0.5,
        )

    # 有限半径では平衡角 ≠ γ（転がり足の本質）。
    assert abs(equilibrium_angle(0.3) - P.gamma) > 1e-2
    # R→0 で点足の平衡角 γ に縮退。
    assert abs(equilibrium_angle(1e-5) - P.gamma) < 1e-4


def test_rolling_phase_conserves_energy():
    """転がり相は保存系（エネルギー drift < 1e-7）。"""
    x0 = np.array([0.2, -0.25, -1.0, 0.3])
    sol = solve_ivp(
        lambda t, x: dynamics(t, x, P),
        (0.0, 0.4),
        x0,
        rtol=1e-11,
        atol=1e-13,
    )
    e0 = energy(sol.y[:, 0], P)
    drift = max(abs(energy(sol.y[:, k], P) - e0) for k in range(sol.y.shape[1]))
    assert drift < 1e-7 * abs(e0)
