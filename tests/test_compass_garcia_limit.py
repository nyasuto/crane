import numpy as np

from crane import references as ref
from crane.models import simplest
from crane.models.compass import CompassParams, make_compass
from crane.search import find_limit_cycle

P_DEG = CompassParams(m=1e-9, m_h=1.0, a=0.0, b=1.0, gamma=ref.GAMMA_REF, g=1.0)
MODEL = make_compass(P_DEG)
P_GARCIA = simplest.SimplestParams(gamma=ref.GAMMA_REF)


def to_garcia(x_abs):
    """compass 絶対角 → Garcia (θ, φ) 座標。"""
    return np.array([x_abs[0], x_abs[0] - x_abs[1], x_abs[2], x_abs[2] - x_abs[3]])


def to_abs(x_g):
    return np.array([x_g[0], x_g[0] - x_g[1], x_g[2], x_g[2] - x_g[3]])


def test_dynamics_reduces_to_garcia():
    """m→0, a=0 で連続相が Garcia EOM に一致（ランダム20状態）。"""
    rng = np.random.default_rng(0)
    max_err = 0.0
    for _ in range(20):
        x_g = rng.uniform([-0.3, -0.6, -2.0, -2.0], [0.3, 0.6, 2.0, 2.0])
        xdot_g = np.asarray(simplest.dynamics(0.0, x_g, P_GARCIA))
        xdot_abs = np.asarray(MODEL.dynamics(0.0, to_abs(x_g)))
        # 加速度を Garcia 座標に変換: θ̈=ẅ_st, φ̈=ẅ_st−ẅ_sw
        got = np.array([xdot_abs[2], xdot_abs[2] - xdot_abs[3]])
        err = np.max(np.abs(got - xdot_g[2:]))
        max_err = max(max_err, err)
        assert np.allclose(got, xdot_g[2:], rtol=0, atol=1e-5), (
            f"dynamics mismatch at x_g={x_g}: got={got}, expected={xdot_g[2:]}, err={err:.2e}"
        )
    print(f"\n[test_dynamics_reduces_to_garcia] max acceleration error = {max_err:.2e}")


def test_impact_reduces_to_garcia():
    """m→0, a=0 で衝突写像が Garcia の cos(2θ) 写像に一致。"""
    rng = np.random.default_rng(1)
    max_err = 0.0
    for _ in range(10):
        theta = rng.uniform(-0.3, -0.1)
        theta_dot = rng.uniform(-0.5, -0.1)
        phi_dot = rng.uniform(-0.5, 0.5)
        x_g_pre = np.array([theta, 2 * theta, theta_dot, phi_dot])
        expected = simplest.heelstrike_map(x_g_pre)
        got = to_garcia(MODEL.impact(to_abs(x_g_pre)))
        err = np.max(np.abs(got - expected))
        max_err = max(max_err, err)
        assert np.allclose(got, expected, rtol=0, atol=1e-5), (
            f"impact mismatch at theta={theta:.3f}: got={got}, expected={expected}, err={err:.2e}"
        )
    print(f"\n[test_impact_reduces_to_garcia] max impact error = {max_err:.2e}")


def test_fixed_point_reduces_to_garcia():
    """退化 compass の不動点・固有値が Phase 1 の結果に一致。"""
    x_g = simplest.lift(ref.LONG_PERIOD_THETA, ref.LONG_PERIOD_THETA_DOT)
    x_abs = to_abs(x_g)
    y_guess = np.array([x_abs[0], x_abs[2], x_abs[3]])
    fp = find_limit_cycle(MODEL, y_guess)
    assert fp.converged
    # Phase 1 実測の真の不動点 (0.2003109, -0.1998325) と照合
    err_theta = abs(fp.y[0] - 0.2003109)
    err_thetadot = abs(fp.y[1] - (-0.1998325))
    print(f"\n[test_fixed_point_reduces_to_garcia] fp.y = {fp.y}")
    print(f"  err_theta = {err_theta:.2e}, err_thetadot = {err_thetadot:.2e}")
    assert np.isclose(fp.y[0], 0.2003109, atol=1e-5), (
        f"theta* mismatch: {fp.y[0]:.7f} vs 0.2003109, err={err_theta:.2e}"
    )
    assert np.isclose(fp.y[1], -0.1998325, atol=1e-5), (
        f"theta_dot* mismatch: {fp.y[1]:.7f} vs -0.1998325, err={err_thetadot:.2e}"
    )
    mags = np.sort(np.abs(fp.eigenvalues))
    print(f"  |eigenvalues| (sorted) = {mags}")
    # 3D 断面の multiplier: 上位2つは Garcia の複素対 |λ|=0.5891、
    # 第3方向は m→0 で swing 速度が slave されるため ≈ 0
    err_eig = np.max(np.abs(mags[1:] - 0.5891))
    print(f"  err_eigenvalues = {err_eig:.2e}")
    assert np.allclose(mags[1:], [0.5891, 0.5891], atol=2e-3), (
        f"eigenvalue mismatch: {mags[1:]} vs [0.5891, 0.5891], err={err_eig:.2e}"
    )
    assert mags[0] < 1e-3, f"smallest eigenvalue too large: {mags[0]:.2e}"
