# tests/test_derive_lagrange.py
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

from crane.derive.impact import angular_momentum
from crane.derive.lagrange import derive_qdd


def test_pendulum_eom():
    """単振子: θ̈ = −(g/l) sin θ を解析解として照合。"""
    th, w, m, l, g = sp.symbols("th w m l g", positive=True)  # noqa: E741
    T = m * l**2 * w**2 / 2
    V = -m * g * l * sp.cos(th)
    qdd = derive_qdd([th], [w], T, V)
    expected = -(g / l) * sp.sin(th)
    assert sp.simplify(qdd[0] - expected) == 0


def test_free_particle_2d():
    """自由粒子: ẍ=0, ÿ=−g。"""
    x, y, vx, vy, m, g = sp.symbols("x y vx vy m g", positive=True)
    T = m * (vx**2 + vy**2) / 2
    V = m * g * y
    qdd = derive_qdd([x, y], [vx, vy], T, V)
    assert sp.simplify(qdd[0]) == 0
    assert sp.simplify(qdd[1] + g) == 0


def test_lambdify_numeric():
    """lambdify した EOM が数値評価できる。"""
    th, w, g, l = sp.symbols("th w g l", positive=True)  # noqa: E741
    T = l**2 * w**2 / 2
    V = -g * l * sp.cos(th)
    qdd = derive_qdd([th], [w], T, V)
    f = sp.lambdify((th, w, g, l), qdd[0], "numpy")
    assert np.isclose(f(0.3, 0.0, 1.0, 1.0), -np.sin(0.3))


def test_double_pendulum_matches_textbook():
    """連立2自由度: 平面二重振り子の教科書 EOM と数値一致（結合系の検証）。

    y 上向き正フレームで V = m g y。EOM を derive_qdd で導出し、
    教科書の M(q)q̈ = f(q, q̇) を独立実装して 3 点で数値比較する。
    """
    th1, th2, w1, w2 = sp.symbols("th1 th2 w1 w2")

    # 数値パラメータ（固定）
    m1_v, m2_v, l1_v, l2_v, g_v = 1.0, 1.5, 1.0, 0.8, 9.81

    # 位置 (y 上向き正): p1 = l1(sin th1, -cos th1), p2 = p1 + l2(sin th2, -cos th2)
    p1 = sp.Matrix([l1_v * sp.sin(th1), -l1_v * sp.cos(th1)])
    p2 = p1 + sp.Matrix([l2_v * sp.sin(th2), -l2_v * sp.cos(th2)])

    # 速度
    v1 = p1.jacobian([th1, th2]) * sp.Matrix([w1, w2])
    v2 = p2.jacobian([th1, th2]) * sp.Matrix([w1, w2])

    T_sym = sp.Rational(1, 2) * (m1_v * v1.dot(v1) + m2_v * v2.dot(v2))
    V_sym = m1_v * g_v * p1[1] + m2_v * g_v * p2[1]

    qdd_sym = derive_qdd([th1, th2], [w1, w2], T_sym, V_sym)
    f_qdd = sp.lambdify((th1, th2, w1, w2), qdd_sym, "numpy")

    def textbook_qdd(th1_v: float, th2_v: float, w1_v: float, w2_v: float) -> np.ndarray:
        """教科書の閉形式 M q̈ = f を numpy で解く。独立実装。"""
        d = th1_v - th2_v
        M11 = (m1_v + m2_v) * l1_v**2
        M12 = m2_v * l1_v * l2_v * np.cos(d)
        M22 = m2_v * l2_v**2
        f1 = -m2_v * l1_v * l2_v * np.sin(d) * w2_v**2 - (m1_v + m2_v) * g_v * l1_v * np.sin(th1_v)
        f2 = m2_v * l1_v * l2_v * np.sin(d) * w1_v**2 - m2_v * g_v * l2_v * np.sin(th2_v)
        M = np.array([[M11, M12], [M12, M22]])
        rhs = np.array([f1, f2])
        return np.linalg.solve(M, rhs)

    rng = np.random.default_rng(42)
    for _ in range(3):
        state = rng.uniform([-0.5, -0.5, -1.0, -1.0], [0.5, 0.5, 1.0, 1.0])
        th1_v, th2_v, w1_v, w2_v = state
        got = np.array(f_qdd(th1_v, th2_v, w1_v, w2_v)).ravel()
        expected = textbook_qdd(th1_v, th2_v, w1_v, w2_v)
        np.testing.assert_allclose(
            got,
            expected,
            atol=1e-12,
            err_msg=f"state={state}: derive_qdd={got} textbook={expected}",
        )


def test_angular_momentum_conserved_without_gravity():
    """Noether: 無重力の二重振り子で pivot 回り角運動量が積分軌道上で保存。

    lagrange.py と impact.py の相互検証。V=0 にすると方位角の循環座標から
    pivot 回り角運動量が保存量になる。solve_ivp で軌道を積分し最大ドリフトを確認。
    """
    th1, th2, w1, w2 = sp.symbols("th1 th2 w1 w2")

    l1_v, l2_v = 1.0, 0.8
    m1_v, m2_v = 1.0, 1.5

    p1 = sp.Matrix([l1_v * sp.sin(th1), -l1_v * sp.cos(th1)])
    p2 = p1 + sp.Matrix([l2_v * sp.sin(th2), -l2_v * sp.cos(th2)])

    v1 = p1.jacobian([th1, th2]) * sp.Matrix([w1, w2])
    v2 = p2.jacobian([th1, th2]) * sp.Matrix([w1, w2])

    T_sym = sp.Rational(1, 2) * (m1_v * v1.dot(v1) + m2_v * v2.dot(v2))
    V_sym = sp.Integer(0)  # 無重力

    qdd_sym = derive_qdd([th1, th2], [w1, w2], T_sym, V_sym)
    f_ode = sp.lambdify((th1, th2, w1, w2), qdd_sym, "numpy")

    origin = sp.Matrix([0, 0])
    L_sym = angular_momentum(
        [(m1_v, p1), (m2_v, p2)],
        [th1, th2],
        [w1, w2],
        origin,
    )
    f_L = sp.lambdify((th1, th2, w1, w2), L_sym, "numpy")

    def ode(t: float, y: np.ndarray) -> list[float]:
        th1v, th2v, w1v, w2v = y
        qdd = np.array(f_ode(th1v, th2v, w1v, w2v)).ravel()
        return [w1v, w2v, float(qdd[0]), float(qdd[1])]

    y0 = [0.3, -0.2, 0.5, -0.4]
    sol = solve_ivp(ode, [0.0, 1.0], y0, rtol=1e-12, atol=1e-14, dense_output=True)

    t_eval = np.linspace(0.0, 1.0, 500)
    ys = sol.sol(t_eval)  # shape (4, N)

    L_vals = np.array(
        [float(f_L(ys[0, i], ys[1, i], ys[2, i], ys[3, i])) for i in range(t_eval.size)]
    )
    drift = np.max(np.abs(L_vals - L_vals[0]))
    assert drift < 1e-9, f"Angular momentum drift = {drift:.3e} > 1e-9"
