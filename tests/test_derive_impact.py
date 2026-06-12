import numpy as np
import sympy as sp

from crane.derive.impact import angular_momentum, cross2d


def test_cross2d():
    a = sp.Matrix([1, 0])
    b = sp.Matrix([0, 2])
    assert cross2d(a, b) == 2


def test_angular_momentum_point_mass():
    """質点の pivot 回り角運動量 L = m l² θ̇。"""
    th, w, m, leg_len = sp.symbols("th w m l", positive=True)
    r = leg_len * sp.Matrix([sp.sin(th), sp.cos(th)])
    L = angular_momentum([(m, r)], [th], [w], sp.Matrix([0, 0]))
    assert sp.simplify(L - (-m * leg_len**2 * w)) == 0 or sp.simplify(L - m * leg_len**2 * w) == 0
    # 符号は回転規約依存。絶対値が m l² w であることを数値で確認:
    f = sp.lambdify((th, w, m, leg_len), L, "numpy")
    assert np.isclose(abs(f(0.3, 1.0, 2.0, 1.5)), 2.0 * 1.5**2)


def test_pivot_transfer_analytic():
    """pivot 乗り移り衝突の解析解 ψ̇ = −0.5 を再現。"""
    th, w = sp.symbols("th w")
    m_val = 2.0
    r = sp.Matrix([sp.sin(th), sp.cos(th)])  # l=1
    P = sp.Matrix([1, 0])
    L_pre = angular_momentum([(m_val, r)], [th], [w], P)
    # post: P 回りの回転、角度 psi、半径 rho
    psi, wp = sp.symbols("psi wp")
    rho = sp.sqrt(2)
    r_post = P + rho * sp.Matrix([-sp.sin(psi), sp.cos(psi)])
    L_post = angular_momentum([(m_val, r_post)], [psi], [wp], P)
    # 保存則を wp について解く（trigsimp で sin²+cos²=1 を展開してから解く）
    sol = sp.solve(sp.Eq(sp.trigsimp(L_post), sp.trigsimp(L_pre)), wp)[0]
    f = sp.lambdify((th, w, psi), sol, "numpy")
    # th=0 (r=(0,1)), w=1; post 配置は同一点: P + sqrt(2)(−sin ψ, cos ψ) = (0,1) → ψ=π/4
    assert np.isclose(f(0.0, 1.0, np.pi / 4), -0.5)
