# tests/test_derive_lagrange.py
import numpy as np
import sympy as sp

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
