# src/crane/derive/impact.py
"""点質量系の角運動量と、保存則に基づく plastic impact 写像の導出部品。"""

import sympy as sp


def cross2d(a: sp.Matrix, b: sp.Matrix) -> sp.Expr:
    """2D 外積（スカラー）。"""
    return a[0] * b[1] - a[1] * b[0]


def angular_momentum(
    bodies: list[tuple[sp.Expr, sp.Matrix]],
    q: list[sp.Symbol],
    qd: list[sp.Symbol],
    about: sp.Matrix,
) -> sp.Expr:
    """点質量系の点 about 回りの角運動量。bodies = [(mass, r(q))]。

    v_i = (∂r_i/∂q)·q̇ で速度を導出する。
    """
    total = sp.S.Zero
    for mass, r in bodies:
        v = r.jacobian(q) * sp.Matrix(qd)
        total += mass * cross2d(r - about, v)
    return sp.expand(total)
