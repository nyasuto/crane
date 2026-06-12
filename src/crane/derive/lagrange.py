# src/crane/derive/lagrange.py
"""Euler-Lagrange 自動導出（scleronomic 系、q/q̇ は plain symbols）。"""

import sympy as sp


def derive_qdd(
    q: list[sp.Symbol],
    qd: list[sp.Symbol],
    T: sp.Expr,
    V: sp.Expr,
    *,
    simplify: bool = True,
) -> sp.Matrix:
    """M(q)·q̈ = ∂L/∂q − (∂p/∂q)·q̇ を解いて q̈ の式ベクトルを返す。

    p = ∂L/∂q̇（一般化運動量）。時間陽依存なし（scleronomic）前提:
    d/dt p = (∂p/∂q)·q̇ + (∂p/∂q̇)·q̈ と展開できる。
    """
    L = T - V
    p_vec = sp.Matrix([sp.diff(L, v) for v in qd])
    M = p_vec.jacobian(qd)  # = Hessian_q̇(T)
    rhs = sp.Matrix([sp.diff(L, qi) for qi in q]) - p_vec.jacobian(q) * sp.Matrix(qd)
    qdd = M.LUsolve(rhs)
    return sp.simplify(qdd) if simplify else qdd
