# src/crane/models/rocker_compass.py
"""Rocker-foot（円弧足）2D compass。derive レイヤーで記号導出し lambdify。

状態 x = [θ_st, θ_sw, θ̇_st, θ̇_sw]（slope 法線基準の絶対角、進行 +x）。
各脚 = 一般剛体を「hip から距離 c±ρ の 2 質点 (m/2)」で等価表現（総質量 m・
CoM 位置 c・慣性 mρ² を再現。ρ→0 で質点に縮退）。stance 円弧足は曲率中心
C=[−R·θ_st, R] で転がる。hip = C − (L−R)·down(θ_st)。R→0 で点足 compass に縮退。

strike 面 g(x)=θ_st+θ_sw=0（両足同半径 R の幾何対称性から compass と同一）。
衝突は新接触点 P_sw 回りの角運動量保存 + 脚交換（compass と同じ swap 構造）。
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import sympy as sp

from crane.derive.impact import angular_momentum
from crane.derive.lagrange import derive_qdd


@dataclass(frozen=True)
class RockerCompassParams:
    m: float  # 脚質量 [kg]
    m_h: float  # hip 質量 [kg]
    c: float  # hip → 脚 CoM 距離 [m]
    rho: float  # 脚 CoM 回りの回転半径 [m]
    L: float  # 脚長 hip→直立時接地点 [m]
    R: float  # 足半径 [m]
    gamma: float  # 斜面角 [rad]
    g: float = 9.81


@lru_cache(maxsize=1)
def _build():
    th_st, th_sw, w_st, w_sw = sp.symbols("th_st th_sw w_st w_sw")
    m, m_h, c, rho, L, R, gamma, g = sp.symbols("m m_h c rho L R gamma g", positive=True)
    q, qd = [th_st, th_sw], [w_st, w_sw]

    def down(theta):
        return sp.Matrix([sp.sin(theta), -sp.cos(theta)])

    C_st = sp.Matrix([-R * th_st, R])  # stance 円弧足の曲率中心（高さ R、転がりで C_x=-R*th_st）
    hip = C_st - (L - R) * down(th_st)  # R→0 で [-L sinθ, L cosθ] = 点足 compass の hip
    # 各脚 = hip から距離 c±ρ の 2 質点（m/2 ずつ）。総質量 m・CoM 位置 c・慣性 mρ² を再現
    p_st1 = hip + (c - rho) * down(th_st)
    p_st2 = hip + (c + rho) * down(th_st)
    p_sw1 = hip + (c - rho) * down(th_sw)
    p_sw2 = hip + (c + rho) * down(th_sw)
    C_sw = hip + (L - R) * down(th_sw)  # swing 円弧足の曲率中心
    foot_sw_contact = sp.Matrix([C_sw[0], 0])  # strike 時の swing 接触点（衝突 pivot）
    P_st = sp.Matrix([-R * th_st, 0])  # stance 接触点

    bodies = [
        (m_h, hip),
        (m / 2, p_st1),
        (m / 2, p_st2),
        (m / 2, p_sw1),
        (m / 2, p_sw2),
    ]

    def vel(r):
        return r.jacobian(q) * sp.Matrix(qd)

    T = sum(mass * vel(r).dot(vel(r)) for mass, r in bodies) / 2
    g_dir = sp.Matrix([sp.sin(gamma), -sp.cos(gamma)])
    V = -g * sum(mass * r.dot(g_dir) for mass, r in bodies)

    qdd = derive_qdd(q, qd, T, V, simplify=False)
    params = (m, m_h, c, rho, L, R, gamma, g)
    args = (th_st, th_sw, w_st, w_sw, *params)
    f_qdd = sp.lambdify(args, [qdd[0], qdd[1]], "numpy")
    f_energy = sp.lambdify(args, T + V, "numpy")
    f_kinetic = sp.lambdify(args, T, "numpy")

    st_masses = [(m / 2, p_st1), (m / 2, p_st2)]
    sw_masses = [(m / 2, p_sw1), (m / 2, p_sw2)]
    L_sys_stance = angular_momentum(bodies, q, qd, P_st)
    L_sys_swfoot = angular_momentum(bodies, q, qd, foot_sw_contact)
    L_st_hip = angular_momentum(st_masses, q, qd, hip)
    L_sw_hip = angular_momentum(sw_masses, q, qd, hip)

    wp_st, wp_sw = sp.symbols("wp_st wp_sw")
    swap = {th_st: th_sw, th_sw: th_st, w_st: wp_st, w_sw: wp_sw}
    eqs = [
        L_sys_stance.subs(swap, simultaneous=True) - L_sys_swfoot,
        L_sw_hip.subs(swap, simultaneous=True) - L_st_hip,
    ]
    A, rhs = sp.linear_eq_to_matrix(eqs, [wp_st, wp_sw])
    qd_post = A.LUsolve(rhs)
    f_impact = sp.lambdify(args, [qd_post[0], qd_post[1]], "numpy")

    return f_qdd, f_energy, f_kinetic, f_impact


def _args(x, p: RockerCompassParams):
    return (*x, p.m, p.m_h, p.c, p.rho, p.L, p.R, p.gamma, p.g)


def dynamics(t: float, x, p: RockerCompassParams):
    f_qdd = _build()[0]
    qdd = f_qdd(*_args(x, p))
    return [x[2], x[3], qdd[0], qdd[1]]


def energy(x, p: RockerCompassParams) -> float:
    return float(_build()[1](*_args(x, p)))


def kinetic_energy(x, p: RockerCompassParams) -> float:
    return float(_build()[2](*_args(x, p)))


def heelstrike_map(x: np.ndarray, p: RockerCompassParams) -> np.ndarray:
    """衝突写像 + 脚ラベル交換。post 速度は post ラベルで返る。"""
    wp = _build()[3](*_args(x, p))
    return np.array([x[1], x[0], float(wp[0]), float(wp[1])])
