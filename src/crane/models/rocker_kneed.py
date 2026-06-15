# src/crane/models/rocker_kneed.py
"""Rocker-foot kneed walker（McGeer 1990b 原機械）。derive レイヤーで記号導出し lambdify。

Phase 3 の点足 kneed (kneed.py) に「stance 足の転がり」だけを足したモデル。
状態 x = [θ_st, θ_th, θ_sh, θ̇_st, θ̇_th, θ̇_sh]、4相機械・質量モデル・knee-strike は
kneed.py と同一。差分は 5 箇所のみ:
  1. hip = C_st − (L−R)·down(θ_st), C_st = [−R·θ_st, R]（stance 円弧足が転がる）
  2. stance 接触点 P_st = [−R·θ_st, 0]
  3. swing 円弧足接触 P_sw = [C_sw_x, 0], C_sw = knee_sw + (l_s−R)·down(θ_sh)
  4. 衝突 pivot: knee-strike/heel-strike の系全体 pivot を [0,0]→P_st, foot_sw→P_sw
  5. パラメータ R
R→0 で kneed.py（点足 kneed）に厳密縮退。strike 面 θ_st+θ_sw=0（両足同半径の対称性）。
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import sympy as sp

from crane.derive.impact import angular_momentum
from crane.derive.lagrange import derive_qdd
from crane.model import HybridModel, PhaseSpec


@dataclass(frozen=True)
class RockerKneedParams:
    m_h: float  # hip 質量 [kg]
    m_t: float  # 大腿質量 [kg]
    m_s: float  # 脛質量 [kg]
    l_t: float  # 大腿長 hip→knee [m]
    l_s: float  # 脛長 knee→foot [m]
    b_t: float  # hip→大腿質量点 [m]
    b_s: float  # knee→脛質量点 [m]
    R: float  # 足半径 [m]
    gamma: float  # 斜面角 [rad]
    g: float = 9.81

    @property
    def l(self) -> float:  # noqa: E743
        return self.l_t + self.l_s


@lru_cache(maxsize=1)
def _build():
    th_st, th_th, th_sh = sp.symbols("th_st th_th th_sh")
    w_st, w_th, w_sh = sp.symbols("w_st w_th w_sh")
    m_h, m_t, m_s = sp.symbols("m_h m_t m_s", positive=True)
    l_t, l_s, b_t, b_s, R = sp.symbols("l_t l_s b_t b_s R", positive=True)
    gamma, g = sp.symbols("gamma g", positive=True)
    length = l_t + l_s

    def down(theta):
        return sp.Matrix([sp.sin(theta), -sp.cos(theta)])

    # (1) stance 脚: 円弧足が転がる。曲率中心 C_st 高さ R、hip は C_st から脚上向き (L−R)
    C_st = sp.Matrix([-R * th_st, R])
    hip = C_st - (length - R) * down(th_st)
    knee_st = hip + l_t * down(th_st)
    p_th_st = hip + b_t * down(th_st)
    p_sh_st = knee_st + b_s * down(th_st)
    # swing 脚
    p_th_sw = hip + b_t * down(th_th)
    knee_sw = hip + l_t * down(th_th)
    p_sh_sw = knee_sw + b_s * down(th_sh)
    # (2)(3) 転がり接触点
    P_st = sp.Matrix([-R * th_st, 0])
    C_sw = knee_sw + (l_s - R) * down(th_sh)  # swing 円弧足の曲率中心
    P_sw = sp.Matrix([C_sw[0], 0])

    q3, qd3 = [th_st, th_th, th_sh], [w_st, w_th, w_sh]
    bodies = [(m_h, hip), (m_t, p_th_st), (m_s, p_sh_st), (m_t, p_th_sw), (m_s, p_sh_sw)]

    def vel(r, q, qd):
        return r.jacobian(q) * sp.Matrix(qd)

    T3 = sum(mass * vel(r, q3, qd3).dot(vel(r, q3, qd3)) for mass, r in bodies) / 2
    g_dir = sp.Matrix([sp.sin(gamma), -sp.cos(gamma)])
    V = -g * sum(mass * r.dot(g_dir) for mass, r in bodies)

    params = (m_h, m_t, m_s, l_t, l_s, b_t, b_s, R, gamma, g)
    args6 = (th_st, th_th, th_sh, w_st, w_th, w_sh, *params)

    qdd3 = derive_qdd(q3, qd3, T3, V, simplify=False)
    f_qdd3 = sp.lambdify(args6, [qdd3[0], qdd3[1], qdd3[2]], "numpy")
    f_energy = sp.lambdify(args6, T3 + V, "numpy")
    f_kinetic = sp.lambdify(args6, T3, "numpy")

    th_sw, w_sw = sp.symbols("th_sw w_sw")
    lock = {th_th: th_sw, th_sh: th_sw, w_th: w_sw, w_sh: w_sw}
    T2 = T3.subs(lock, simultaneous=True)
    V2 = V.subs(lock, simultaneous=True)
    qdd2 = derive_qdd([th_st, th_sw], [w_st, w_sw], T2, V2, simplify=False)
    f_qdd2 = sp.lambdify((th_st, th_sw, w_st, w_sw, *params), [qdd2[0], qdd2[1]], "numpy")

    # --- knee-strike（脚交換なし、3→2 速度）: (4) 系全体 pivot を P_st に ---
    L_sys_pivot3 = angular_momentum(bodies, q3, qd3, P_st)
    L_swleg_hip3 = angular_momentum([(m_t, p_th_sw), (m_s, p_sh_sw)], q3, qd3, hip)
    wp_st, wp_sw = sp.symbols("wp_st wp_sw")
    post_ks = {w_st: wp_st, w_th: wp_sw, w_sh: wp_sw}
    eqs_ks = [
        L_sys_pivot3.subs(post_ks, simultaneous=True) - L_sys_pivot3,
        L_swleg_hip3.subs(post_ks, simultaneous=True) - L_swleg_hip3,
    ]
    A_ks, rhs_ks = sp.linear_eq_to_matrix(eqs_ks, [wp_st, wp_sw])
    qd_post_ks = A_ks.LUsolve(rhs_ks)
    f_impact_ks = sp.lambdify(args6, [qd_post_ks[0], qd_post_ks[1]], "numpy")

    # --- heel-strike（脚交換、両脚剛体、2→2 速度）: (4) pivot を P_st / P_sw に ---
    wq_st, wq_leg = sp.symbols("wq_st wq_leg")
    swap_hs = {
        th_st: th_th,
        th_th: th_st,
        th_sh: th_st,
        w_st: wq_st,
        w_th: wq_leg,
        w_sh: wq_leg,
    }
    L_stleg_hip3 = angular_momentum([(m_t, p_th_st), (m_s, p_sh_st)], q3, qd3, hip)
    L_sys_swfoot3 = angular_momentum(bodies, q3, qd3, P_sw)
    eqs_hs = [
        L_sys_pivot3.subs(swap_hs, simultaneous=True) - L_sys_swfoot3,
        L_swleg_hip3.subs(swap_hs, simultaneous=True) - L_stleg_hip3,
    ]
    A_hs, rhs_hs = sp.linear_eq_to_matrix(eqs_hs, [wq_st, wq_leg])
    qd_post_hs = A_hs.LUsolve(rhs_hs)
    f_impact_hs = sp.lambdify(args6, [qd_post_hs[0], qd_post_hs[1]], "numpy")

    return f_qdd3, f_qdd2, f_energy, f_kinetic, f_impact_ks, f_impact_hs


def _args(x, p: RockerKneedParams):
    return (*x, p.m_h, p.m_t, p.m_s, p.l_t, p.l_s, p.b_t, p.b_s, p.R, p.gamma, p.g)


def dynamics_unlocked(t: float, x, p: RockerKneedParams):
    f_qdd3 = _build()[0]
    qdd = f_qdd3(*_args(x, p))
    return [x[3], x[4], x[5], qdd[0], qdd[1], qdd[2]]


def dynamics_locked(t: float, x, p: RockerKneedParams):
    f_qdd2 = _build()[1]
    qdd = f_qdd2(
        x[0], x[1], x[3], x[4], p.m_h, p.m_t, p.m_s, p.l_t, p.l_s, p.b_t, p.b_s, p.R, p.gamma, p.g
    )
    return [x[3], x[4], x[5], qdd[0], qdd[1], qdd[1]]


def energy(x, p: RockerKneedParams) -> float:
    return float(_build()[2](*_args(x, p)))


def kinetic_energy(x, p: RockerKneedParams) -> float:
    return float(_build()[3](*_args(x, p)))


def kneestrike_map(x: np.ndarray, p: RockerKneedParams) -> np.ndarray:
    """knee-strike: 脚交換なし。swing 膝ロックで 3→2 速度。"""
    wp = _build()[4](*_args(x, p))
    return np.array([x[0], x[1], x[2], float(wp[0]), float(wp[1]), float(wp[1])])


def heelstrike_map(x: np.ndarray, p: RockerKneedParams) -> np.ndarray:
    """heel-strike: 脚交換。衝突中は両脚剛体。post で θ̇_th⁺=θ̇_sh⁺。"""
    wq = _build()[5](*_args(x, p))
    return np.array([x[1], x[0], x[0], float(wq[0]), float(wq[1]), float(wq[1])])


def make_rocker_kneed(p: RockerKneedParams) -> HybridModel:
    """相機械: unlocked → knee-strike → locked → heel-strike（kneed と同一構造）。"""
    unlocked = PhaseSpec(
        dynamics=lambda t, x: dynamics_unlocked(t, x, p),
        event_value=lambda x: x[2] - x[1],  # θ_sh − θ_th = 0 で knee-strike
        event_accept=lambda x: True,
        impact=lambda x: kneestrike_map(x, p),
    )
    locked = PhaseSpec(
        dynamics=lambda t, x: dynamics_locked(t, x, p),
        event_value=lambda x: x[0] + x[1],  # θ_st + θ_sw = 0 で heel-strike
        event_accept=lambda x: x[0] < 0.0 and (x[3] + x[4]) < 0.0,
        impact=lambda x: heelstrike_map(x, p),
    )
    return HybridModel(
        phases=(unlocked, locked),
        lift=lambda y: np.array([y[0], -y[0], -y[0], y[1], y[2], y[2]]),
        project=lambda x: np.array([x[0], x[3], x[4]]),
    )
